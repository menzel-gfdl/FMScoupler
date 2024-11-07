program main
use aerosol_diagnostics, only: AerosolDiagnostics
use am4, only: Atmosphere_t, create_atmosphere, destroy_atmosphere, h2o, o3, &
               read_time_slice
use block_control_mod, only: block_control_type
use cloud_diagnostics, only: CloudDiagnostics
use constants_mod, only: pi
use diag_manager_mod, only: diag_axis_init, diag_manager_end, diag_manager_init, &
                            diag_manager_set_time_end, diag_send_complete
use field_manager_mod, only: model_atmos
use flux_diagnostics, only: FluxDiagnostics
use fms_mod, only: check_nml_error, clock_loop, error_mesg, fatal, fms_end, fms_init, input_nml_file, &
                   mpp_clock_begin, mpp_clock_end, mpp_clock_id, mpp_pe, mpp_root_pe, stdlog
use get_cal_time_mod, only: get_cal_time
use gfdl_fluxes, only: BroadbandFluxes
use grid2_mod, only: grid_end, grid_init
use mpp_domains_mod, only: domain2d, mpp_get_compute_domain
use mpp_mod, only: mpp_clock_set_grain
use physics_radiation_exch_mod, only: clouds_from_moist_block_type
use radiation_context, only: all_, clean, clean_clear, clear, RadiationContext
use random_numbers_mod, only: randomnumberstream
use random_number_streams_mod, only: random_number_streams_init, &
                                     random_number_streams_end, &
                                     get_random_number_streams
use solar_constant, only: SolarConstant
use solar_spectrum, only: SolarSpectrum
use time_interp_external2_mod, only: time_interp_external_init
use time_manager_mod, only: get_date, julian, print_time, set_calendar_type, time_manager_init, &
                            time_type, operator(+), operator(-)
use tracer_manager_mod, only: get_number_tracers, get_tracer_index, &
                              tracer_manager_end, tracer_manager_init
use utilities, only: catch_error, integrate
implicit none


type(AerosolDiagnostics), dimension(:), allocatable :: aerosol_species_diags
type(Atmosphere_t), dimension(:), allocatable, target :: atm
integer, dimension(4) :: axes
integer :: block_
type(CloudDiagnostics) :: cloud_diags
integer, dimension(6) :: date
real :: dt
type(FluxDiagnostics), dimension(4) :: flux_diags
integer :: i
integer :: id_lat
integer :: id_lon
integer :: id_pfull
integer :: id_phalf
integer :: ios
integer :: j
integer :: k
integer :: last_infrared_band !Index of the last "nir" band.
integer :: logfile_handle
integer :: longwave_axis_id
real, dimension(:, :), allocatable :: longwave_band_limits
integer :: num_bands
integer :: num_blocks
integer, dimension(:), allocatable :: num_columns
integer :: num_lat
integer :: num_layers
integer :: num_levels
integer :: num_lon
type(RadiationContext) :: radiation_context
type(time_type) :: time
type(time_type) :: time_next
type(time_type) :: timestep
integer :: shortwave_axis_id
real, dimension(:,:), allocatable :: shortwave_band_limits
type(SolarSpectrum) :: solar_flux_spectrum
type(SolarConstant) :: solar_flux_constant
real :: surface_albedo_weight !Weighting needed to combine "nir" and "vis" albedo values in
                              !in the band that contains the infrared cut-off.
integer :: t

!MPP timers.
integer :: aerosol_optics_clock
integer :: cloud_optics_clock
integer :: flux_solver_clock
integer :: gas_optics_clock
integer :: radiation_driver_clock

!Physics driver parameters.
real, parameter :: max_diam_drop = 50.e-6
real, parameter :: min_diam_ice = 10.e-6
real, parameter :: min_diam_drop = 2.e-6
real, parameter :: n_min = 1.e6
real, parameter :: dcs = 200.e-6
real, parameter :: qcvar = 1.
real, parameter :: qmin = 1.e-10

!Lon-lat column blocking.
type(block_control_type) :: column_blocking

!Runtime options.
integer :: near_infrared_cutoff = 14600 !Wavenumber [cm-1] that distinguishes the visible from the near-infrared.
integer :: nxblocks = 1
integer :: nyblocks = 1
character(len=256) :: solar_constant_path = ""
character(len=256) :: solar_spectrum_path = ""
namelist /standalone_radiation_nml/ near_infrared_cutoff, &
                                    nxblocks, &
                                    nyblocks, &
                                    solar_constant_path, &
                                    solar_spectrum_path

!Start up mpi.
call fms_init()
call mpp_clock_set_grain(clock_loop)
call time_manager_init()
call tracer_manager_init()
call grid_init()
call time_interp_external_init()
call random_number_streams_init()

!Set up timers.
radiation_driver_clock = mpp_clock_id("  Radiation: radiation driver")
aerosol_optics_clock = mpp_clock_id("    Radiation: aerosol optical properties")
cloud_optics_clock = mpp_clock_id("    Radiation: cloud optical properties")
flux_solver_clock = mpp_clock_id("    Radiation: flux solver")
gas_optics_clock = mpp_clock_id("    Radiation: gas optical properties")

!Read the namelist.
read(input_nml_file, nml=standalone_radiation_nml, iostat=ios)
ios = check_nml_error(ios, "standalone_radiation_nml")
logfile_handle = stdlog()
if (mpp_pe() .eq. mpp_root_pe()) then
  write(logfile_handle, nml=standalone_radiation_nml)
endif

!Allocate space for the input data and create the column blocking.
call create_atmosphere(atm, column_blocking, nxblocks, nyblocks)

!Set the model time.
if (trim(atm(1)%calendar) .eq. "julian") then
  call set_calendar_type(julian)
else
  call error_mesg("main", "only julian calendar supported.", fatal)
endif
if (atm(1)%num_times .gt. 1) then
  dt = atm(1)%time(2) - atm(1)%time(1)
  do t = 3, atm(1)%num_times
    !Check to see if the input dataset has a constant timestep.
    if (abs((atm(1)%time(t) - atm(1)%time(t - 1)) - dt) .gt. 1.e-10) then
      call error_mesg("main", "timestep is not constant in the input datasets.", fatal)
    endif
  enddo
else
  !If the input dataset only has one time level, we can't determine the
  !model radiatoin timestep that was used, so set it to zero.
  dt = 0.
endif

!Model diagnostics are output at the end of a timestep, yet calculated
!using the time at the beginning of a timestep.  To mimic that, subtract
!one timestep off of the first time found in the input dataset.
time = get_cal_time(max(atm(1)%time(1) - dt, 0.), atm(1)%time_units, atm(1)%calendar)
timestep = get_cal_time(atm(1)%time(1), atm(1)%time_units, atm(1)%calendar)
timestep = timestep - time
write(*, *) atm(1)%time, dt


!Read in the solar data.
call solar_flux_constant%create("solar_flux", trim(solar_constant_path))
call solar_flux_spectrum%create(trim(solar_spectrum_path))

!Get horizontal domain size.
num_blocks = size(column_blocking%ibs)
allocate(num_columns(num_blocks))
do i = 1, num_blocks
  num_lon = column_blocking%ibe(1) - column_blocking%ibs(1) + 1
  num_lat = column_blocking%jbe(1) - column_blocking%jbs(1) + 1
  num_columns(i) = num_lat*num_lon
enddo
num_layers = atm(1)%num_layers
num_levels = num_layers + 1

!Initialize diag_manager.
call get_date(time, date(1), date(2), date(3), date(4), date(5), date(6))
call diag_manager_init(time_init=date)

!Initialize the radiation object.
call radiation_context%create(num_columns, num_layers, num_blocks, solar_flux_spectrum%grid, &
                              solar_flux_spectrum%flux, atm(1)%longitude_bounds, atm(1)%latitude_bounds, &
                              max_diam_drop, min_diam_ice, min_diam_drop, n_min, dcs, qcvar, qmin)
deallocate(num_columns)

!Initialize the diagnostics.
id_lon = diag_axis_init("lon", atm(1)%longitude, "degrees_east", "X", domain2=atm(1)%domain)
id_lat = diag_axis_init("lat", atm(1)%latitude, "degrees_north", "Y", domain2=atm(1)%domain)
id_phalf = diag_axis_init("phalf", atm(1)%level, "mb", "Z")
id_pfull = diag_axis_init("pfull", atm(1)%layer, "mb", "Z")
axes = [id_lon, id_lat, id_pfull, id_phalf]
call radiation_context%longwave_gas_optics%band_limits(longwave_band_limits)
call radiation_context%shortwave_gas_optics%band_limits(shortwave_band_limits)
longwave_axis_id = diag_axis_init("longwave_band", longwave_band_limits(2, :), "cm-1", "Z")
deallocate(longwave_band_limits)
shortwave_axis_id = diag_axis_init("shortwave_band", shortwave_band_limits(2, :), "cm-1", "Z")
do i = 1, size(flux_diags)
  call flux_diags(i)%create(i, time, axes, longwave_axis_id, shortwave_axis_id)
enddo
call cloud_diags%create(time, axes(:2))
allocate(aerosol_species_diags(size(radiation_context%aerosol_optics%family)))
do i = 1, size(radiation_context%aerosol_optics%family)
  call aerosol_species_diags(i)%create(time, axes, radiation_context%aerosol_optics%family(i)%name, &
                                       radiation_context%aerosol_optics%family(i)%name)
enddo
call diag_manager_set_time_end(get_cal_time(atm(1)%time(atm(1)%num_times), atm(1)%time_units, atm(1)%calendar))

!Calculate variables need for dealing with the land spectral decomposition.
num_bands = size(shortwave_band_limits, 2)
if (real(near_infrared_cutoff) .lt. shortwave_band_limits(1, 1)) then
  !There are no infrared bands.
  last_infrared_band = 0
elseif (real(near_infrared_cutoff) .gt. shortwave_band_limits(2, num_bands)) then
  !There are no visible bands.$
  last_infrared_band = num_bands + 1
else
  do i = 1, num_bands
    if (shortwave_band_limits(1, i) .le. real(near_infrared_cutoff) .and. &
        shortwave_band_limits(2, i) .ge. real(near_infrared_cutoff)) then
      !If the lower bound of the band is greater than or equal to the
      !cut-off, then the previous band is the last near-infrared band.
      last_infrared_band = i
      exit
    endif
  enddo
endif
if (last_infrared_band .ne. 0 .and. last_infrared_band .lt. num_bands) then
  i = minloc(abs(solar_flux_spectrum%grid(:) - shortwave_band_limits(1, last_infrared_band)), 1)
  j = minloc(abs(solar_flux_spectrum%grid(:) - shortwave_band_limits(2, last_infrared_band)), 1)
  k = minloc(abs(solar_flux_spectrum%grid(:) - real(near_infrared_cutoff)), 1)
  surface_albedo_weight = integrate(solar_flux_spectrum%grid(i:k), solar_flux_spectrum%flux(i:k))/ &
                          integrate(solar_flux_spectrum%grid(i:j), solar_flux_spectrum%flux(i:j))
else
  surface_albedo_weight = 1.
endif
deallocate(shortwave_band_limits)

!Main loop.
do t = 1, atm(1)%num_times
  !Calculate the current time.
  time_next = time + timestep
  call print_time(time, "Running timestep: ")

  !Read in the atmospheric properies.
  call read_time_slice(atm, t, column_blocking)

  !Time interpolation.
  call solar_flux_constant%update(time)
  call radiation_context%update(time)

!$omp parallel do private(block_) default(shared)
  do block_ = 1, num_blocks
    call radiation_scheme(radiation_context, atm(block_), column_blocking, num_layers, block_, &
                          aerosol_optics_clock, cloud_optics_clock, flux_solver_clock, &
                          gas_optics_clock, radiation_driver_clock, h2o, o3, last_infrared_band, &
                          aerosol_species_diags, cloud_diags, flux_diags, time, time_next, solar_flux_constant, &
                          surface_albedo_weight, all_, clean, clean_clear, clear)
  enddo

  !Write out diagnostics.
  call diag_send_complete(timestep)
  time = time_next
enddo

!Clean up.
if (allocated(aerosol_species_diags)) then
  do i = 1, size(aerosol_species_diags)
    call aerosol_species_diags(i)%destroy()
  enddo
  deallocate(aerosol_species_diags)
endif
call random_number_streams_end()
call diag_manager_end(time)
call radiation_context%destroy()
do i = 1, num_blocks
  call destroy_atmosphere(atm(i))
enddo
deallocate(atm)
call tracer_manager_end()
call fms_end()


contains


subroutine radiation_scheme(radiation_context, atm, column_blocking, num_layers, block_, &
                            aerosol_optics_clock, cloud_optics_clock, flux_solver_clock, &
                            gas_optics_clock, radiation_driver_clock, h2o, o3, last_infrared_band, &
                            aerosol_species_diags, cloud_diags, flux_diags, time, time_next, solar_flux_constant, &
                            surface_albedo_weight, all_, clean, clean_clear, clear)

  type(RadiationContext), intent(inout) :: radiation_context
  type(Atmosphere_t), intent(in), target :: atm
  type(block_control_type), intent(in) :: column_blocking
  integer, intent(in) :: num_layers
  integer, intent(in) :: block_
  integer, intent(in) :: aerosol_optics_clock
  integer, intent(in) :: cloud_optics_clock
  integer, intent(in) :: flux_solver_clock
  integer, intent(in) :: gas_optics_clock
  integer, intent(in) :: radiation_driver_clock
  integer, intent(in) :: h2o
  integer, intent(in) :: o3
  integer, intent(in) :: last_infrared_band
  type(AerosolDiagnostics), dimension(:), intent(inout) :: aerosol_species_diags
  type(CloudDiagnostics), intent(inout) :: cloud_diags
  type(FluxDiagnostics), dimension(:), intent(inout) :: flux_diags
  type(time_type), intent(in) :: time
  type(time_type), intent(in) :: time_next
  type(SolarConstant), intent(in) :: solar_flux_constant
  real, intent(in) :: surface_albedo_weight
  integer, intent(in) :: all_
  integer, intent(in) :: clean
  integer, intent(in) :: clean_clear
  integer, intent(in) :: clear

  integer :: band, column, i, n, num_bands, num_columns, num_lat, num_levels, num_lon, s
  real, dimension(:, :, :), allocatable :: aerosol_relative_humidity
  real, dimension(:), allocatable :: average_cloud
  real, dimension(:), allocatable :: average_high_cloud
  real, dimension(:), allocatable :: average_low_cloud
  real, dimension(:), allocatable :: average_mid_cloud
  real, dimension(:, :, :), allocatable :: buffer
  real, dimension(:), allocatable :: cloud_top_pressure
  real, dimension(:), allocatable :: cloud_top_temperature
  real, dimension(:, :), pointer :: convective_droplet_number
  real, dimension(:, :), pointer :: convective_fraction
  real, dimension(:, :), pointer :: convective_ice_content
  real, dimension(:, :), allocatable :: convective_ice_size
  real, dimension(:, :), pointer :: convective_liquid_content
  real, dimension(:, :), allocatable :: convective_liquid_size
  real, dimension(:), pointer :: daylight_fraction
  real, dimension(:), pointer :: diffuse_albedo_nir
  real, dimension(:), pointer :: direct_albedo_nir
  real, dimension(:), pointer :: diffuse_albedo_vis
  real, dimension(:), pointer :: direct_albedo_vis
  real, dimension(:, :), allocatable :: diffuse_surface_albedo
  real, dimension(:, :), allocatable :: direct_surface_albedo
  real, dimension(:, :), allocatable :: flux_ratio
  real, dimension(:), pointer :: land_fraction
  real, dimension(:, :), pointer :: layer_pressure
  real, dimension(:, :), pointer :: layer_temperature
  real, dimension(:, :), pointer :: layer_thickness
  real, dimension(:, :), pointer :: level_pressure
  real, dimension(:, :), pointer :: level_temperature
  type(BroadbandFluxes), dimension(4) :: longwave_broadband_fluxes
  integer, dimension(:, :), allocatable :: longwave_gpoint_limits
  real, dimension(:, :), pointer :: ozone
  real, dimension(:, :), allocatable :: rain_size
  type(BroadbandFluxes), dimension(4) :: shortwave_broadband_fluxes
  integer, dimension(:, :), allocatable :: shortwave_gpoint_limits
  real, dimension(:, :), allocatable :: snow_size
  real, dimension(:, :), pointer :: stratiform_droplet_number
  real, dimension(:, :), pointer :: stratiform_fraction
  real, dimension(:, :), pointer :: stratiform_ice_content
  real, dimension(:, :), pointer :: stratiform_ice_number
  real, dimension(:, :), allocatable :: stratiform_ice_size
  real, dimension(:, :), pointer :: stratiform_liquid_content
  real, dimension(:, :), allocatable :: stratiform_liquid_size
  real, dimension(:, :), pointer :: stratiform_rain
  real, dimension(:, :), pointer :: stratiform_rain_size
  real, dimension(:, :), pointer :: stratiform_snow
  real, dimension(:, :), pointer :: stratiform_snow_size
  type(randomnumberstream), dimension(:, :), allocatable, target :: streams !(lon, lat).
  type(randomnumberstream), dimension(:), pointer :: streams_pointer !(column).
  real, dimension(:, :), allocatable :: surface_emissivity
  real, dimension(:), pointer :: surface_temperature
  real, dimension(:), allocatable :: total_solar_flux
  real, dimension(:, :), pointer :: water_vapor
  real, dimension(:), pointer :: zenith

  num_lon = column_blocking%ibe(block_) - column_blocking%ibs(block_) + 1
  num_lat = column_blocking%jbe(block_) - column_blocking%jbs(block_) + 1
  num_columns = num_lon*num_lat
  num_levels = num_layers + 1

  !Get gpoint limits
  call radiation_context%longwave_gas_optics%gpoint_limits(longwave_gpoint_limits)
  call radiation_context%shortwave_gas_optics%gpoint_limits(shortwave_gpoint_limits)

  !Allocate thread-specific arrays.
  allocate(aerosol_relative_humidity(num_lon, num_lat, num_layers))
  allocate(average_cloud(num_columns))
  allocate(average_high_cloud(num_columns))
  allocate(average_low_cloud(num_columns))
  allocate(average_mid_cloud(num_columns))
  allocate(cloud_top_pressure(num_columns))
  allocate(cloud_top_temperature(num_columns))
  allocate(convective_ice_size(num_columns, num_layers))
  allocate(convective_liquid_size(num_columns, num_layers))
  allocate(rain_size(num_columns, num_layers))
  allocate(snow_size(num_columns, num_layers))
  allocate(stratiform_ice_size(num_columns, num_layers))
  allocate(stratiform_liquid_size(num_columns, num_layers))

  !Get the random numbers streams.
  allocate(streams(num_lon, num_lat))
  call get_random_number_streams(atm%layer_temperature(:, :, 1), streams, 2)
  streams_pointer(1:num_columns) => streams(1:num_lon, 1:num_lat)

  !Copy aerosol relative_humidity for now.
  aerosol_relative_humidity(:, :, :) = atm%ppmv(:, :, :, h2o)

  do i = 1, size(longwave_broadband_fluxes)
    call longwave_broadband_fluxes(i)%create(num_columns, num_levels, .false.)
    call shortwave_broadband_fluxes(i)%create(num_columns, num_levels, i .eq. all_)
  enddo
  allocate(flux_ratio(num_lon, num_lat))
  flux_ratio(:, :) = 1.

  !Set the surface emissivity.
  num_bands = radiation_context%longwave_gas_optics%num_bands()
  allocate(surface_emissivity(num_bands, num_columns))
  surface_emissivity(:, :) = 1.

  !Start the radiation timer.
  call mpp_clock_begin(radiation_driver_clock)

  !Update gas concentrations.
  water_vapor(1:num_columns, 1:num_layers) => atm%ppmv(1:num_lon, 1:num_lat, 1:num_layers, h2o)
  ozone(1:num_columns, 1:num_layers) => atm%ppmv(1:num_lon, 1:num_lat, 1:num_layers, o3)
  call radiation_context%update_concentrations(water_vapor, ozone, block_)

  !Calculate gas optics.
  layer_pressure(1:num_columns, 1:num_layers) => atm%layer_pressure(1:num_lon, 1:num_lat, 1:num_layers)
  level_pressure(1:num_columns, 1:num_levels) => atm%level_pressure(1:num_lon, 1:num_lat, 1:num_levels)
  layer_temperature(1:num_columns, 1:num_layers) => atm%layer_temperature(1:num_lon, 1:num_lat, 1:num_layers)
  level_temperature(1:num_columns, 1:num_levels) => atm%level_temperature(1:num_lon, 1:num_lat, 1:num_levels)
  surface_temperature(1:num_columns) => atm%surface_temperature(1:num_lon, 1:num_lat)
  zenith(1:num_columns) => atm%solar_zenith_angle(1:num_lon, 1:num_lat)
  call mpp_clock_begin(gas_optics_clock)
  call radiation_context%calculate_gas_optics(layer_pressure, level_pressure, &
                                              layer_temperature, level_temperature, &
                                              surface_temperature, block_, zenith)
  call mpp_clock_end(gas_optics_clock)

  !Calculate cloud optics.
  land_fraction(1:num_columns) => atm%land_fraction(1:num_lon, 1:num_lat)
  convective_droplet_number(1:num_columns, 1:num_layers) => atm%shallow_droplet_number(1:num_lon, 1:num_lat, 1:num_layers)
  convective_fraction(1:num_columns, 1:num_layers) => atm%shallow_cloud_fraction(1:num_lon, 1:num_lat, 1:num_layers)
  convective_liquid_content(1:num_columns, 1:num_layers) => atm%shallow_cloud_liquid_content(1:num_lon, 1:num_lat, 1:num_layers)
  convective_ice_content(1:num_columns, 1:num_layers) => atm%shallow_cloud_ice_content(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_droplet_number(1:num_columns, 1:num_layers) => atm%stratiform_droplet_number(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_fraction(1:num_columns, 1:num_layers) => atm%stratiform_cloud_fraction(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_liquid_content(1:num_columns, 1:num_layers) => atm%stratiform_cloud_liquid_content(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_ice_content(1:num_columns, 1:num_layers) => atm%stratiform_cloud_ice_content(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_ice_number(1:num_columns, 1:num_layers) => atm%stratiform_ice_number(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_rain(1:num_columns, 1:num_layers) => atm%stratiform_rain(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_rain_size(1:num_columns, 1:num_layers) => atm%stratiform_rain_size(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_snow(1:num_columns, 1:num_layers) => atm%stratiform_snow(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_snow_size(1:num_columns, 1:num_layers) => atm%stratiform_snow_size(1:num_lon, 1:num_lat, 1:num_layers)
  layer_thickness(1:num_columns, 1:num_layers) => atm%layer_thickness(1:num_lon, 1:num_lat, 1:num_layers)
  daylight_fraction(1:num_columns) => atm%daylight_fraction(1:num_lon, 1:num_lat)
  call mpp_clock_begin(cloud_optics_clock)
  call radiation_context%calculate_cloud_optics(layer_pressure, level_pressure, layer_temperature, &
                                                level_temperature, land_fraction, &
                                                layer_thickness, convective_droplet_number, &
                                                convective_fraction, convective_ice_content, &
                                                convective_liquid_content, stratiform_droplet_number, &
                                                stratiform_fraction, stratiform_ice_content, stratiform_ice_number, &
                                                stratiform_liquid_content, stratiform_rain, stratiform_rain_size, &
                                                stratiform_snow, stratiform_snow_size, block_, &
                                                average_high_cloud, average_mid_cloud, average_low_cloud, &
                                                average_cloud, cloud_top_pressure, cloud_top_temperature, &
                                                convective_liquid_size, convective_ice_size, &
                                                stratiform_liquid_size, stratiform_ice_size, &
                                                rain_size, snow_size, streams_pointer)
  streams_pointer => null()
  deallocate(streams)
  call mpp_clock_end(cloud_optics_clock)
  call cloud_diags%save_data(average_high_cloud, average_mid_cloud, average_low_cloud, average_cloud, &
                             cloud_top_pressure, cloud_top_temperature, stratiform_liquid_size, stratiform_ice_size, &
                             rain_size, snow_size, num_lon, num_lat, &
                             column_blocking%ibs(block_) - column_blocking%isc + 1, &
                             column_blocking%jbs(block_) - column_blocking%jsc + 1, time_next)

  !Calculate aerosol optics.
  call mpp_clock_begin(aerosol_optics_clock)
  call radiation_context%calculate_aerosol_optics(atm%aerosols, aerosol_relative_humidity, &
                                                  atm%level_pressure, layer_thickness, time, &
                                                  column_blocking%ibs(block_) - column_blocking%isc + 1, &
                                                  column_blocking%ibe(block_) - column_blocking%isc + 1, &
                                                  column_blocking%jbs(block_) - column_blocking%jsc + 1, &
                                                  column_blocking%jbe(block_) - column_blocking%jsc + 1, &
                                                  block_)
  call mpp_clock_end(aerosol_optics_clock)
  allocate(buffer(num_lon, num_lat, num_layers))
  do i = 1, size(radiation_context%aerosol_optics%family)
    buffer = 0.
    do j = 1, size(radiation_context%aerosol_optics%family(i)%indices)
      n = radiation_context%aerosol_optics%family(i)%indices(j)
      buffer(:, :, :) = buffer(:, :, :) + atm%aerosols(1:num_lon, 1:num_lat, :, n)
    enddo
    call aerosol_species_diags(i)%save_data(buffer, layer_thickness, &
                                            column_blocking%ibs(block_) - column_blocking%isc + 1, &
                                            column_blocking%jbs(block_) - column_blocking%jsc + 1, &
                                            1, time_next)
  enddo
  deallocate(buffer)

  !Calculate the total aerosol optics.
  s = size(radiation_context%aerosol_optics%longwave_species) + 2
  radiation_context%longwave_aerosol_optical_properties(s, block_)%g = 0.
  radiation_context%longwave_aerosol_optical_properties(s, block_)%ssa = 0.
  radiation_context%longwave_aerosol_optical_properties(s, block_)%tau = 0.
  radiation_context%shortwave_aerosol_optical_properties(s, block_)%g = 0.
  radiation_context%shortwave_aerosol_optical_properties(s, block_)%ssa = 0.
  radiation_context%shortwave_aerosol_optical_properties(s, block_)%tau = 0.
  do i = 1, size(radiation_context%aerosol_optics%longwave_species) + 1
    call catch_error(radiation_context%longwave_aerosol_optical_properties(i, block_)%increment( &
      radiation_context%longwave_aerosol_optical_properties(s, block_)))
    call catch_error(radiation_context%shortwave_aerosol_optical_properties(i, block_)%increment( &
      radiation_context%shortwave_aerosol_optical_properties(s, block_)))
  enddo

  !Calculate the surface albedo.
  diffuse_albedo_nir(1:num_columns) => atm%surface_albedo_diffuse_ir(1:num_lon, 1:num_lat)
  direct_albedo_nir(1:num_columns) => atm%surface_albedo_direct_ir(1:num_lon, 1:num_lat)
  diffuse_albedo_vis(1:num_columns) => atm%surface_albedo_diffuse_uv(1:num_lon, 1:num_lat)
  direct_albedo_vis(1:num_columns) => atm%surface_albedo_direct_uv(1:num_lon, 1:num_lat)
  num_bands = radiation_context%shortwave_gas_optics%num_bands()
  allocate(diffuse_surface_albedo(num_bands, num_columns))
  allocate(direct_surface_albedo(num_bands, num_columns))
  do i = 1, num_bands
    if (i .lt. last_infrared_band) then
      diffuse_surface_albedo(i, :) = diffuse_albedo_nir(:)
      direct_surface_albedo(i, :) = direct_albedo_nir(:)
    elseif (i .gt. last_infrared_band) then
      diffuse_surface_albedo(i, :) = diffuse_albedo_vis(:)
      direct_surface_albedo(i, :) = direct_albedo_vis(:)
    else
      diffuse_surface_albedo(i, :) = surface_albedo_weight*diffuse_albedo_nir(:) + &
                                     (1 - surface_albedo_weight)*diffuse_albedo_vis(:)
      direct_surface_albedo(i, :) = surface_albedo_weight*direct_albedo_nir(:) + &
                                    (1 - surface_albedo_weight)*direct_albedo_vis(:)
    endif
  enddo

  !Calculate the incident solar flux.
  allocate(total_solar_flux(num_columns))
  total_solar_flux = sum(radiation_context%toa(block_)%flux, 2)
  do j = 1, num_lat
    do i = 1, num_lon
      column = (j - 1)*num_lon + i
      radiation_context%toa(block_)%flux(column, :) = radiation_context%toa(block_)%flux(column, :)* &
                                                      atm%earth_sun_distance_fraction* &
                                                      atm%daylight_fraction(i, j)*solar_flux_constant%flux/ &
                                                      total_solar_flux(column)
    enddo
  enddo
  deallocate(total_solar_flux)

  !Calculate the fluxes.
  call mpp_clock_begin(flux_solver_clock)
  !Clean-clear sky.
  call radiation_context%calculate_longwave_fluxes(surface_emissivity, clean_clear, block_)
  call radiation_context%longwave_fluxes(block_)%integrated(longwave_broadband_fluxes(clean_clear))
  call longwave_broadband_fluxes(clean_clear)%heating_rates(level_pressure)
  call radiation_context%calculate_shortwave_fluxes(zenith, direct_surface_albedo, &
                                                    diffuse_surface_albedo, clean_clear, block_)
  call radiation_context%shortwave_fluxes(block_)%integrated(shortwave_broadband_fluxes(clean_clear))
  call shortwave_broadband_fluxes(clean_clear)%heating_rates(level_pressure)
  call flux_diags(clean_clear)%save_data(longwave_gpoint_limits, &
                                         longwave_broadband_fluxes(clean_clear), &
                                         shortwave_gpoint_limits, &
                                         shortwave_broadband_fluxes(clean_clear), &
                                         level_pressure, num_lon, num_lat, &
                                         column_blocking%ibs(block_) - column_blocking%isc + 1, &
                                         column_blocking%jbs(block_) - column_blocking%jsc + 1, &
                                         1, time_next, flux_ratio)

  !Clean sky.
  call radiation_context%calculate_longwave_fluxes(surface_emissivity, clean, block_)
  call radiation_context%longwave_fluxes(block_)%integrated(longwave_broadband_fluxes(clean))
  call longwave_broadband_fluxes(clean)%heating_rates(level_pressure)
  call radiation_context%calculate_shortwave_fluxes(zenith, direct_surface_albedo, &
                                                    diffuse_surface_albedo, clean, block_)
  call radiation_context%shortwave_fluxes(block_)%integrated(shortwave_broadband_fluxes(clean))
  call shortwave_broadband_fluxes(clean)%heating_rates(level_pressure)
  call flux_diags(clean)%save_data(longwave_gpoint_limits, &
                                   longwave_broadband_fluxes(clean), &
                                   shortwave_gpoint_limits, &
                                   shortwave_broadband_fluxes(clean), &
                                   level_pressure, num_lon, num_lat, &
                                   column_blocking%ibs(block_) - column_blocking%isc + 1, &
                                   column_blocking%jbs(block_) - column_blocking%jsc + 1, &
                                   1, time_next, flux_ratio)

  !Clear sky.
  call radiation_context%calculate_longwave_fluxes(surface_emissivity, clear, block_)
  call radiation_context%longwave_fluxes(block_)%integrated(longwave_broadband_fluxes(clear))
  call longwave_broadband_fluxes(clear)%heating_rates(level_pressure)
  call radiation_context%calculate_shortwave_fluxes(zenith, direct_surface_albedo, &
                                                    diffuse_surface_albedo, clear, block_)
  call radiation_context%shortwave_fluxes(block_)%integrated(shortwave_broadband_fluxes(clear))
  call shortwave_broadband_fluxes(clear)%heating_rates(level_pressure)
  call flux_diags(clear)%save_data(longwave_gpoint_limits, &
                                   longwave_broadband_fluxes(clear), &
                                   shortwave_gpoint_limits, &
                                   shortwave_broadband_fluxes(clear), &
                                   level_pressure, num_lon, num_lat, &
                                   column_blocking%ibs(block_) - column_blocking%isc + 1, &
                                   column_blocking%jbs(block_) - column_blocking%jsc + 1, &
                                   1, time_next, flux_ratio)

  !All sky.
  call radiation_context%calculate_longwave_fluxes(surface_emissivity, all_, block_)
  call radiation_context%longwave_fluxes(block_)%integrated(longwave_broadband_fluxes(all_))
  call longwave_broadband_fluxes(all_)%heating_rates(level_pressure)
  call radiation_context%calculate_shortwave_fluxes(zenith, direct_surface_albedo, &
                                                    diffuse_surface_albedo, all_, block_)
  call radiation_context%shortwave_fluxes(block_)%integrated(shortwave_broadband_fluxes(all_))
  call shortwave_broadband_fluxes(all_)%heating_rates(level_pressure)
  call flux_diags(all_)%save_data(longwave_gpoint_limits, &
                                  longwave_broadband_fluxes(all_), &
                                  shortwave_gpoint_limits, &
                                  shortwave_broadband_fluxes(all_), &
                                  level_pressure, num_lon, num_lat, &
                                  column_blocking%ibs(block_) - column_blocking%isc + 1, &
                                  column_blocking%jbs(block_) - column_blocking%jsc + 1, &
                                  1, time_next, flux_ratio)
  call mpp_clock_end(flux_solver_clock)
  deallocate(diffuse_surface_albedo, direct_surface_albedo)
  deallocate(longwave_gpoint_limits, shortwave_gpoint_limits)
  deallocate(average_cloud)
  deallocate(average_high_cloud)
  deallocate(average_low_cloud)
  deallocate(average_mid_cloud)
  deallocate(cloud_top_pressure)
  deallocate(cloud_top_temperature)
  deallocate(convective_ice_size)
  deallocate(convective_liquid_size)
  deallocate(rain_size)
  deallocate(snow_size)
  deallocate(stratiform_ice_size)
  deallocate(stratiform_liquid_size)
  deallocate(aerosol_relative_humidity)
  deallocate(flux_ratio)
  deallocate(surface_emissivity)
  do i = 1, size(longwave_broadband_fluxes)
    call longwave_broadband_fluxes(i)%destroy()
    call shortwave_broadband_fluxes(i)%destroy()
  enddo
  call mpp_clock_end(radiation_driver_clock)
end subroutine radiation_scheme


end program main
