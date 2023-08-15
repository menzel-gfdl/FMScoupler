program main
use am4, only: Atmosphere_t, create_atmosphere, destroy_atmosphere, h2o, o3, &
               read_time_slice
use constants_mod, only: pi
use diag_manager_mod, only: diag_axis_init, diag_manager_end, diag_manager_init, &
                            diag_manager_set_time_end, diag_send_complete
use field_manager_mod, only: model_atmos
use flux_diagnostics, only: FluxDiagnostics
use fms_mod, only: check_nml_error, error_mesg, fatal, fms_end, fms_init, input_nml_file, &
                   mpp_pe, mpp_root_pe, stdlog
use get_cal_time_mod, only: get_cal_time
use grid2_mod, only: grid_end, grid_init
use mpp_domains_mod, only: domain2d, mpp_get_compute_domain, mpp_get_io_domain
use physics_radiation_exch_mod, only: clouds_from_moist_block_type
use radiation_context, only: all_, clean, clean_clear, clear, RadiationContext
use solar_constant, only: SolarConstant
use solar_spectrum, only: SolarSpectrum
use time_interp_external2_mod, only: time_interp_external_init
use time_manager_mod, only: get_date, julian, print_time, set_calendar_type, time_manager_init, &
                            time_type
use tracer_manager_mod, only: get_number_tracers, get_tracer_index, &
                              tracer_manager_end, tracer_manager_init
use utilities, only: integrate
implicit none


real, dimension(:,:,:), allocatable :: aerosol_relative_humidity
type(Atmosphere_t), target :: atm
real, dimension(:), allocatable :: average_cloud
real, dimension(:), allocatable :: average_high_cloud
real, dimension(:), allocatable :: average_low_cloud
real, dimension(:), allocatable :: average_mid_cloud
integer, dimension(4) :: axes
integer :: column
real, dimension(:,:), pointer :: convective_droplet_number
real, dimension(:,:), pointer :: convective_fraction
real, dimension(:,:), pointer :: convective_ice_content
real, dimension(:,:), pointer :: convective_liquid_content
integer, dimension(6) :: date
real :: dt
real, dimension(:), pointer :: diffuse_albedo_nir
real, dimension(:), pointer :: direct_albedo_nir
real, dimension(:), pointer :: diffuse_albedo_vis
real, dimension(:), pointer :: direct_albedo_vis
real, dimension(:,:), allocatable :: diffuse_surface_albedo
real, dimension(:,:), allocatable :: direct_surface_albedo
type(FluxDiagnostics), dimension(4) :: flux_diags
integer :: i
integer :: ios
integer :: id_lat
integer :: id_lon
integer :: id_pfull
integer :: id_phalf
type(domain2d), pointer :: io_domain
integer :: j
integer :: k
real, dimension(:), pointer :: land_fraction
integer :: last_infrared_band !Index of the last "nir" band.
real, dimension(:,:), pointer :: layer_pressure
real, dimension(:,:), pointer :: layer_temperature
real, dimension(:,:), pointer :: layer_thickness
real, dimension(:,:), pointer :: level_pressure
real, dimension(:,:), pointer :: level_temperature
integer :: logfile_handle
integer :: longwave_axis_id
real, dimension(:,:), allocatable :: longwave_band_limits
integer, dimension(:,:), allocatable :: longwave_gpoint_limits
integer :: num_bands
integer, dimension(1) :: num_columns
integer :: num_lat
integer :: num_layers
integer :: num_levels
integer :: num_lon
real, dimension(:,:), pointer :: ozone
type(RadiationContext) :: radiation_context
real, dimension(:,:), allocatable :: surface_emissivity !(lon, lat).
type(time_type) :: time
type(time_type) :: time_next
integer :: shortwave_axis_id
real, dimension(:,:), allocatable :: shortwave_band_limits
integer, dimension(:,:), allocatable :: shortwave_gpoint_limits
type(SolarSpectrum) :: solar_flux_spectrum
type(SolarConstant) :: solar_flux_constant
real, dimension(:,:), pointer :: stratiform_droplet_number
real, dimension(:,:), pointer :: stratiform_fraction
real, dimension(:,:), pointer :: stratiform_ice_content
real, dimension(:,:), pointer :: stratiform_liquid_content
real :: surface_albedo_weight !Weighting needed to combine "nir" and "vis" albedo values in
                              !in the band that contains the infrared cut-off.
real, dimension(:), pointer :: surface_temperature
integer :: t
real, dimension(:), allocatable :: total_solar_flux
real, dimension(:,:), pointer :: water_vapor
real, dimension(:), pointer :: zenith

!Runtime options.
integer :: near_infrared_cutoff = 14600 !Wavenumber [cm-1] that distinguishes the visible from the near-infrared.
character(len=256) :: solar_constant_path = ""
character(len=256) :: solar_spectrum_path = ""
namelist /standalone_radiation_nml/ near_infrared_cutoff, &
                                    solar_constant_path, &
                                    solar_spectrum_path

!Start up mpi.
call fms_init()
call time_manager_init()
call tracer_manager_init()
call grid_init()
call time_interp_external_init()

!Read the namelist.
read(input_nml_file, nml=standalone_radiation_nml, iostat=ios)
ios = check_nml_error(ios, "standalone_radiation_nml")
logfile_handle = stdlog()
if (mpp_pe() .eq. mpp_root_pe()) then
  write(logfile_handle, nml=standalone_radiation_nml)
endif

!Read in input data.
call create_atmosphere(atm)

!Set the model time.
if (trim(atm%calendar) .eq. "julian") then
  call set_calendar_type(julian)
else
  call error_mesg("main", "only julian calendar supported.", fatal)
endif
time = get_cal_time(atm%time(1), atm%time_units, atm%calendar)
if (atm%num_times .gt. 1) then
  dt = atm%time(2) - atm%time(1)
else
  dt = 24.*3000.
endif

!Read in the solar data.
call solar_flux_constant%create("solar_flux", trim(solar_constant_path))
call solar_flux_spectrum%create(trim(solar_spectrum_path))

!Get horizontal domain size.
io_domain => mpp_get_io_domain(atm%domain)
call mpp_get_compute_domain(io_domain, xsize=num_lon, ysize=num_lat)
num_columns(1) = num_lat*num_lon
num_layers = atm%num_layers
num_levels = num_layers + 1

!Allocate cloud arrays.
allocate(average_cloud(num_columns(1)))
allocate(average_high_cloud(num_columns(1)))
allocate(average_low_cloud(num_columns(1)))
allocate(average_mid_cloud(num_columns(1)))

!Allocate aerosol arrays.
allocate(aerosol_relative_humidity(num_lon, num_lat, num_layers))

!Initialize diag_manager.
call get_date(time, date(1), date(2), date(3), date(4), date(5), date(6))
call diag_manager_init(time_init=date)

!Initialize the radiation object.
call radiation_context%create(num_columns, num_layers, size(num_columns), solar_flux_spectrum%grid, &
                              solar_flux_spectrum%flux, atm%longitude_bounds, atm%latitude_bounds)

!Initialize the diagnostics.
id_lon = diag_axis_init("lon", atm%longitude, "degrees_east", "X", domain2=atm%domain)
id_lat = diag_axis_init("lat", atm%latitude, "degrees_north", "Y", domain2=atm%domain)
id_phalf = diag_axis_init("phalf", atm%level, "mb", "Z")
id_pfull = diag_axis_init("pfull", atm%layer, "mb", "Z")
axes = [id_lon, id_lat, id_pfull, id_phalf]
call radiation_context%longwave_gas_optics%band_limits(longwave_band_limits)
call radiation_context%shortwave_gas_optics%band_limits(shortwave_band_limits)
longwave_axis_id = diag_axis_init("longwave_band", longwave_band_limits(2, :), "cm-1", "Z")
shortwave_axis_id = diag_axis_init("shortwave_band", shortwave_band_limits(2, :), "cm-1", "Z")
do i = 1, size(flux_diags)
  call flux_diags(i)%create(i, time, axes, longwave_axis_id, shortwave_axis_id)
enddo

!Set the surface emissivity.
num_bands = radiation_context%longwave_gas_optics%num_bands()
allocate(surface_emissivity(num_bands, num_columns(1)))
surface_emissivity(:, :) = 1.

!Calculate variables need for dealing with the land spectral decomposition.
call radiation_context%shortwave_gas_optics%band_limits(shortwave_band_limits)
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

!Main loop.
do t = 1, atm%num_times
  !Calculate the current time.
  time = get_cal_time(atm%time(t), atm%time_units, atm%calendar)
  time_next = get_cal_time(atm%time(t) + dt, atm%time_units, atm%calendar)
  call print_time(time, "Running timestep: ")

  !Read in the atmospheric properies.
  call read_time_slice(atm, t)
  call solar_flux_constant%update(time)

  !Copy aerosol relative_humidity for now.
  aerosol_relative_humidity(:,:,:) = atm%ppmv(:,:,:,h2o)

  !Time interpolation.
  call radiation_context%update(time)

  !Update gas concentrations.
  water_vapor(1:num_columns(1), 1:num_layers) => atm%ppmv(1:num_lon, 1:num_lat, 1:num_layers, h2o)
  ozone(1:num_columns(1), 1:num_layers) => atm%ppmv(1:num_lon, 1:num_lat, 1:num_layers, o3)
  call radiation_context%update_concentrations(water_vapor, ozone, 1)

  !Calculate gas optics.
  layer_pressure(1:num_columns(1), 1:num_layers) => atm%layer_pressure(1:num_lon, 1:num_lat, 1:num_layers)
  level_pressure(1:num_columns(1), 1:num_levels) => atm%level_pressure(1:num_lon, 1:num_lat, 1:num_levels)
  layer_temperature(1:num_columns(1), 1:num_layers) => atm%layer_temperature(1:num_lon, 1:num_lat, 1:num_layers)
  level_temperature(1:num_columns(1), 1:num_levels) => atm%level_temperature(1:num_lon, 1:num_lat, 1:num_levels)
  surface_temperature(1:num_columns(1)) => atm%surface_temperature(1:num_lon, 1:num_lat)
  zenith(1:num_columns(1)) => atm%solar_zenith_angle(1:num_lon, 1:num_lat)
  call radiation_context%calculate_gas_optics(layer_pressure, level_pressure, &
                                              layer_temperature, level_temperature, &
                                              surface_temperature, 1, zenith)

  !Calculate cloud optics.
  land_fraction(1:num_columns(1)) => atm%land_fraction(1:num_lon, 1:num_lat)
  convective_droplet_number(1:num_columns(1), 1:num_layers) => atm%shallow_droplet_number(1:num_lon, 1:num_lat, 1:num_layers)
  convective_fraction(1:num_columns(1), 1:num_layers) => atm%shallow_cloud_fraction(1:num_lon, 1:num_lat, 1:num_layers)
  convective_liquid_content(1:num_columns(1), 1:num_layers) => atm%shallow_cloud_liquid_content(1:num_lon, 1:num_lat, 1:num_layers)
  convective_ice_content(1:num_columns(1), 1:num_layers) => atm%shallow_cloud_ice_content(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_droplet_number(1:num_columns(1), 1:num_layers) => atm%stratiform_droplet_number(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_fraction(1:num_columns(1), 1:num_layers) => atm%stratiform_cloud_fraction(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_liquid_content(1:num_columns(1), 1:num_layers) => atm%stratiform_cloud_liquid_content(1:num_lon, 1:num_lat, 1:num_layers)
  stratiform_ice_content(1:num_columns(1), 1:num_layers) => atm%stratiform_cloud_ice_content(1:num_lon, 1:num_lat, 1:num_layers)
  layer_thickness(1:num_columns(1), 1:num_layers) => atm%layer_thickness(1:num_lon, 1:num_lat, 1:num_layers)
  call radiation_context%calculate_cloud_optics(layer_pressure, layer_temperature, land_fraction, &
                                                layer_thickness, convective_droplet_number, &
                                                convective_fraction, convective_ice_content, &
                                                convective_liquid_content, stratiform_droplet_number, &
                                                stratiform_fraction, stratiform_ice_content, &
                                                stratiform_liquid_content, 1, level_pressure, &
                                                average_high_cloud, average_mid_cloud, average_low_cloud, &
                                                average_cloud)

  !Calculate aerosol optics.
  call radiation_context%calculate_aerosol_optics(atm%aerosols, aerosol_relative_humidity, &
                                                  atm%level_pressure, layer_thickness, time, &
                                                  1, num_lon, 1, num_lat, 1)

  !Calculate the surface albedo.
  diffuse_albedo_nir(1:num_columns(1)) => atm%surface_albedo_diffuse_ir(1:num_lon, 1:num_lat)
  direct_albedo_nir(1:num_columns(1)) => atm%surface_albedo_direct_ir(1:num_lon, 1:num_lat)
  diffuse_albedo_vis(1:num_columns(1)) => atm%surface_albedo_diffuse_uv(1:num_lon, 1:num_lat)
  direct_albedo_vis(1:num_columns(1)) => atm%surface_albedo_direct_uv(1:num_lon, 1:num_lat)
  num_bands = radiation_context%shortwave_gas_optics%num_bands()
  allocate(diffuse_surface_albedo(num_bands, num_columns(1)))
  allocate(direct_surface_albedo(num_bands, num_columns(1)))
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
  allocate(total_solar_flux(num_columns(1)))
  total_solar_flux = sum(radiation_context%toa(1)%flux, 2)
  do j = 1, num_lat
    do i = 1, num_lon
      column = (j - 1)*num_lon + i
      radiation_context%toa(1)%flux(column, :) = radiation_context%toa(1)%flux(column, :)* &
                                                 atm%earth_sun_distance_fraction* &
                                                 atm%daylight_fraction(i, j)*solar_flux_constant%flux/ &
                                                 total_solar_flux(column)
    enddo
  enddo
  deallocate(total_solar_flux)

  !Calculate the fluxes.
  !Clean-clear sky.
  call radiation_context%calculate_longwave_fluxes(surface_emissivity, clean_clear, 1)
  call radiation_context%calculate_shortwave_fluxes(zenith, direct_surface_albedo, &
                                                    diffuse_surface_albedo, clean_clear, 1)
  call radiation_context%longwave_gas_optics%gpoint_limits(longwave_gpoint_limits)
  call radiation_context%shortwave_gas_optics%gpoint_limits(shortwave_gpoint_limits)
  call flux_diags(clean_clear)%save_data(longwave_gpoint_limits, &
                                         radiation_context%longwave_fluxes(1), &
                                         shortwave_gpoint_limits, &
                                         radiation_context%shortwave_fluxes(1), &
                                         num_lon, num_lat, 1, 1, 1, time, level_pressure)

  !Clean sky.
  call radiation_context%calculate_longwave_fluxes(surface_emissivity, clean, 1)
  call radiation_context%calculate_shortwave_fluxes(zenith, direct_surface_albedo, &
                                                    diffuse_surface_albedo, clean, 1)
  call flux_diags(clean)%save_data(longwave_gpoint_limits, &
                                   radiation_context%longwave_fluxes(1), &
                                   shortwave_gpoint_limits, &
                                   radiation_context%shortwave_fluxes(1), &
                                   num_lon, num_lat, 1, 1, 1, time, level_pressure)

  !Clear sky.
  call radiation_context%calculate_longwave_fluxes(surface_emissivity, clear, 1)
  call radiation_context%calculate_shortwave_fluxes(zenith, direct_surface_albedo, &
                                                    diffuse_surface_albedo, clear, 1)
  call flux_diags(clear)%save_data(longwave_gpoint_limits, &
                                   radiation_context%longwave_fluxes(1), &
                                   shortwave_gpoint_limits, &
                                   radiation_context%shortwave_fluxes(1), &
                                   num_lon, num_lat, 1, 1, 1, time, level_pressure)

  !All sky.$
  call radiation_context%calculate_longwave_fluxes(surface_emissivity, all_, 1)
  call radiation_context%calculate_shortwave_fluxes(zenith, direct_surface_albedo, &
                                                    diffuse_surface_albedo, all_, 1)
  call flux_diags(all_)%save_data(longwave_gpoint_limits, &
                                  radiation_context%longwave_fluxes(1), &
                                  shortwave_gpoint_limits, &
                                  radiation_context%shortwave_fluxes(1), &
                                  num_lon, num_lat, 1, 1, 1, time, level_pressure)
  deallocate(diffuse_surface_albedo, direct_surface_albedo)
  deallocate(longwave_gpoint_limits, shortwave_gpoint_limits)

  !Write out diagnostics.
  call diag_manager_set_time_end(time)
  call diag_send_complete(time)
  time = time_next
enddo

!Clean up.
deallocate(average_cloud)
deallocate(average_high_cloud)
deallocate(average_low_cloud)
deallocate(average_mid_cloud)
deallocate(aerosol_relative_humidity)
call diag_manager_end(time)
call radiation_context%destroy()
call destroy_atmosphere(atm)
call tracer_manager_end()
call fms_end()

end program main
