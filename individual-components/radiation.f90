program main
use am4, only: Atmosphere_t, create_atmosphere, destroy_atmosphere, h2o, o3, &
               read_time_slice
use constants_mod, only: pi
use diag_manager_mod, only: diag_axis_init, diag_manager_end, diag_manager_init, &
                            diag_manager_set_time_end, diag_send_complete
use field_manager_mod, only: model_atmos
use fms_mod, only: error_mesg, fatal, fms_end, fms_init
use get_cal_time_mod, only: get_cal_time
use grid2_mod, only: grid_end, grid_init
use mpp_domains_mod, only: domain2d, mpp_get_compute_domain, mpp_get_io_domain
use physics_radiation_exch_mod, only: clouds_from_moist_block_type
use radiation_driver2, only: Radiation
use time_manager_mod, only: get_date, julian, print_time, set_calendar_type, time_manager_init, &
                            time_type
use tracer_manager_mod, only: get_number_tracers, get_tracer_index, &
                              tracer_manager_end, tracer_manager_init
implicit none

integer :: i, id_lat, id_lon, id_pfull, id_phalf, j, num_tracers, nx, ny
integer, dimension(6) :: date
real :: dt
type(Atmosphere_t) :: atm
type(Radiation) :: rte
type(time_type) :: time, time_next
real, dimension(:,:,:), allocatable :: aerosol_relative_humidity
real, dimension(:,:,:,:), allocatable :: aerosol_concentration
real, dimension(:,:,:,:), allocatable :: cloud_tracers
type(clouds_from_moist_block_type) :: clouds
type(domain2d), pointer :: io_domain
real, dimension(:,:,:), allocatable :: longwave_all_sky_flux_down
real, dimension(:,:,:), allocatable :: longwave_all_sky_flux_up
real, dimension(:,:,:), allocatable :: longwave_heating_rate
real, dimension(:,:,:), allocatable :: shortwave_all_sky_flux_down
real, dimension(:,:,:), allocatable :: shortwave_all_sky_flux_up
real, dimension(:,:,:), allocatable :: shortwave_all_sky_flux_direct
real, dimension(:,:,:), allocatable :: shortwave_heating_rate
real, dimension(:,:,:), allocatable :: shortwave_flux_up_vis
real, dimension(:,:,:), allocatable :: shortwave_flux_down_vis
real, dimension(:,:,:), allocatable :: shortwave_flux_direct_vis

!Start up mpi.
call fms_init()
call time_manager_init()
call tracer_manager_init()
call grid_init()

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

!Allocate cloud type arrays.
io_domain => mpp_get_io_domain(atm%domain)
call mpp_get_compute_domain(io_domain, xsize=nx, ysize=ny)
clouds%index_uw_conv = 1
clouds%index_strat = 2
allocate(clouds%cloud_data(2))
do i = 1, size(clouds%cloud_data)
  allocate(clouds%cloud_data(i)%droplet_number(nx, ny, atm%num_layers))
  allocate(clouds%cloud_data(i)%cloud_area(nx, ny, atm%num_layers))
  allocate(clouds%cloud_data(i)%liquid_amt(nx, ny, atm%num_layers))
  allocate(clouds%cloud_data(i)%ice_amt(nx, ny, atm%num_layers))
enddo
call get_number_tracers(model_atmos, num_tracers)
allocate(cloud_tracers(nx, ny, atm%num_layers, num_tracers))
cloud_tracers(:,:,:,:) = 0.

!Allocate aerosol arrays.
allocate(aerosol_relative_humidity(nx, ny, atm%num_layers))
allocate(aerosol_concentration(nx, ny, atm%num_layers, num_tracers))
aerosol_concentration(:,:,:,:) = 0.

!Set up the diagnostics.
call get_date(time, date(1), date(2), date(3), date(4), date(5), date(6))
call diag_manager_init(time_init=date)
id_lon = diag_axis_init("lon", atm%longitude, "degrees_east", "X", domain2=atm%domain)
id_lat = diag_axis_init("lat", atm%latitude, "degrees_north", "Y", domain2=atm%domain)
id_phalf = diag_axis_init("phalf", atm%level, "mb", "Z")
id_pfull = diag_axis_init("pfull", atm%layer, "mb", "Z")

!Initialize the radiation object.
call rte%create(atm%num_layers, atm%longitude_bounds, atm%latitude_bounds)
call rte%initialize_diagnostics(id_lon, id_lat, id_pfull, id_phalf, time)

!Allocate output arrays.
allocate(longwave_all_sky_flux_down(nx, ny, atm%num_layers+1))
allocate(longwave_all_sky_flux_up(nx, ny, atm%num_layers+1))
allocate(longwave_heating_rate(nx, ny, atm%num_layers))
allocate(shortwave_all_sky_flux_up(nx, ny, atm%num_layers+1))
allocate(shortwave_all_sky_flux_down(nx, ny, atm%num_layers+1))
allocate(shortwave_all_sky_flux_direct(nx, ny, atm%num_layers+1))
allocate(shortwave_heating_rate(nx, ny, atm%num_layers))
allocate(shortwave_flux_up_vis(nx, ny, atm%num_layers+1))
allocate(shortwave_flux_down_vis(nx, ny, atm%num_layers+1))
allocate(shortwave_flux_direct_vis(nx, ny, atm%num_layers+1))

do i = 1, atm%num_times
  !Calculate the current time.
  time = get_cal_time(atm%time(i), atm%time_units, atm%calendar)
  time_next = get_cal_time(atm%time(i) + dt, atm%time_units, atm%calendar)
  call print_time(time, "Running timestep: ")

  !Read in the atmospheric properies.
  call read_time_slice(atm, i)

  !Populate the clouds derived type.
  j = clouds%index_uw_conv
  clouds%cloud_data(j)%droplet_number(:,:,:) = atm%shallow_droplet_number(:,:,:)
  clouds%cloud_data(j)%cloud_area(:,:,:) = atm%shallow_cloud_fraction(:,:,:)
  clouds%cloud_data(j)%liquid_amt(:,:,:) = atm%shallow_cloud_liquid_content(:,:,:)
  clouds%cloud_data(j)%ice_amt(:,:,:) = atm%shallow_cloud_ice_content(:,:,:)
  j = get_tracer_index(model_atmos, "liq_drp")
  cloud_tracers(:,:,:,j) = atm%stratiform_droplet_number(:,:,:)
  j = get_tracer_index(model_atmos, "cld_amt")
  cloud_tracers(:,:,:,j) = atm%stratiform_cloud_fraction(:,:,:)
  j = get_tracer_index(model_atmos, "liq_wat")
  cloud_tracers(:,:,:,j) = atm%stratiform_cloud_liquid_content(:,:,:)
  j = get_tracer_index(model_atmos, "ice_wat")
  cloud_tracers(:,:,:,j) = atm%stratiform_cloud_ice_content(:,:,:)
  j = clouds%index_strat
  clouds%cloud_data(j)%droplet_number(:,:,:) = atm%stratiform_droplet_number(:,:,:)
  clouds%cloud_data(j)%cloud_area(:,:,:) = atm%stratiform_cloud_fraction(:,:,:)
  clouds%cloud_data(j)%liquid_amt(:,:,:) = atm%stratiform_cloud_liquid_content(:,:,:)
  clouds%cloud_data(j)%ice_amt(:,:,:) = atm%stratiform_cloud_ice_content(:,:,:)

  !Copy aerosol relative_humidity for now.
  aerosol_relative_humidity(:,:,:) = atm%ppmv(:,:,:,h2o)

  call rte%update(time)
  call rte%calculate(atm%layer_pressure(:,:,:), &
                     atm%level_pressure(:,:,:), &
                     atm%layer_temperature(:,:,:), &
                     atm%level_temperature(:,:,:), &
                     atm%surface_temperature(:,:), &
                     aerosol_relative_humidity, &
                     atm%ppmv(:,:,:,h2o), &
                     atm%ppmv(:,:,:,o3), &
                     atm%layer_thickness(:,:,:), &
                     atm%solar_zenith_angle(:,:), &
                     atm%surface_albedo_direct_uv(:,:), &
                     atm%surface_albedo_diffuse_uv(:,:), &
                     atm%surface_albedo_direct_ir(:,:), &
                     atm%surface_albedo_diffuse_ir(:,:), &
                     atm%daylight_fraction(:,:), &
                     atm%earth_sun_distance_fraction, &
                     aerosol_concentration, &
                     clouds, &
                     cloud_tracers, &
                     atm%land_fraction(:,:), &
                     time, &
                     time_next, &
                     longwave_all_sky_flux_down, &
                     longwave_all_sky_flux_up, &
                     longwave_heating_rate, &
                     shortwave_all_sky_flux_down, &
                     shortwave_all_sky_flux_up, &
                     shortwave_all_sky_flux_direct, &
                     shortwave_heating_rate, &
                     shortwave_flux_up_vis, &
                     shortwave_flux_down_vis, &
                     shortwave_flux_direct_vis)
  call diag_manager_set_time_end(time)
  call diag_send_complete(time)
  time = time_next
enddo

!Clean up.
deallocate(aerosol_relative_humidity, aerosol_concentration, cloud_tracers)
call diag_manager_end(time)
do i = 1, size(clouds%cloud_data)
  deallocate(clouds%cloud_data(i)%droplet_number)
  deallocate(clouds%cloud_data(i)%cloud_area)
  deallocate(clouds%cloud_data(i)%liquid_amt)
  deallocate(clouds%cloud_data(i)%ice_amt)
enddo
deallocate(clouds%cloud_data)
deallocate(longwave_all_sky_flux_down)
deallocate(longwave_all_sky_flux_up)
deallocate(shortwave_all_sky_flux_down)
deallocate(shortwave_all_sky_flux_up)
deallocate(longwave_heating_rate)
deallocate(shortwave_heating_rate)
deallocate(shortwave_flux_up_vis)
deallocate(shortwave_flux_down_vis)
deallocate(shortwave_flux_direct_vis)
call rte%destroy()
call destroy_atmosphere(atm)
call tracer_manager_end()
call fms_end()

end program main
