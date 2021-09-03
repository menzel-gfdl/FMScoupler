program main
use am4, only: Atmosphere_t, create_atmosphere, destroy_atmosphere, h2o, o3
use constants_mod, only: pi
use diag_manager_mod, only: diag_axis_init, diag_manager_end, diag_manager_init
use field_manager_mod, only: model_atmos
use fms_mod, only: error_mesg, fatal, fms_end, fms_init
use get_cal_time_mod, only: get_cal_time
use mpp_domains_mod, only: domain2d, mpp_get_compute_domain, mpp_get_io_domain
use physics_radiation_exch_mod, only: clouds_from_moist_block_type
use radiation_driver2, only: Radiation
use time_manager_mod, only: julian, set_calendar_type, time_manager_init, time_type
use tracer_manager_mod, only: get_number_tracers, get_tracer_index, &
                              tracer_manager_end, tracer_manager_init
implicit none

integer :: i, id_lat, id_lon, id_pfull, id_phalf, iec, isc, &
           j, jec, jsc, num_tracers, nx, ny
real :: dlat, dlon, dt
type(Atmosphere_t) :: atm
type(Radiation) :: rte
type(time_type) :: time, time_next

real, dimension(:,:), allocatable :: latitude, longitude
real, dimension(:,:,:), allocatable :: aerosol_relative_humidity
real, dimension(:,:,:,:), allocatable :: aerosol_concentration
real, dimension(:,:,:,:), allocatable :: cloud_tracers
type(clouds_from_moist_block_type) :: clouds
type(domain2d), pointer :: io_domain

!Start up mpi.
call fms_init()
call time_manager_init()
call tracer_manager_init()

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
  dt = 0.
endif

!Calculate longitude and latitude grid box corners.
io_domain => mpp_get_io_domain(atm%domain)
call mpp_get_compute_domain(io_domain, xbegin=isc, xend=iec, xsize=nx, &
                            ybegin=jsc, yend=jec, ysize=ny)
if (nx .le. 1 .or. ny .le. 1) then
  call error_mesg("main", "must have more than 1 point in the x and y directions.", &
                  fatal)
endif
dlat = abs(atm%latitude(2) - atm%latitude(1))/2.
allocate(latitude(nx + 1, ny + 1))
latitude(:,1) = (2.*pi/real(ny + 1))*(atm%latitude(jsc) - dlat)
do i = 1, ny
  latitude(:,i + 1) = (2.*pi/real(ny + 1))*(atm%latitude(jsc + i - 1) + dlat)
enddo
dlon = abs(atm%longitude(2) - atm%longitude(1))/2.
allocate(longitude(nx + 1, ny + 1))
longitude(1,:) = (2.*pi/real(nx + 1))*(atm%longitude(isc) - dlon)
do i = 1, nx
  longitude(i + 1,:) = (2.*pi/real(nx + 1))*(atm%longitude(isc + i - 1) + dlon)
enddo

!Allocate cloud type arrays.
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
aerosol_relative_humidity(:,:,:) = 0.
allocate(aerosol_concentration(nx, ny, atm%num_layers, 15))
aerosol_concentration(:,:,:,:) = 0.

!Set up the diagnostics.
call diag_manager_init()
!call diag_manager_init(time_init=(/time%year, time%month, time%day, time%hour, &
!                                   time%minute, time%second/))
id_lon = diag_axis_init("lon", atm%longitude, "degrees_east", "X", domain2=atm%domain)
id_lat = diag_axis_init("lat", atm%latitude, "degrees_north", "Y", domain2=atm%domain)
id_phalf = diag_axis_init("phalf", atm%level, "mb", "Z")
id_pfull = diag_axis_init("pfull", atm%layer, "mb", "Z")

!Initialize the radiation object.
call rte%create(atm%num_layers, longitude, latitude)
call rte%initialize_diagnostics(id_lon, id_lat, id_pfull, id_phalf, time)

do i = 1, atm%num_times
  !Calculate the current time.
  time = get_cal_time(atm%time(i), atm%time_units, atm%calendar)
  time_next = get_cal_time(atm%time(i) + dt, atm%time_units, atm%calendar)

  !Populate the clouds derived type.
  j = clouds%index_uw_conv
  clouds%cloud_data(j)%droplet_number(:,:,:) = atm%shallow_droplet_number(:,:,:,i)
  clouds%cloud_data(j)%cloud_area(:,:,:) = atm%shallow_cloud_fraction(:,:,:,i)
  clouds%cloud_data(j)%liquid_amt(:,:,:) = atm%shallow_cloud_liquid_content(:,:,:,i)
  clouds%cloud_data(j)%ice_amt(:,:,:) = atm%shallow_cloud_ice_content(:,:,:,i)
  j = get_tracer_index(model_atmos, "liq_drp")
  cloud_tracers(:,:,:,j) = atm%stratiform_droplet_number(:,:,:,i)
  j = get_tracer_index(model_atmos, "cld_amt")
  cloud_tracers(:,:,:,j) = atm%stratiform_cloud_fraction(:,:,:,i)
  j = get_tracer_index(model_atmos, "liq_wat")
  cloud_tracers(:,:,:,j) = atm%stratiform_cloud_liquid_content(:,:,:,i)
  j = get_tracer_index(model_atmos, "ice_wat")
  cloud_tracers(:,:,:,j) = atm%stratiform_cloud_ice_content(:,:,:,i)
  j = clouds%index_strat
  clouds%cloud_data(j)%droplet_number(:,:,:) = atm%stratiform_droplet_number(:,:,:,i)
  clouds%cloud_data(j)%cloud_area(:,:,:) = atm%stratiform_cloud_fraction(:,:,:,i)
  clouds%cloud_data(j)%liquid_amt(:,:,:) = atm%stratiform_cloud_liquid_content(:,:,:,i)
  clouds%cloud_data(j)%ice_amt(:,:,:) = atm%stratiform_cloud_ice_content(:,:,:,i)

  call rte%update(time)
  call rte%calculate(atm%layer_pressure(:,:,:,i), &
                     atm%level_pressure(:,:,:,i), &
                     atm%layer_temperature(:,:,:,i), &
                     atm%level_temperature(:,:,:,i), &
                     atm%surface_temperature(:,:,i), &
                     aerosol_relative_humidity, &
                     atm%ppmv(:,:,:,i,h2o), &
                     atm%ppmv(:,:,:,i,o3), &
                     atm%layer_thickness(:,:,:,i), &
                     atm%solar_zenith_angle(:,:,i), &
                     atm%surface_albedo_direct_uv(:,:,i), &
                     atm%surface_albedo_diffuse_uv(:,:,i), &
                     atm%surface_albedo_direct_ir(:,:,i), &
                     atm%surface_albedo_diffuse_ir(:,:,i), &
                     atm%daylight_fraction(:,:,i), &
                     atm%earth_sun_distance_fraction(i), &
                     aerosol_concentration, &
                     clouds, &
                     cloud_tracers, &
                     atm%land_fraction(:,:,i), &
                     time, &
                     time_next)
  time = time_next
enddo

!Clean up.
deallocate(latitude, longitude, aerosol_relative_humidity)
deallocate(aerosol_concentration, cloud_tracers)
call diag_manager_end(time)
do i = 1, size(clouds%cloud_data)
  deallocate(clouds%cloud_data(i)%droplet_number)
  deallocate(clouds%cloud_data(i)%cloud_area)
  deallocate(clouds%cloud_data(i)%liquid_amt)
  deallocate(clouds%cloud_data(i)%ice_amt)
enddo
deallocate(clouds%cloud_data)
call rte%destroy()
call destroy_atmosphere(atm)
call tracer_manager_end()
call fms_end()

end program main
