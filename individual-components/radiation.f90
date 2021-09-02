program main
use am4, only: Atmosphere_t, create_atmosphere, destroy_atmosphere
use radiation_driver2, only: Radiation
use time_manager_mod, only: operator(+), time_type

integer :: i, num_layers, num_times
type(Atmosphere_t) :: atm
type(Radiation) :: rte
type(time_type) :: dt, time, time_next

!Start up mpi.
call mpp_init()

!Read in input data.
call create_atmosphere(atm)

call rte%create(atm%num_layers, longitude, latitude)

do i = 1, atm%num_times

  call rte%update(time)
  call rte%calculate(atm%layer_pressure(:,:,:,i), &
                     atm%level_pressure(:,:,:,i), &
                     atm%layer_temperature(:,:,:,i), &
                     atm%level_temperature(:,:,:,i), &
                     atm%surface_temperature(:,:,i), &
                     relative_humidity, &
                     atm%water_vapor(:,:,:,i), &
                     atm%ozone(:,:,:,i), &
                     atm%layer_thickness(:,:,:,i), &
                     atm%cosine_zenith(:,:,i), &
                     atm%direct_surface_albedo_in_visible(:,:,i), &
                     atm%diffuse_surface_albedo_in_visible(:,:,i), &
                     atm%direct_surface_albedo_in_infrared(:,:,i), &
                     atm%diffuse_surface_albedo_in_infrared(:,:,i), &
                     atm%dailylight_fraction(:,:,i), &
                     atm%earth_sun_distance_fraction(i), &
                     aerosol_concentration, &
                     clouds, &
                     cloud_tracers, &
                     atm%land_fraction(:,:,i), &
                     time, &
                     time_next)
  time = time_next
enddo

call rte%destroy()
call destroy_atmosphere(atm)
