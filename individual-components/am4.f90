!> @brief Reads in data from AM4 history files.
module am4
use fms_mod, only: check_nml_error, error_mesg, fatal, input_nml_file, mpp_npes
use fms2_io_mod, only: close_file, FmsNetcdfDomainFile_t, FmsNetcdfFile_t, &
                       get_dimension_size, get_variable_attribute, get_variable_size, &
                       open_file, read_data, register_axis
use mo_rte_kind, only: wp
use mpp_domains_mod, only: domain2d, mpp_define_mosaic, mpp_define_io_domain, &
                           mpp_get_compute_domain, mpp_get_io_domain

implicit none
private


type, public :: Atmosphere_t
  character(len=128) :: calendar !< Calendar type for the time axis.
  real(kind=wp), dimension(:,:,:), allocatable :: daylight_fraction !< Daylight correction factor (block_size, num_blocks, time).
  type(domain2d) :: domain !< 2d domain.
  real(kind=wp), dimension(:), allocatable :: earth_sun_distance_fraction !< Earth sun distance fraction.
  real(kind=wp), dimension(:,:,:), allocatable :: land_fraction !< Land fraction (block_size, num_blocks, time).
  real(kind=wp), dimension(:), allocatable :: latitude !< Latitude [degrees].
  real(kind=wp), dimension(:), allocatable :: layer !< Pressure [mb].
  real(kind=wp), dimension(:,:,:,:), allocatable :: layer_pressure !< Pressure [Pa] (block_size, layer, num_blocks, time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: layer_temperature !< Temperature [K] (block_size, layer, num_blocks, time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: layer_thickness !< Thickness [m] (block_size, layer, num_blocks, time).
  real(kind=wp), dimension(:), allocatable :: level !< Pressure [mb].
  real(kind=wp), dimension(:,:,:,:), allocatable :: level_pressure !< Pressure [Pa] (blocks_size level, num_blocks, time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: level_temperature !< Temperature [K] (block_size, level, num_blocks, time).
  real(kind=wp), dimension(:), allocatable :: longitude !< Longitude [degrees].
  integer :: num_layers !< Number of layers.
  integer :: num_levels !< Number of levels.
  integer :: num_times !< Number of times.
  real(kind=wp), dimension(:,:,:,:,:), allocatable :: ppmv !< Molecular abundancee (block_size, level, num_blocks, time, molecule).
  real(kind=wp), dimension(:,:,:,:), allocatable :: shallow_cloud_fraction !< Saturation volume fraction (block_size, layer, num_blocks, time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: shallow_cloud_ice_content !< Cloud ice water content [g m-3]  (block_size, layer, num_blocks, time)
  real(kind=wp), dimension(:,:,:,:), allocatable :: shallow_cloud_liquid_content !< Cloud liquid water content [g m-3] (block_size, layer, num_blocks, time)
  real(kind=wp), dimension(:,:,:,:), allocatable :: shallow_droplet_number !< Cloud liquid droplet number [km-1] (block_size, layer, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: solar_zenith_angle !< Solar zenith angle [degrees] (block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: stratiform_cloud_fraction !< Saturation volume fraction (block_size, layer, num_blocks, time).
  real(kind=wp), dimension(:,:,:,:), allocatable :: stratiform_cloud_ice_content !< Cloud ice water content [g m-3]  (block_size, layer, num_blocks, time)
  real(kind=wp), dimension(:,:,:,:), allocatable :: stratiform_cloud_liquid_content !< Cloud liquid water content [g m-3] (block_size, layer, num_blocks, time)
  real(kind=wp), dimension(:,:,:,:), allocatable :: stratiform_droplet_number !< Cloud liquid droplet number [km-1] (block_size, layer, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: surface_albedo_diffuse_ir !< Surface albedo for infrared diffuse beam (block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: surface_albedo_diffuse_uv !< Surface albedo for ultraviolet diffuse beam (block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: surface_albedo_direct_ir !< Surface albedo for infrared direct beam (block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: surface_albedo_direct_uv !< Surface albedo for ultraviolet direct beam (block_size, num_blocks, time).
  real(kind=wp), dimension(:,:,:), allocatable :: surface_temperature !< Surface temperature [K] (block_size, num_blocks, time).
  real(kind=wp), dimension(:), allocatable :: time !< Time [hours or days or months].
  character(len=128) :: time_units !< Time units (i.e. days since 0000-00-00 00:00:00)
  real(kind=wp), dimension(:), allocatable :: total_solar_irradiance !< Total solar irradiance [W m-2] (time).
end type Atmosphere_t

public :: create_atmosphere
public :: destroy_atmosphere
integer, parameter, public :: h2o = 1
integer, parameter, public :: o3 = 2

character(len=256) :: atmos_path = "none"
logical :: clearsky = .false.
namelist /am4_nml/ atmos_path, clearsky


contains


!> @brief Reserve memory and read in atmospheric data.
subroutine create_atmosphere(atm)

  type(Atmosphere_t), intent(inout) :: atm

  character(len=256) :: attr
  integer, dimension(4) :: dim_sizes
  integer :: err, i, ierr, nx, ny
  real, dimension(:,:), allocatable :: buffer2d
  type(FmsNetcdfDomainFile_t) :: dataset
  type(FmsNetcdfFile_t) :: tilefile
  type(domain2d), pointer :: io_domain

  read(input_nml_file, am4_nml, iostat=ierr)
  err = check_nml_error(ierr, "am4_nml")

  !Open the *.tile1.nc file and get the longitude and latitude sizes.
  i = len_trim(atmos_path)
  if (.not. open_file(tilefile, atmos_path(:i - 3)//".tile1.nc", "read")) then
    call error_mesg("create_atmosphere", &
                    "cannot find the "//trim(atmos_path)//" tile files.", fatal)
  endif
  call get_dimension_size(tilefile, "grid_xt", nx)
  call get_dimension_size(tilefile, "grid_yt", ny)
  call close_file(tilefile)

  !Create a 2d domain.
  call create_atmosphere_domain(nx, ny, atm%domain)

  !Open dataset.
  if (.not. open_file(dataset, atmos_path, "read", atm%domain)) then
    call error_mesg("create_atmosphere", "cannot open file "//trim(atmos_path)//".", &
                    fatal)
  endif

  !Get domain sizes.
  io_domain => mpp_get_io_domain(atm%domain)
  if (.not. associated(io_domain)) then
    call error_mesg("create_atmosphere", "I/O domain doesn't exist.", fatal)
  endif
  call mpp_get_compute_domain(io_domain, xsize=nx, ysize=ny)

  !Identify longitude and latitude dimensions.
  call register_axis(dataset, "grid_xt", "x")
  call register_axis(dataset, "grid_yt", "y")

  !Read in axes so they can be passed to diag_manager.
  call get_dimension_size(dataset, "grid_xt", i)
  allocate(atm%longitude(i))
  call read_data(dataset, "grid_xt", atm%longitude)
  call get_dimension_size(dataset, "grid_yt", i)
  allocate(atm%latitude(i))
  call read_data(dataset, "grid_yt", atm%latitude)
  call get_dimension_size(dataset, "pfull", i)
  allocate(atm%layer(i))
  call read_data(dataset, "pfull", atm%layer)
  call get_dimension_size(dataset, "phalf", i)
  allocate(atm%level(i))
  call read_data(dataset, "phalf", atm%level)

  !Read in the times.
  call get_variable_size(dataset, "time", dim_sizes(1:1))
  atm%num_times = dim_sizes(1)
  allocate(atm%time(atm%num_times))
  call read_data(dataset, "time", atm%time)
  call get_variable_attribute(dataset, "time", "units", atm%time_units)
  call get_variable_attribute(dataset, "time", "calendar", atm%calendar)

  !Surface temperature.
  allocate(atm%surface_temperature(nx, ny, atm%num_times))
  call read_data(dataset, "surface_temperature", atm%surface_temperature)

  !Land fraction.
  allocate(atm%land_fraction(nx, ny, atm%num_times))
  call read_data(dataset, "land_fraction", atm%land_fraction)

  !Solar zenith angle.
  allocate(atm%solar_zenith_angle(nx, ny, atm%num_times))
  call read_data(dataset, "cosine_zenith", atm%solar_zenith_angle)

  !Solar constant.
  allocate(atm%total_solar_irradiance(atm%num_times))
  allocate(buffer2d(1, atm%num_times))
  call read_data(dataset, "solar_constant", buffer2d)
  atm%total_solar_irradiance(:) = buffer2d(1,:)
  deallocate(buffer2d)

  !Daylight fraction.
  allocate(atm%daylight_fraction(nx, ny, atm%num_times))
  call read_data(dataset, "daylight_fraction", atm%daylight_fraction)

  !Earth sun distance fraction.
  allocate(atm%earth_sun_distance_fraction(atm%num_times))
  allocate(buffer2d(1, atm%num_times))
  call read_data(dataset, "earth_sun_distance_fraction", buffer2d)
  atm%earth_sun_distance_fraction(:) = buffer2d(1,:)
  deallocate(buffer2d)

  !Get the number of layers and levels.
  call get_dimension_size(dataset, "pfull", atm%num_layers)
  atm%num_levels = atm%num_layers + 1

  !Pressure.
  allocate(atm%layer_pressure(nx, ny, atm%num_layers, atm%num_times))
  call read_data(dataset, "layer_pressure", atm%layer_pressure)
  allocate(atm%level_pressure(nx, ny, atm%num_levels, atm%num_times))
  call read_data(dataset, "level_pressure", atm%level_pressure)

  !Temperature.
  allocate(atm%layer_temperature(nx, ny, atm%num_layers, atm%num_times))
  call read_data(dataset, "layer_temperature", atm%layer_temperature)
  allocate(atm%level_temperature(nx, ny, atm%num_levels, atm%num_times))
  call read_data(dataset, "level_temperature", atm%level_temperature)

  !Molecular abundances.
  allocate(atm%ppmv(nx, ny, atm%num_layers, atm%num_times, 2))

  !Water abundance.
  call read_data(dataset, "water_vapor", atm%ppmv(:,:,:,:,h2o))

  !Read water vapor and ozone.
  call read_data(dataset, "ozone", atm%ppmv(:,:,:,:,o3))

  !Read in the surface albedos.
  allocate(atm%surface_albedo_direct_uv(nx, ny, atm%num_times))
  call read_data(dataset, "visible_direct_albedo", atm%surface_albedo_direct_uv)
  allocate(atm%surface_albedo_diffuse_uv(nx, ny, atm%num_times))
  call read_data(dataset, "visible_diffuse_albedo", atm%surface_albedo_diffuse_uv)
  allocate(atm%surface_albedo_direct_ir(nx, ny, atm%num_times))
  call read_data(dataset, "infrared_direct_albedo", atm%surface_albedo_direct_ir)
  allocate(atm%surface_albedo_diffuse_ir(nx, ny, atm%num_times))
  call read_data(dataset, "infrared_diffuse_albedo", atm%surface_albedo_diffuse_ir)

  if (.not. clearsky) then
    !Layer thickness.
    allocate(atm%layer_thickness(nx, ny, atm%num_layers, atm%num_times))
    call read_data(dataset, "layer_thickness", atm%layer_thickness)

    !Stratiform cloud fraction.
    allocate(atm%stratiform_cloud_fraction(nx, ny, atm%num_layers, atm%num_times))
    call read_data(dataset, "stratiform_cloud_fraction", atm%stratiform_cloud_fraction)

    !Shallow cloud fraction.
    allocate(atm%shallow_cloud_fraction(nx, ny, atm%num_layers, atm%num_times))
    call read_data(dataset, "shallow_cloud_fraction", atm%shallow_cloud_fraction)

    !Stratiform ice content.
    allocate(atm%stratiform_cloud_ice_content(nx, ny, atm%num_layers, atm%num_times))
    call read_data(dataset, "stratiform_ice_content", atm%stratiform_cloud_ice_content)

    !Shallow ice content.
    allocate(atm%shallow_cloud_ice_content(nx, ny, atm%num_layers, atm%num_times))
    call read_data(dataset, "shallow_ice_content", atm%shallow_cloud_ice_content)

    !Stratiform liquid content.
    allocate(atm%stratiform_cloud_liquid_content(nx, ny, atm%num_layers, atm%num_times))
    call read_data(dataset, "stratiform_liquid_content", atm%stratiform_cloud_liquid_content)

    !Shallow liquid content.
    allocate(atm%shallow_cloud_liquid_content(nx, ny, atm%num_layers, atm%num_times))
    call read_data(dataset, "shallow_liquid_content", atm%shallow_cloud_liquid_content)

    !Stratiform droplet number.
    allocate(atm%stratiform_droplet_number(nx, ny, atm%num_layers, atm%num_times))
    call read_data(dataset, "stratiform_droplet_number", atm%stratiform_droplet_number)

    !Shallow droplet number.
    allocate(atm%shallow_droplet_number(nx, ny, atm%num_layers, atm%num_times))
    call read_data(dataset, "shallow_droplet_number", atm%shallow_droplet_number)
  endif
  call close_file(dataset)
end subroutine create_atmosphere


!> @brief Free memory for atmosphere.
subroutine destroy_atmosphere(atm)

  type(Atmosphere_t), intent(inout) :: atm

  if (allocated(atm%shallow_cloud_fraction)) deallocate(atm%shallow_cloud_fraction)
  if (allocated(atm%stratiform_cloud_fraction)) deallocate(atm%stratiform_cloud_fraction)
  if (allocated(atm%shallow_cloud_ice_content)) deallocate(atm%shallow_cloud_ice_content)
  if (allocated(atm%stratiform_cloud_ice_content)) deallocate(atm%stratiform_cloud_ice_content)
  if (allocated(atm%shallow_cloud_liquid_content)) deallocate(atm%shallow_cloud_liquid_content)
  if (allocated(atm%stratiform_cloud_liquid_content)) deallocate(atm%stratiform_cloud_liquid_content)
  if (allocated(atm%daylight_fraction)) deallocate(atm%daylight_fraction)
  if (allocated(atm%latitude)) deallocate(atm%latitude)
  if (allocated(atm%layer_pressure)) deallocate(atm%layer_pressure)
  if (allocated(atm%layer_temperature)) deallocate(atm%layer_temperature)
  if (allocated(atm%layer_thickness)) deallocate(atm%layer_thickness)
  if (allocated(atm%level)) deallocate(atm%level)
  if (allocated(atm%level_pressure)) deallocate(atm%level_pressure)
  if (allocated(atm%level_temperature)) deallocate(atm%level_temperature)
  if (allocated(atm%longitude)) deallocate(atm%longitude)
  if (allocated(atm%ppmv)) deallocate(atm%ppmv)
  if (allocated(atm%solar_zenith_angle)) deallocate(atm%solar_zenith_angle)
  if (allocated(atm%surface_albedo_diffuse_ir)) deallocate(atm%surface_albedo_diffuse_ir)
  if (allocated(atm%surface_albedo_diffuse_uv)) deallocate(atm%surface_albedo_diffuse_uv)
  if (allocated(atm%surface_albedo_direct_ir)) deallocate(atm%surface_albedo_direct_ir)
  if (allocated(atm%surface_albedo_direct_uv)) deallocate(atm%surface_albedo_direct_uv)
  if (allocated(atm%surface_temperature)) deallocate(atm%surface_temperature)
  if (allocated(atm%time)) deallocate(atm%time)
  if (allocated(atm%total_solar_irradiance)) deallocate(atm%total_solar_irradiance)
  if (allocated(atm%stratiform_droplet_number)) deallocate(atm%stratiform_droplet_number)
  if (allocated(atm%shallow_droplet_number)) deallocate(atm%shallow_droplet_number)
  if (allocated(atm%land_fraction)) deallocate(atm%land_fraction)
end subroutine destroy_atmosphere


!> @brief Create a 2d domain.
subroutine create_atmosphere_domain(nx, ny, domain)

  integer, intent(in) :: nx !< Number of points per tile in the x-direction.
  integer, intent(in) :: ny !< Number of points per tile in the y-direction.
  type(domain2d), intent(inout) :: domain !< A cubed-sphere domain.

  integer, dimension(4,6) :: global_indices
  integer, dimension(2,6) :: layout
  integer, dimension(6) :: ni, nj, pe_start, pe_end
  integer, dimension(2) :: io_layout, msize
  integer, dimension(12) :: tile1, tile2, istart1, iend1, jstart1, jend1, istart2, &
                            iend2, jstart2, jend2
  integer :: i, npes
  integer, parameter :: ntiles = 6
  integer, parameter :: num_contact = 12
  integer, parameter :: ehalo = 2
  integer, parameter :: whalo = 2
  integer, parameter :: nhalo = 2
  integer, parameter :: shalo = 2

  npes = mpp_npes()
  if (mod(npes, ntiles) .ne. 0) then
    call error_mesg("create_atmosphere_domain", &
                    "number or ranks must be a multiple of 6.", fatal)
  endif
  do i = 1, ntiles
    ni(i) = nx
    nj(i) = ny
    global_indices(:,i) = (/1, nx, 1, ny/)
    layout(:,i) = (/1, npes/ntiles/)
    pe_start(i) = (i - 1)*npes/ntiles
    pe_end(i) = i*npes/ntiles - 1
  enddo
  io_layout(:) = 1

  !Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
  tile1(1) = 1
  tile2(1) = 2
  istart1(1) = ni(1)
  iend1(1) = ni(1)
  jstart1(1) = 1
  jend1(1) = nj(1)
  istart2(1) = 1
  iend2(1) = 1
  jstart2(1) = 1
  jend2(1) = nj(2)

  !Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
  tile1(2) = 1
  tile2(2) = 3
  istart1(2) = 1
  iend1(2) = ni(1)
  jstart1(2) = nj(1)
  jend1(2) = nj(1)
  istart2(2) = 1
  iend2(2) = 1
  jstart2(2) = nj(3)
  jend2(2) = 1

  !Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
  tile1(3) = 1
  tile2(3) = 5
  istart1(3) = 1
  iend1(3) = 1
  jstart1(3) = 1
  jend1(3) = nj(1)
  istart2(3) = ni(5)
  iend2(3) = 1
  jstart2(3) = nj(5)
  jend2(3) = nj(5)

  !Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
  tile1(4) = 1
  tile2(4) = 6
  istart1(4) = 1
  iend1(4) = ni(1)
  jstart1(4) = 1
  jend1(4) = 1
  istart2(4) = 1
  iend2(4) = ni(6)
  jstart2(4) = nj(6)
  jend2(4) = nj(6)

  !Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
  tile1(5) = 2
  tile2(5) = 3
  istart1(5) = 1
  iend1(5) = ni(2)
  jstart1(5) = nj(2)
  jend1(5) = nj(2)
  istart2(5) = 1
  iend2(5) = ni(3)
  jstart2(5) = 1
  jend2(5) = 1

  !Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
  tile1(6) = 2
  tile2(6) = 4
  istart1(6) = ni(2)
  iend1(6) = ni(2)
  jstart1(6) = 1
  jend1(6) = nj(2)
  istart2(6) = ni(4)
  iend2(6) = 1
  jstart2(6) = 1
  jend2(6) = 1

  !Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
  tile1(7) = 2
  tile2(7) = 6
  istart1(7) = 1
  iend1(7) = ni(2)
  jstart1(7) = 1
  jend1(7) = 1
  istart2(7) = ni(6)
  iend2(7) = ni(6)
  jstart2(7) = nj(6)
  jend2(7) = 1

  !Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
  tile1(8) = 3
  tile2(8) = 4
  istart1(8) = ni(3)
  iend1(8) = ni(3)
  jstart1(8) = 1
  jend1(8) = nj(3)
  istart2(8) = 1
  iend2(8) = 1
  jstart2(8) = 1
  jend2(8) = nj(4)

  !Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
  tile1(9) = 3
  tile2(9) = 5
  istart1(9) = 1
  iend1(9) = ni(3)
  jstart1(9) = nj(3)
  jend1(9) = nj(3)
  istart2(9) = 1
  iend2(9) = 1
  jstart2(9) = nj(5)
  jend2(9) = 1

  !Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
  tile1(10) = 4
  tile2(10) = 5
  istart1(10) = 1
  iend1(10) = ni(4)
  jstart1(10) = nj(4)
  jend1(10) = nj(4)
  istart2(10) = 1
  iend2(10) = ni(5)
  jstart2(10) = 1
  jend2(10) = 1

  !Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
  tile1(11) = 4
  tile2(11) = 6
  istart1(11) = ni(4)
  iend1(11) = ni(4)
  jstart1(11) = 1
  jend1(11) = nj(4)
  istart2(11) = ni(6)
  iend2(11) = 1
  jstart2(11) = 1
  jend2(11) = 1

  !Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
  tile1(12) = 5
  tile2(12) = 6
  istart1(12) = ni(5)
  iend1(12) = ni(5)
  jstart1(12) = 1
  jend1(12) = nj(5)
  istart2(12) = 1
  iend2(12) = 1
  jstart2(12) = 1
  jend2(12) = nj(6)
  msize(1) = maxval(ni(:)/layout(1,:)) + whalo + ehalo + 1
  msize(2) = maxval(nj(:)/layout(2,:)) + shalo + nhalo + 1
  call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contact, tile1, &
                         tile2, istart1, iend1, jstart1, jend1, istart2, iend2, &
                         jstart2, jend2, pe_start, pe_end, symmetry=.true., &
                         whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                         name=trim("Cubed-sphere"), memory_size=msize)
  call mpp_define_io_domain(domain, io_layout)
end subroutine create_atmosphere_domain


end module am4
