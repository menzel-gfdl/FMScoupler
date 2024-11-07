!> @brief Reads in data from AM4 history files.
module am4
use block_control_mod, only: block_control_type, define_blocks
use constants_mod, only: pi
use fms_mod, only: check_nml_error, error_mesg, fatal, input_nml_file, mpp_npes
use fms2_io_mod, only: close_file, FmsNetcdfDomainFile_t, FmsNetcdfFile_t, &
                       get_dimension_size, get_variable_attribute, get_variable_size, &
                       open_file, read_data, register_axis
use grid2_mod, only: get_grid_cell_vertices
use mo_rte_kind, only: wp
use mpp_domains_mod, only: domain2d, mpp_define_mosaic, mpp_define_io_domain, &
                           mpp_get_compute_domain, mpp_get_current_ntile, &
                           mpp_get_io_domain, mpp_get_tile_id
implicit none
private


type, public :: Atmosphere_t
  character(len=128) :: calendar !< Calendar type for the time axis.
  real(kind=wp), dimension(:,:), allocatable :: daylight_fraction !< Daylight correction factor (lon, lat).
  type(domain2d) :: domain !< 2d domain.
  real(kind=wp) :: earth_sun_distance_fraction !< Earth sun distance fraction.
  real(kind=wp), dimension(:,:), allocatable :: land_fraction !< Land fraction (lon, lat).
  real(kind=wp), dimension(:), allocatable :: latitude !< Y dimension data (lat).
  real(kind=wp), dimension(:,:), allocatable :: latitude_bounds !< Latitude of cell vertices [degrees] (lon, lat).
  real(kind=wp), dimension(:), allocatable :: layer !< Pressure [mb] (layer).
  real(kind=wp), dimension(:,:,:), allocatable :: layer_pressure !< Pressure [Pa] (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: layer_temperature !< Temperature [K] (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: layer_thickness !< Thickness [m] (lon, lat, layer).
  real(kind=wp), dimension(:), allocatable :: level !< Pressure [mb] (level).
  real(kind=wp), dimension(:,:,:), allocatable :: level_pressure !< Pressure [Pa] (lon, lat, level).
  real(kind=wp), dimension(:,:,:), allocatable :: level_temperature !< Temperature [K] (lon, lat, level).
  real(kind=wp), dimension(:), allocatable :: longitude !< X dimension data (lon).
  real(kind=wp), dimension(:,:), allocatable :: longitude_bounds !< Longitude of cell vertices [degrees] (lon, lat).
  integer :: num_layers !< Number of layers.
  integer :: num_levels !< Number of levels.
  integer :: num_times !< Number of times.
  real(kind=wp), dimension(:,:,:,:), allocatable :: ppmv !< Molecular abundancee (lon, lat, layer, molecule).
  real(kind=wp), dimension(:,:,:), allocatable :: shallow_cloud_fraction !< Saturation volume fraction (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: shallow_cloud_ice_content !< Cloud ice water content [g m-3]  (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: shallow_cloud_liquid_content !< Cloud liquid water content [g m-3] (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: shallow_droplet_number !< Cloud liquid droplet number [km-1] (lon, lat, layer).
  real(kind=wp), dimension(:,:), allocatable :: solar_zenith_angle !< Solar zenith angle [degrees] (on, lat).
  real(kind=wp), dimension(:,:,:), allocatable :: stratiform_cloud_fraction !< Saturation volume fraction (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: stratiform_cloud_ice_content !< Cloud ice water content [g m-3]  (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: stratiform_cloud_liquid_content !< Cloud liquid water content [g m-3] (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: stratiform_droplet_number !< Cloud liquid droplet number [km-1] (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: stratiform_ice_number !< Cloud ice crystal number [???] (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: stratiform_rain !< Rain droplet number [???] (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: stratiform_rain_size !< Rain droplet size [microns] (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: stratiform_snow !< Snow crystal number [???] (lon, lat, layer).
  real(kind=wp), dimension(:,:,:), allocatable :: stratiform_snow_size !< Snow crystal sizse [microns] (lon, lat, layer).
  real(kind=wp), dimension(:,:), allocatable :: surface_albedo_diffuse_ir !< Surface albedo for infrared diffuse beam (lon, lat).
  real(kind=wp), dimension(:,:), allocatable :: surface_albedo_diffuse_uv !< Surface albedo for ultraviolet diffuse beam (lon, lat).
  real(kind=wp), dimension(:,:), allocatable :: surface_albedo_direct_ir !< Surface albedo for infrared direct beam (lon, lat).
  real(kind=wp), dimension(:,:), allocatable :: surface_albedo_direct_uv !< Surface albedo for ultraviolet direct beam (lon, lat).
  real(kind=wp), dimension(:,:), allocatable :: surface_temperature !< Surface temperature [K] (lon, lat).
  real(kind=wp), dimension(:), allocatable :: time !< Time [hours or days or months].
  character(len=128) :: time_units !< Time units (i.e. days since 0000-00-00 00:00:00)
  real(kind=wp) :: total_solar_irradiance !< Total solar irradiance [W m-2].
  real(kind=wp), dimension(:,:,:,:), allocatable :: aerosols !< Aerosol concentrations (lon, lat, layer, species).
end type Atmosphere_t


public :: create_atmosphere
public :: read_time_slice
public :: destroy_atmosphere
integer, parameter, public :: h2o = 1
integer, parameter, public :: o3 = 2
character(len=64), dimension(46) :: valid_names = [ &
  "solar_constant             ", &
  "earth_sun_distance_fraction", &
  "surface_temperature        ", &
  "land_fraction              ", &
  "cosine_zenith              ", &
  "daylight_fraction          ", &
  "visible_direct_albedo      ", &
  "visible_diffuse_albedo     ", &
  "infrared_direct_albedo     ", &
  "infrared_diffuse_albedo    ", &
  "layer_pressure             ", &
  "layer_temperature          ", &
  "water_vapor                ", &
  "ozone                      ", &
  "level_pressure             ", &
  "level_temperature          ", &
  "layer_thickness            ", &
  "stratiform_cloud_fraction  ", &
  "shallow_cloud_fraction     ", &
  "stratiform_ice_content     ", &
  "shallow_ice_content        ", &
  "stratiform_liquid_content  ", &
  "shallow_liquid_content     ", &
  "stratiform_droplet_number  ", &
  "shallow_droplet_number     ", &
  "stratiform_ice_number      ", &
  "stratiform_rain            ", &
  "stratiform_rain_size       ", &
  "stratiform_snow            ", &
  "stratiform_snow_size       ", &
  "soa_concentration          ", &
  "dust1_concentration        ", &
  "dust2_concentration        ", &
  "dust3_concentration        ", &
  "dust4_concentration        ", &
  "dust5_concentration        ", &
  "sulfate_concentration      ", &
  "ssalt1_concentration       ", &
  "ssalt2_concentration       ", &
  "ssalt3_concentration       ", &
  "ssalt4_concentration       ", &
  "ssalt5_concentration       ", &
  "bcphobic_concentration     ", &
  "bcphilic_concentration     ", &
  "omphobic_concentration     ", &
  "omphilic_concentration     " &
]


character(len=256) :: atmos_path = "none"
logical :: cleansky = .false.
logical :: clearsky = .false.
integer, dimension(2) :: io_layout = (/1, 1/)
integer, dimension(2) :: layout = (/0, 0/)
character(len=256) :: override_path = "none"
character(len=64), dimension(20) :: override_variables = &
  ["none", "none", "none", "none", "none", "none", "none", "none", "none", "none", &
   "none", "none", "none", "none", "none", "none", "none", "none", "none", "none"]
integer :: override_z_lower = -1
integer :: override_z_upper = -1
namelist /am4_nml/ atmos_path, cleansky, clearsky, io_layout, layout, &
                   override_path, override_variables, override_z_lower, &
                   override_z_upper


interface override_data
  module procedure override_data_1d
  module procedure override_data_2d
  module procedure override_data_3d
end interface override_data


contains


!> @brief Reserve memory and read in atmospheric data.
subroutine create_atmosphere(atm, column_blocking, nxblocks, nyblocks)

  type(Atmosphere_t), dimension(:), allocatable, intent(inout) :: atm
  type(block_control_type), intent(inout) :: column_blocking
  integer, intent(in) :: nxblocks
  integer, intent(in) :: nyblocks

  character(len=256) :: attr
  type(FmsNetcdfDomainFile_t) :: dataset
  integer, dimension(4) :: dim_sizes
  logical :: block_flag
  integer :: err, i, iec, ierr, isc, j, jec, jsc, num_blocks, nx, ny
  type(domain2d), pointer :: io_domain
  type(FmsNetcdfFile_t) :: tilefile

  read(input_nml_file, am4_nml, iostat=ierr)
  err = check_nml_error(ierr, "am4_nml")

  !Sanity checks.
  if (trim(override_path) .ne. "none") then
    do j = 1, size(override_variables)
      if (trim(override_variables(j)) .ne. "none") then
        do i = 1, size(valid_names)
          if (trim(override_variables(j)) .eq. trim(valid_names(i))) then
            exit
          endif
        enddo
        if (i .gt. size(valid_names)) then
          call error_mesg("create_atmosphere", &
                          "cannot find override variable "//trim(override_variables(j))//".", fatal)
        endif
      endif
    enddo
  endif

  !Open the *.tile1.nc file and get the longitude and latitude sizes.
  i = len_trim(atmos_path)
  if (.not. open_file(tilefile, atmos_path(:i - 3)//".tile1.nc", "read")) then
    call error_mesg("create_atmosphere", &
                    "cannot find the "//trim(atmos_path)//" tile files.", fatal)
  endif
  call get_dimension_size(tilefile, "grid_xt", nx)
  call get_dimension_size(tilefile, "grid_yt", ny)
  call close_file(tilefile)

  !Create a 2d domain and define the column blocking.
  num_blocks = nxblocks*nyblocks
  allocate(atm(num_blocks))
  call create_atmosphere_domain(nx, ny, atm(1)%domain)
  call mpp_get_compute_domain(atm(1)%domain, xbegin=isc, xend=iec, ybegin=jsc, yend=jec)
  block_flag = .true.
  call define_blocks("atmosphere", column_blocking, isc, iec, jsc, jec, 1, &
                     nxblocks, nyblocks, block_flag)

  !Get the longitude and latitude bounds for each domain grid cell.
  call lon_and_lat_bounds(atm(1))

  !Open dataset.
  if (.not. open_file(dataset, atmos_path, "read", atm(1)%domain)) then
    call error_mesg("create_atmosphere", "cannot open file "//trim(atmos_path)//".", &
                    fatal)
  endif

  !Get domain sizes.
  io_domain => mpp_get_io_domain(atm(1)%domain)
  if (.not. associated(io_domain)) then
    call error_mesg("create_atmosphere", "I/O domain doesn't exist.", fatal)
  endif
  call mpp_get_compute_domain(io_domain, xsize=nx, ysize=ny)

  !Identify longitude and latitude dimensions.
  call register_axis(dataset, "grid_xt", "x")
  call register_axis(dataset, "grid_yt", "y")

  !Read in axes so they can be passed to diag_manager.
  call get_dimension_size(dataset, "grid_xt", i)
  allocate(atm(1)%longitude(i))
  call read_data(dataset, "grid_xt", atm(1)%longitude)
  call get_dimension_size(dataset, "grid_yt", i)
  allocate(atm(1)%latitude(i))
  call read_data(dataset, "grid_yt", atm(1)%latitude)
  call get_dimension_size(dataset, "pfull", i)
  allocate(atm(1)%layer(i))
  call read_data(dataset, "pfull", atm(1)%layer)
  call get_dimension_size(dataset, "phalf", i)
  allocate(atm(1)%level(i))
  call read_data(dataset, "phalf", atm(1)%level)

  !Read in the times.
  call get_variable_size(dataset, "time", dim_sizes(1:1))
  atm(1)%num_times = dim_sizes(1)
  allocate(atm(1)%time(atm(1)%num_times))
  call read_data(dataset, "time", atm(1)%time)
  call get_variable_attribute(dataset, "time", "units", atm(1)%time_units)
  call get_variable_attribute(dataset, "time", "calendar", atm(1)%calendar)

  !Get the number of layers and levels.
  call get_dimension_size(dataset, "pfull", atm(1)%num_layers)
  atm(1)%num_levels = atm(1)%num_layers + 1
  call close_file(dataset)

  !Allocate the input data buffers.
  do i = 1, num_blocks
    nx = column_blocking%ibe(i) - column_blocking%ibs(i) + 1
    ny = column_blocking%jbe(i) - column_blocking%jbs(i) + 1
    allocate(atm(i)%surface_temperature(nx, ny))
    allocate(atm(i)%land_fraction(nx, ny))
    allocate(atm(i)%solar_zenith_angle(nx, ny))
    allocate(atm(i)%daylight_fraction(nx, ny))
    allocate(atm(i)%layer_pressure(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%level_pressure(nx, ny, atm(1)%num_levels))
    allocate(atm(i)%layer_temperature(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%level_temperature(nx, ny, atm(1)%num_levels))
    allocate(atm(i)%ppmv(nx, ny, atm(1)%num_layers, 2))
    allocate(atm(i)%surface_albedo_direct_uv(nx, ny))
    allocate(atm(i)%surface_albedo_diffuse_uv(nx, ny))
    allocate(atm(i)%surface_albedo_direct_ir(nx, ny))
    allocate(atm(i)%surface_albedo_diffuse_ir(nx, ny))
    allocate(atm(i)%layer_thickness(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%stratiform_cloud_fraction(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%shallow_cloud_fraction(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%stratiform_cloud_ice_content(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%shallow_cloud_ice_content(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%stratiform_cloud_liquid_content(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%shallow_cloud_liquid_content(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%stratiform_droplet_number(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%shallow_droplet_number(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%stratiform_ice_number(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%stratiform_rain(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%stratiform_rain_size(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%stratiform_snow(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%stratiform_snow_size(nx, ny, atm(1)%num_layers))
    allocate(atm(i)%aerosols(nx, ny, atm(1)%num_layers, 16))
  enddo
end subroutine create_atmosphere


subroutine read_time_slice(atm, time_level, column_blocking)

  type(Atmosphere_t), dimension(:), intent(inout) :: atm
  integer, intent(in) :: time_level
  type(block_control_type), intent(inout) :: column_blocking

  real, dimension(1, 1) :: buffer_scalar
  real, dimension(:, :), allocatable :: buffer2d
  real, dimension(:, :, :), allocatable :: buffer3d_lay
  real, dimension(:, :, :), allocatable :: buffer3d_lev
  type(FmsNetcdfDomainFile_t) :: dataset
  integer :: i, num_blocks, nx, ny
  type(domain2d), pointer :: io_domain

  !Open dataset.
  if (.not. open_file(dataset, atmos_path, "read", atm(1)%domain)) then
    call error_mesg("create_atmosphere", "cannot open file "//trim(atmos_path)//".", &
                    fatal)
  endif

  !Identify longitude and latitude dimensions.
  call register_axis(dataset, "grid_xt", "x")
  call register_axis(dataset, "grid_yt", "y")

  !Allocate buffers.
  io_domain => mpp_get_io_domain(atm(1)%domain)
  call mpp_get_compute_domain(io_domain, xsize=nx, ysize=ny)
  allocate(buffer2d(nx, ny))
  allocate(buffer3d_lay(nx, ny, atm(1)%num_layers))
  allocate(buffer3d_lev(nx, ny, atm(1)%num_levels))
  num_blocks = size(column_blocking%ibs, 1)

  !Read in the data and block it.
  call read_data(dataset, "solar_constant", buffer_scalar(:, 1), unlim_dim_level=time_level)
  call override_data("solar_constant", buffer_scalar(:, 1), time_level, atm(1)%domain)
  do i = 1, num_blocks
    atm(i)%total_solar_irradiance = buffer_scalar(1, 1)
  enddo

  call read_data(dataset, "earth_sun_distance_fraction", buffer_scalar(:, 1), unlim_dim_level=time_level)
  call override_data("earth_sun_distance_fraction", buffer_scalar(:, 1), time_level, atm(1)%domain)
  do i = 1, num_blocks
    atm(i)%earth_sun_distance_fraction = buffer_scalar(1, 1)
  enddo

  call read_data(dataset, "surface_temperature", buffer2d, unlim_dim_level=time_level)
  call override_data("surface_temperature", buffer2d, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_2d(buffer2d, atm(i)%surface_temperature, column_blocking, i)
  enddo

  call read_data(dataset, "land_fraction", buffer2d, unlim_dim_level=time_level)
  call override_data("land_fraction", buffer2d, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_2d(buffer2d, atm(i)%land_fraction, column_blocking, i)
  enddo

  call read_data(dataset, "cosine_zenith", buffer2d, unlim_dim_level=time_level)
  call override_data("cosine_zenith", buffer2d, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_2d(buffer2d, atm(i)%solar_zenith_angle, column_blocking, i)
  enddo

  call read_data(dataset, "daylight_fraction", buffer2d, unlim_dim_level=time_level)
  call override_data("daylight_fraction", buffer2d, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_2d(buffer2d, atm(i)%daylight_fraction, column_blocking, i)
  enddo

  call read_data(dataset, "visible_direct_albedo", buffer2d, unlim_dim_level=time_level)
  call override_data("visible_direct_albedo", buffer2d, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_2d(buffer2d, atm(i)%surface_albedo_direct_uv, column_blocking, i)
  enddo

  call read_data(dataset, "visible_diffuse_albedo", buffer2d, unlim_dim_level=time_level)
  call override_data("visible_diffuse_albedo", buffer2d, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_2d(buffer2d, atm(i)%surface_albedo_diffuse_uv, column_blocking, i)
  enddo

  call read_data(dataset, "infrared_direct_albedo", buffer2d, unlim_dim_level=time_level)
  call override_data("infrared_direct_albedo", buffer2d, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_2d(buffer2d, atm(i)%surface_albedo_direct_ir, column_blocking, i)
  enddo

  call read_data(dataset, "infrared_diffuse_albedo", buffer2d, unlim_dim_level=time_level)
  call override_data("infrared_diffuse_albedo", buffer2d, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_2d(buffer2d, atm(i)%surface_albedo_diffuse_ir, column_blocking, i)
  enddo

  call read_data(dataset, "layer_pressure", buffer3d_lay, unlim_dim_level=time_level)
  call override_data("layer_pressure", buffer3d_lay, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_3d(buffer3d_lay, atm(i)%layer_pressure, column_blocking, i)
  enddo

  call read_data(dataset, "layer_temperature", buffer3d_lay, unlim_dim_level=time_level)
  call override_data("layer_temperature", buffer3d_lay, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_3d(buffer3d_lay, atm(i)%layer_temperature, column_blocking, i)
  enddo

  call read_data(dataset, "water_vapor", buffer3d_lay, unlim_dim_level=time_level)
  call override_data("water_vapor", buffer3d_lay, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_3d(buffer3d_lay, atm(i)%ppmv(:, :, :, h2o), column_blocking, i)
  enddo

  call read_data(dataset, "ozone", buffer3d_lay, unlim_dim_level=time_level)
  call override_data("ozone", buffer3d_lay, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_3d(buffer3d_lay, atm(i)%ppmv(:, :, :, o3), column_blocking, i)
  enddo

  call read_data(dataset, "level_pressure", buffer3d_lev, unlim_dim_level=time_level)
  call override_data("level_pressure", buffer3d_lev, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_3d(buffer3d_lev, atm(i)%level_pressure, column_blocking, i)
  enddo

  call read_data(dataset, "level_temperature", buffer3d_lev, unlim_dim_level=time_level)
  call override_data("level_temperature", buffer3d_lev, time_level, atm(1)%domain)
  do i = 1, num_blocks
    call block_data_3d(buffer3d_lev, atm(i)%level_temperature, column_blocking, i)
  enddo

  if (clearsky) then
    do i = 1, num_blocks
      atm(i)%layer_thickness = 0.
      atm(i)%stratiform_cloud_fraction = 0.
      atm(i)%shallow_cloud_fraction = 0.
      atm(i)%stratiform_cloud_ice_content = 0.
      atm(i)%shallow_cloud_ice_content = 0.
      atm(i)%stratiform_cloud_liquid_content = 0.
      atm(i)%shallow_cloud_liquid_content = 0.
      atm(i)%stratiform_droplet_number = 0.
      atm(i)%shallow_droplet_number = 0.
      atm(i)%stratiform_ice_number = 0.
      atm(i)%stratiform_rain = 0.
      atm(i)%stratiform_rain_size = 0.
      atm(i)%stratiform_snow = 0.
      atm(i)%stratiform_snow_size = 0.
    enddo
  else
    call read_data(dataset, "layer_thickness", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("layer_thickness", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%layer_thickness, column_blocking, i)
    enddo

    call read_data(dataset, "stratiform_cloud_fraction", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("stratiform_cloud_fraction", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%stratiform_cloud_fraction, column_blocking, i)
    enddo

    call read_data(dataset, "shallow_cloud_fraction", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("shallow_cloud_fraction", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%shallow_cloud_fraction, column_blocking, i)
    enddo

    call read_data(dataset, "stratiform_ice_content", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("stratiform_ice_content", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%stratiform_cloud_ice_content, column_blocking, i)
    enddo

    call read_data(dataset, "shallow_ice_content", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("shallow_ice_content", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%shallow_cloud_ice_content, column_blocking, i)
    enddo

    call read_data(dataset, "stratiform_liquid_content", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("stratiform_liquid_content", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%stratiform_cloud_liquid_content, column_blocking, i)
    enddo

    call read_data(dataset, "shallow_liquid_content", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("shallow_liquid_content", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%shallow_cloud_liquid_content, column_blocking, i)
    enddo

    call read_data(dataset, "stratiform_droplet_number", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("stratiform_droplet_number", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%stratiform_droplet_number, column_blocking, i)
    enddo

    call read_data(dataset, "shallow_droplet_number", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("shallow_droplet_number", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%shallow_droplet_number, column_blocking, i)
    enddo

    call read_data(dataset, "stratiform_ice_number", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("stratiform_ice_number", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%stratiform_ice_number, column_blocking, i)
    enddo

    call read_data(dataset, "stratiform_rain", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("stratiform_rain", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%stratiform_rain, column_blocking, i)
    enddo

    call read_data(dataset, "stratiform_rain_size", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("stratiform_rain_size", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%stratiform_rain_size, column_blocking, i)
    enddo

    call read_data(dataset, "stratiform_snow", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("stratiform_snow", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%stratiform_snow, column_blocking, i)
    enddo

    call read_data(dataset, "stratiform_snow_size", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("stratiform_snow_size", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%stratiform_snow_size, column_blocking, i)
    enddo
  endif

  if (cleansky) then
    do i = 1, num_blocks
      atm(i)%aerosols = 0.
    enddo
  else
    call read_data(dataset, "soa_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("soa_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 1), column_blocking, i)
    enddo

    call read_data(dataset, "dust1_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("dust1_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 2), column_blocking, i)
    enddo

    call read_data(dataset, "dust2_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("dust2_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 3), column_blocking, i)
    enddo

    call read_data(dataset, "dust3_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("dust3_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 4), column_blocking, i)
    enddo

    call read_data(dataset, "dust4_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("dust4_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 5), column_blocking, i)
    enddo

    call read_data(dataset, "dust5_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("dust5_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 6), column_blocking, i)
    enddo

    call read_data(dataset, "sulfate_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("sulfate_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 7), column_blocking, i)
    enddo

    call read_data(dataset, "ssalt1_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("ssalt1_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 8), column_blocking, i)
    enddo

    call read_data(dataset, "ssalt2_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("ssalt2_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 9), column_blocking, i)
    enddo

    call read_data(dataset, "ssalt3_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("ssalt3_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 10), column_blocking, i)
    enddo

    call read_data(dataset, "ssalt4_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("ssalt4_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 11), column_blocking, i)
    enddo

    call read_data(dataset, "ssalt5_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("ssalt5_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 12), column_blocking, i)
    enddo

    call read_data(dataset, "bcphobic_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("bcphobic_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 13), column_blocking, i)
    enddo

    call read_data(dataset, "bcphilic_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("bcphilic_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 14), column_blocking, i)
    enddo

    call read_data(dataset, "omphobic_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("omphobic_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 15), column_blocking, i)
    enddo

    call read_data(dataset, "omphilic_concentration", buffer3d_lay, unlim_dim_level=time_level)
    call override_data("omphilic_concentration", buffer3d_lay, time_level, atm(1)%domain)
    do i = 1, num_blocks
      call block_data_3d(buffer3d_lay, atm(i)%aerosols(:, :, :, 16), column_blocking, i)
    enddo
  endif
  call close_file(dataset)
  deallocate(buffer2d, buffer3d_lay, buffer3d_lev)
end subroutine read_time_slice


subroutine override_data_1d(name, buffer, time_level, domain)

  character(len=*), intent(in) :: name
  real, dimension(:), intent(inout) :: buffer
  integer, intent(in) :: time_level
  type(domain2d), intent(in) :: domain !< 2d domain.

  type(FmsNetcdfDomainFile_t) :: dataset
  integer :: i

  do i = 1, size(override_variables)
    if (trim(override_variables(i)) .eq. "none") then
      !Return once the first "none" is found.
       return
    endif
    if (trim(override_variables(i)) .eq. trim(name)) then
      exit
    endif
  enddo
  if (i .gt. size(override_variables)) then
    return
  endif
  if (.not. open_file(dataset, override_path, "read", domain)) then
    call error_mesg("override_data_1d", "cannot find "//trim(override_path)//".", fatal)
  endif
  call register_axis(dataset, "grid_xt", "x")
  call register_axis(dataset, "grid_yt", "y")
  call read_data(dataset, name, buffer(:), unlim_dim_level=time_level)
  call close_file(dataset)
end subroutine override_data_1d


subroutine override_data_2d(name, buffer, time_level, domain)

  character(len=*), intent(in) :: name
  real, dimension(:, :), intent(inout) :: buffer
  integer, intent(in) :: time_level
  type(domain2d), intent(in) :: domain !< 2d domain.

  type(FmsNetcdfDomainFile_t) :: dataset
  integer :: i

  do i = 1, size(override_variables)
    if (trim(override_variables(i)) .eq. "none") then
      !Return once the first "none" is found.
       return
    endif
    if (trim(override_variables(i)) .eq. trim(name)) then
      exit
    endif
  enddo
  if (i .gt. size(override_variables)) then
    return
  endif
  if (.not. open_file(dataset, override_path, "read", domain)) then
    call error_mesg("override_data_2d", "cannot find "//trim(override_path)//".", fatal)
  endif
  call register_axis(dataset, "grid_xt", "x")
  call register_axis(dataset, "grid_yt", "y")
  call read_data(dataset, name, buffer, unlim_dim_level=time_level)
  call close_file(dataset)
end subroutine override_data_2d


subroutine override_data_3d(name, buffer, time_level, domain)

  character(len=*), intent(in) :: name
  real, dimension(:, :, :), intent(inout) :: buffer
  integer, intent(in) :: time_level
  type(domain2d), intent(in) :: domain !< 2d domain.

  real, dimension(:, :, :), allocatable :: override_buffer
  integer, dimension(3) :: count_
  type(FmsNetcdfDomainFile_t) :: dataset
  integer :: i
  integer, dimension(3) :: start

  do i = 1, size(override_variables)
    if (trim(override_variables(i)) .eq. "none") then
      !Return once the first "none" is found.
       return
    endif
    if (trim(override_variables(i)) .eq. trim(name)) then
      exit
    endif
  enddo
  if (i .gt. size(override_variables)) then
    return
  endif
  if (override_z_lower .lt. 0 .or. override_z_upper .lt. 0) then
    call error_mesg("override_data_3d", "you must set the overrize z limits.", fatal)
  endif
  if (.not. open_file(dataset, override_path, "read", domain)) then
    call error_mesg("override_data_3d", "cannot find "//trim(override_path)//".", fatal)
  endif
  call register_axis(dataset, "grid_xt", "x")
  call register_axis(dataset, "grid_yt", "y")
  allocate(override_buffer(size(buffer, 1), size(buffer, 2), size(buffer, 3)))
  call read_data(dataset, name, override_buffer, unlim_dim_level=time_level)
  call close_file(dataset)
  buffer(:, :, override_z_lower:override_z_upper) = override_buffer(:, :, override_z_lower:override_z_upper)
  deallocate(override_buffer)
end subroutine override_data_3d


subroutine block_data_2d(src, dest, column_blocking, block_)

  real, dimension(:, :), intent(in) :: src
  real, dimension(:, :), intent(inout) :: dest
  type(block_control_type), intent(in) :: column_blocking
  integer, intent(in) :: block_

  integer :: is, ie, js, je, num_blocks

  is = column_blocking%ibs(block_) - column_blocking%isc + 1
  ie = column_blocking%ibe(block_) - column_blocking%isc + 1
  js = column_blocking%jbs(block_) - column_blocking%jsc + 1
  je = column_blocking%jbe(block_) - column_blocking%jsc + 1
  dest(:, :) = src(is:ie, js:je)
end subroutine block_data_2d


subroutine block_data_3d(src, dest, column_blocking, block_)

  real, dimension(:, :, :), intent(in) :: src
  real, dimension(:, :, :), intent(inout) :: dest
  type(block_control_type), intent(in) :: column_blocking
  integer, intent(in) :: block_

  integer :: is, ie, js, je

  is = column_blocking%ibs(block_) - column_blocking%isc + 1
  ie = column_blocking%ibe(block_) - column_blocking%isc + 1
  js = column_blocking%jbs(block_) - column_blocking%jsc + 1
  je = column_blocking%jbe(block_) - column_blocking%jsc + 1
  dest(:, :, :) = src(is:ie, js:je, :)
end subroutine block_data_3d


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
  if (allocated(atm%latitude_bounds)) deallocate(atm%latitude_bounds)
  if (allocated(atm%layer_pressure)) deallocate(atm%layer_pressure)
  if (allocated(atm%layer_temperature)) deallocate(atm%layer_temperature)
  if (allocated(atm%layer_thickness)) deallocate(atm%layer_thickness)
  if (allocated(atm%level)) deallocate(atm%level)
  if (allocated(atm%level_pressure)) deallocate(atm%level_pressure)
  if (allocated(atm%level_temperature)) deallocate(atm%level_temperature)
  if (allocated(atm%longitude)) deallocate(atm%longitude)
  if (allocated(atm%longitude_bounds)) deallocate(atm%longitude_bounds)
  if (allocated(atm%ppmv)) deallocate(atm%ppmv)
  if (allocated(atm%solar_zenith_angle)) deallocate(atm%solar_zenith_angle)
  if (allocated(atm%surface_albedo_diffuse_ir)) deallocate(atm%surface_albedo_diffuse_ir)
  if (allocated(atm%surface_albedo_diffuse_uv)) deallocate(atm%surface_albedo_diffuse_uv)
  if (allocated(atm%surface_albedo_direct_ir)) deallocate(atm%surface_albedo_direct_ir)
  if (allocated(atm%surface_albedo_direct_uv)) deallocate(atm%surface_albedo_direct_uv)
  if (allocated(atm%surface_temperature)) deallocate(atm%surface_temperature)
  if (allocated(atm%time)) deallocate(atm%time)
  if (allocated(atm%stratiform_droplet_number)) deallocate(atm%stratiform_droplet_number)
  if (allocated(atm%stratiform_ice_number)) deallocate(atm%stratiform_ice_number)
  if (allocated(atm%stratiform_rain)) deallocate(atm%stratiform_rain)
  if (allocated(atm%stratiform_rain_size)) deallocate(atm%stratiform_rain_size)
  if (allocated(atm%stratiform_snow)) deallocate(atm%stratiform_snow)
  if (allocated(atm%stratiform_snow_size)) deallocate(atm%stratiform_snow_size)
  if (allocated(atm%shallow_droplet_number)) deallocate(atm%shallow_droplet_number)
  if (allocated(atm%land_fraction)) deallocate(atm%land_fraction)
  if (allocated(atm%aerosols)) deallocate(atm%aerosols)
end subroutine destroy_atmosphere


!> @brief Create a 2d domain.
subroutine create_atmosphere_domain(nx, ny, domain)

  integer, intent(in) :: nx !< Number of points per tile in the x-direction.
  integer, intent(in) :: ny !< Number of points per tile in the y-direction.
  type(domain2d), intent(inout) :: domain !< A cubed-sphere domain.

  integer, dimension(4,6) :: global_indices
  integer, dimension(2,6) :: layout_
  integer, dimension(6) :: ni, nj, pe_start, pe_end
  integer, dimension(2) :: msize
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
    layout_(:,i) = layout(:)
    if (layout_(1,i) .gt. nx .or. layout_(2,i) .gt. ny) then
      call error_mesg("create_atmosphere_domain", &
                      "layout is greater than dimension length.", fatal)
    endif
    pe_start(i) = (i - 1)*npes/ntiles
    pe_end(i) = i*npes/ntiles - 1
  enddo

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
  msize(1) = maxval(ni(:)/layout_(1,:)) + whalo + ehalo + 1
  msize(2) = maxval(nj(:)/layout_(2,:)) + shalo + nhalo + 1
  call mpp_define_mosaic(global_indices, layout_, domain, ntiles, num_contact, tile1, &
                         tile2, istart1, iend1, jstart1, jend1, istart2, iend2, &
                         jstart2, jend2, pe_start, pe_end, symmetry=.true., &
                         whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                         name=trim("Cubed-sphere"), memory_size=msize)
  call mpp_define_io_domain(domain, io_layout)
end subroutine create_atmosphere_domain


subroutine lon_and_lat_bounds(atm)

  type(Atmosphere_t), intent(inout) :: atm

  integer :: nx, ny
  integer, dimension(1) :: tile

  call mpp_get_compute_domain(atm%domain, xsize=nx, ysize=ny)
  if (mpp_get_current_ntile(atm%domain) .ne. size(tile)) then
    call error_mesg("lon_and_lat_bounds", "More than one domain tile per rank.", fatal)
  endif
  tile = mpp_get_tile_id(atm%domain)
  allocate(atm%longitude_bounds(nx + 1, ny + 1))
  allocate(atm%latitude_bounds(nx + 1, ny + 1))
  call get_grid_cell_vertices("ATM", tile(1), atm%longitude_bounds, &
                              atm%latitude_bounds, atm%domain)
  atm%longitude_bounds(:,:) = atm%longitude_bounds(:,:)*(2.*pi/360.)
  atm%latitude_bounds(:,:) = atm%latitude_bounds(:,:)*(2.*pi/360.)
end subroutine lon_and_lat_bounds


end module am4
