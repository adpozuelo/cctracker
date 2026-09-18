! ==============================================================================
! MODULE: mod_common (CPU)
! ==============================================================================
! DESCRIPTION:
!   Global definitions, data types, and shared state variables for CPU cluster
!   center analysis.
!   Defines floating-point precision constants (sp, dp), simulation box geometry
!   parameters with periodic boundary boundaries (side, side_div_2, low, high),
!   derived types for per-timestep snapshot configurations and tracked cluster
!   centers, and dynamically allocated global arrays for lifetimes and histograms.
!
! DATA STRUCTURES:
!   - configuration: Stores snapshot ID, number of centers, and coordinate/velocity
!                    arrays for all cluster centers in a given time frame.
!   - center: Represents a single cluster center point with 3D coordinates and
!             velocity vectors.
!
! AUTHORS:
!   Antonio Díaz Pozuelo   - adiaz@iqf.csic.es
!   Enrique Lomba          - enrique.lomba@csic.es
!   Instituto de Química Física Blas Cabrera (IQF-CSIC)
!
! DATE:
!   September 2026
! ==============================================================================
module mod_common
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   implicit none
   character(:), allocatable :: input_filename
   character(:), allocatable :: centers_filename
   character(:), allocatable :: trajectories_filename
   character(:), allocatable :: histogram_filename
   integer, parameter :: n_dim = 3
   integer :: time_steps, max_3n_centers, max_n_centers_final, n_clusters_atl
   integer :: min_n_centers = 1000000, max_n_centers = 0
   integer, allocatable :: clusters_time_life(:, :)
   integer, allocatable :: clusters_time_steps(:)
   integer, allocatable :: histogram(:)
   real(sp), parameter :: histogram_dt = 5
   real(sp) :: min_distance = 1000000_sp
   real(sp) :: side(n_dim), side_div_2(n_dim),high(n_dim),low(n_dim)
   real(sp), allocatable :: host_distances(:)
   type :: configuration
      integer :: id
      integer :: n_centers
      real(sp), allocatable :: centers(:, :)
      real(sp), allocatable :: vel(:, :)
   end type
   type(configuration), allocatable :: host_configurations(:)
   type :: center
      real(sp) :: coordinate(n_dim) = 0.0_sp
      real(sp) :: velocity(n_dim) = 0.0_sp
   end type
   type(center), allocatable :: host_trajectories(:, :)
   save
contains

end module mod_common
