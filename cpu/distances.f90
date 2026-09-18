! ==============================================================================
! MODULE: mod_distances (CPU)
! ==============================================================================
! DESCRIPTION:
!   Spatial distance calculation routines for cluster centers under 3D periodic
!   boundary conditions (minimum image convention).
!
! AUTHORS:
!   Antonio Díaz Pozuelo   - adiaz@iqf.csic.es
!   Enrique Lomba          - enrique.lomba@csic.es
!   Instituto de Química Física Blas Cabrera (IQF-CSIC)
!
! DATE:
!   September 2026
! ==============================================================================
module mod_distances
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use mod_common
   implicit none
contains

! ==============================================================================
! SUBROUTINE: cpu_centers_min_distance
! ==============================================================================
! DESCRIPTION:
!   Calculates the global minimum distance between any pair of cluster centers
!   across all time steps in the simulation, applying periodic boundary conditions
!   (minimum image convention). The resulting minimum distance serves as a baseline
!   scale to establish the spatial threshold for tracking cluster identities across
!   consecutive time frames.
!
! ARGUMENTS:
!   None. Operates on shared module variables from mod_common:
!     - time_steps: Total number of time steps (snapshots)
!     - host_configurations: Array of configuration types containing center positions
!     - side(n_dim), side_div_2(n_dim): Simulation box dimensions for periodic wrapping
!     - min_distance: Global minimum distance output variable (updated in-place)
!
! METHOD:
!   Iterates through each time frame and loops over all unique pairs of cluster
!   centers (j < k). For each Cartesian coordinate (x, y, z), applies minimum
!   image convention using side and side/2. Evaluates the Euclidean distance
!   and updates the running minimum across all frames.
!
! AUTHORS:
!   Antonio Díaz Pozuelo   - adiaz@iqf.csic.es
!   Enrique Lomba          - enrique.lomba@csic.es
!   Instituto de Química Física Blas Cabrera (IQF-CSIC)
!
! DATE:
!   September 2026
! ==============================================================================
   subroutine cpu_centers_min_distance()
      integer :: i, j, k, s
      real(sp) :: distance
      real(sp) :: distance_arr(n_dim)
      print *, 'Proccesing minimal distances...'
      do i = 1, time_steps
         do j = 1, host_configurations(i)%n_centers - 1
            do k = j + 1, host_configurations(i)%n_centers
               do s = 1, n_dim
                  distance_arr(s) = host_configurations(i)%centers(s, j) - host_configurations(i)%centers(s, k)
                  if (distance_arr(s) > side_div_2(s)) distance_arr(s) = distance_arr(s) - side(s)
                  if (distance_arr(s) < -side_div_2(s)) distance_arr(s) = distance_arr(s) + side(s)
                  distance_arr(s) = distance_arr(s)*distance_arr(s)
               end do
               distance = sqrt(sum(distance_arr))
               if (distance < min_distance) min_distance = distance
            end do
         end do
      end do
   end subroutine cpu_centers_min_distance
end module mod_distances
