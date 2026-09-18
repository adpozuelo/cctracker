! ==============================================================================
! PROGRAM: centers (CPU)
! ==============================================================================
! DESCRIPTION:
!   Cluster center trajectory tracking and lifetime analysis for molecular
!   dynamics and Monte Carlo simulations (CPU implementation).
!   Reads LAMMPS-formatted cluster center coordinates, calculates the minimum
!   inter-center distance across time frames under periodic boundary conditions
!   (minimum image convention), links cluster trajectories across consecutive
!   snapshots, isolates All-Time-Living (ATL) persistent clusters, and computes
!   cluster lifetime distributions.
!
! FEATURES:
!   - Parses simulation parameters and LAMMPS trajectory dumps (centers.lammpstrj)
!   - Computes minimum pairwise distances under 3D periodic boundary conditions
!   - Adaptive distance thresholding for robust temporal cluster linking
!   - Continuous cluster trajectory reconstruction across frames
!   - Extraction of All-Time-Living (ATL) clusters to LAMMPS trajectory format
!   - Binned histogram generation of cluster lifetime survival times
!
! AUTHORS:
!   Antonio Díaz Pozuelo   - adiaz@iqf.csic.es
!   Enrique Lomba          - enrique.lomba@csic.es
!   Instituto de Química Física Blas Cabrera (IQF-CSIC)
!
! DATE:
!   September 2026
! ==============================================================================
program centers
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use mod_common
   use mod_io
   use mod_distances
   use mod_trajectories
   implicit none
   integer :: i, istat
   real :: distance_adj
   character(:), allocatable :: directory

   directory = '.'
   input_filename = directory//'/'//'Input'
   centers_filename = directory//'/'//'centers.lammpstrj'
   trajectories_filename = directory//'/'//'centers_atl.lammpstrj'
   histogram_filename = directory//'/'//'centers_histogram.dat'
   call read_input_file()
   call read_centers_file()
   print *, 'time_steps = ', time_steps
   print *, 'side = ', side
   print *, 'min_n_centers = ', min_n_centers
   print *, 'max_n_centers = ', max_n_centers
   call cpu_centers_min_distance()
   print *, 'cal_min_distance = ', min_distance
   distance_adj = 0.5
   min_distance = min_distance*distance_adj
   print *, 'distance_corrector = ', distance_adj
   print *, 'est_min_distance = ', min_distance
   call cpu_trajectories_analysis()
   print *, 'clusters_processed = ', max_n_centers_final
   call write_trajectories_atl()
   print *, 'n_clusters_atl = ', n_clusters_atl
   write (*, '(A,f7.2)') '[% cluster_atl] (cluster_atl / clusters_proc) = ', &
      (n_clusters_atl/real(max_n_centers_final))*100
   call write_histogram()

end program centers
