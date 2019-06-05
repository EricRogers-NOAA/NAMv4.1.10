module ensemble_fastread_commvar
!$$$ module documentation block
!           .      .    .                                       .
! module:   ensemble_fastread_commvar
!   prgmmr: parrish          org: np22                date: 2017-01-30
!
! abstract:  This module is used as temporary storage for gfs nems format ensemble fields.
!             This is a special application of parallel read for using gfs ensemble with
!             regional NAM model analysis.
!
! program history log:
!   2017-01-30  parrish - initial documentation
!
! subroutines included:
!   sub init_nmmb_to_a
!   sub nmmb_h_to_a
!   sub nmmb_h_to_a8
!   sub nmmb_v_to_a
!   sub nmmb_a_to_h
!   sub nmmb_a_to_v
!   sub b_to_a_interpolate
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block

!     temporary storage for ensemble fastread (modify so no oz, cw to save space, since oz, cw
!                                              not used for ensembles in this application)


   use kinds, only: r_kind,r_single,i_kind

   implicit none

! set default to private
   private
! set subroutines to public
   public :: init_ens_fast_read
! set passed variables to public
   public :: ens_fast_read
   public :: en_loc3,clons,slons,m_cvars2d,m_cvars3d
   public :: n2d,n3d,vars2d,vars3d,filenames
   public :: lat1_gfs,lon1_gfs,lat2_gfs,lon2_gfs
   public :: istart_gfs,jstart_gfs
   public :: nlon_gfs,nlat_gfs,nsig_gfs,num_fields_gfs

   real(r_kind),allocatable,dimension(:):: clons(:),slons(:)  !  temporary arrays for later computation
                                                              !    of u,v pole rows
   real(r_single),allocatable,dimension(:,:,:,:,:):: en_loc3  !  temporary array to hold ensembles
   integer(i_kind),allocatable:: istart_gfs(:),jstart_gfs(:)
   integer(i_kind),allocatable :: m_cvars2d(:),m_cvars3d(:)   !  map to vars in en_loc3

   integer(i_kind) :: ens_fast_read     ! =0   -- slow read
                                        ! =-1  -- fast read in get_gefs_from_???
                                        ! = 1  -- fast read in gsi_mod (uses less memory ??)
   integer(i_kind) :: n2d                                     !  number 2d vars (psfc and zsfc)
   integer(i_kind) :: n3d                                     !  number 3d vars (u,v,tv,q)
   character(len=4),allocatable :: vars2d(:),vars3d(:)        !  var names
   character(len=255),allocatable :: filenames(:)             !  ensemble file names
   integer(i_kind) lat1_gfs,lon1_gfs,lat2_gfs,lon2_gfs        !  gfs subdomain dims
   integer(i_kind) nlon_gfs,nlat_gfs,nsig_gfs,num_fields_gfs  !  gfs global dims

contains

subroutine init_ens_fast_read 
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    init_ens_fast_read  initialize ens_fast_read = 0 (slow ensemble read)
!   prgmmr: parrish          org: np22                date: 2017-01-30
!
! abstract: Initialize value of input parameter ens_fast_read = 0 (slow ensemble read)
!
! program history log:
!   2017-01-30  parrish, initial documentation
!
!   input argument list:
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

   implicit none

   ens_fast_read=0

end subroutine init_ens_fast_read

end module ensemble_fastread_commvar
