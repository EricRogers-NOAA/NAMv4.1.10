subroutine get_user_ens_fastread_(ntindex,iret)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    get_user_ens_           pretend atmos bkg is the ensemble
!   prgmmr: mahajan          org: emc/ncep            date: 2016-06-30
!
! abstract: Read in GFS ensemble members in to GSI ensemble.  This is the
!            version which reads all ensemble members simultaneously in
!            parallel to n_ens processors.  This is followed by a scatter
!            to subdomains on all processors.  This version will only work
!            if n_ens <= npe, where npe is the total number of processors
!            available.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      NOTE:  In this version, just copy pole row values to halo rows beyond
!      pole.  This replicates what is done in curent GSI.  Postpone using proper
!      halo values beyond poles to later versions when genex is more extensively used.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! program history log:
!   2016-06-30  mahajan  - initial code
!   2016-10-11  parrish  - create fast parallel code
!
!   input argument list:
!     ntindex  - time index for ensemble
!     ens_atm_bundle - atm bundle w/ fields for ensemble
!
!   output argument list:
!     iret           - return code, 0 for successful read.
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

    use mpimod, only: mpi_comm_world,ierror,mype,npe,mpi_real8,mpi_integer4,mpi_max
    use kinds, only: i_kind,r_single,r_kind
    use constants, only: zero
    use hybrid_ensemble_parameters, only: n_ens
    use hybrid_ensemble_parameters, only: ensemble_path
    use ensemble_fastread_commvar, only: en_loc3
    use ensemble_fastread_commvar, only: vars2d,vars3d,n2d,n3d,m_cvars2d,m_cvars3d
    use ensemble_fastread_commvar, only: lat1_gfs,lon1_gfs,lat2_gfs,lon2_gfs
    use ensemble_fastread_commvar, only: istart_gfs,jstart_gfs,filenames
    use ensemble_fastread_commvar, only: nlon_gfs,nlat_gfs,nsig_gfs,num_fields_gfs
    use genex_mod, only: genex_info,genex_create_info,genex,genex_destroy_info

    implicit none

    ! Declare passed variables
    integer(i_kind),     intent(in   ) :: ntindex
    integer(i_kind),     intent(  out) :: iret

    ! Declare internal variables
    character(len=*),parameter :: myname='get_user_ens_fastread_'
    character(len=255) :: filename
    integer(i_kind) :: i,ii,j,jj,k,n
    integer(i_kind) :: io_pe,n_io_pe_s,n_io_pe_e,n_io_pe_em,i_ens
    integer(i_kind) :: ip,ips,ipe,jps,jpe
    integer(i_kind) :: ias,iae,iasm,iaem,iaemz,jas,jae,jasm,jaem,jaemz
    integer(i_kind) :: kas,kae,kasm,kaem,kaemz,mas,mae,masm,maem,maemz
    integer(i_kind) :: ibs,ibe,ibsm,ibem,ibemz,jbs,jbe,jbsm,jbem,jbemz
    integer(i_kind) :: kbs,kbe,kbsm,kbem,kbemz,mbs,mbe,mbsm,mbem,mbemz
    integer(i_kind) :: nlat,nlon,lat1,lon1,lat2,lon2,nsig
    type(genex_info) :: s_a2b
    integer(i_kind),allocatable :: istart(:),jstart(:)
    real(r_single),allocatable :: en_full(:,:,:,:)
    real(r_single),allocatable :: en_loc(:,:,:,:)
      integer(i_kind) nlonh
    integer(i_kind),allocatable :: m_cvars2dw(:),m_cvars3dw(:)

    lat1=lat1_gfs ; lon1=lon1_gfs ; lat2=lat2_gfs ; lon2=lon2_gfs
    nlat=nlat_gfs ; nlon=nlon_gfs ; nsig=nsig_gfs
    allocate(istart(npe),jstart(npe))
    istart=istart_gfs ; jstart=jstart_gfs

    iret = 0

    nlonh=nlon/2

!.......  preliminary code

!  set up partition of available processors for parallel read

      if(n_ens > npe) then
         write(6,*) ' CANNOT READ ENSEMBLE -- n_ens > npe, increase npe >= n_ens '
         iret= 99
         return
      endif
      call ens_io_partition(n_ens,io_pe,n_io_pe_s,n_io_pe_e,n_io_pe_em,npe,i_ens,mype)

! setup communicator for scatter to subdomains:

!      first, define gsi subdomain boundaries in global units:

       ip=1   !  halo width is hardwired at 1
       ips=istart(mype+1)
       ipe=ips+lat1-1
       jps=jstart(mype+1)
       jpe=jps+lon1-1

      ias=1 ; iae=0 ; jas=1 ; jae=0 ; kas=1 ; kae=0 ; mas=1 ; mae=0
      if(mype==io_pe) then
         ias=1 ; iae=nlat
         jas=1 ; jae=nlon
         kas=1 ; kae=n3d*nsig+n2d
         mas=n_io_pe_s ; mae=n_io_pe_em
      endif
      iasm=ias ; iaem=iae ; jasm=jas ; jaem=jae ; kasm=kas ; kaem=kae ; masm=mas ; maem=mae

      ibs =ips    ; ibe =ipe    ; jbs =jps    ; jbe =jpe
      ibsm=ibs-ip ; ibem=ibe+ip ; jbsm=jbs-ip ; jbem=jbe+ip
      kbs =1   ; kbe =n3d*nsig+n2d ; mbs =1   ; mbe =n_ens
      kbsm=kbs ; kbem=kbe ; mbsm=mbs ; mbem=mbe
      iaemz=max(iasm,iaem) ; jaemz=max(jasm,jaem)
      kaemz=max(kasm,kaem) ; maemz=max(masm,maem)
      ibemz=max(ibsm,ibem) ; jbemz=max(jbsm,jbem)
      kbemz=max(kbsm,kbem) ; mbemz=max(mbsm,mbem)
      call genex_create_info(s_a2b,ias ,iae ,jas ,jae ,kas ,kae ,mas ,mae , &
                                   ibs ,ibe ,jbs ,jbe ,kbs ,kbe ,mbs ,mbe , &
                                   iasm,iaem,jasm,jaem,kasm,kaem,masm,maem, &
                                   ibsm,ibem,jbsm,jbem,kbsm,kbem,mbsm,mbem)

!!  read ensembles
 
     allocate(en_full(iasm:iaemz,jasm:jaemz,kasm:kaemz,masm:maemz))

     filename='xxxxxx'
     if(mas==mae) filename=trim(adjustl(filenames(mas)))

!   allocate m_cvars2dw,m_cvars3dw

     allocate(m_cvars2dw(n2d))
     allocate(m_cvars3dw(n3d))

     m_cvars2dw=-999
     m_cvars3dw=-999

     call mpi_barrier(mpi_comm_world,ierror)
     if(mas==mae) &
        call parallel_read_nemsio_state(en_full,m_cvars2dw,m_cvars3dw, &
                                        ias ,jas ,mas ,mae  , &
                                        iasm,iaemz,jasm,jaemz,kasm,kaemz,masm,maemz, &
                                        filename,.true.)
     call mpi_barrier(mpi_comm_world,ierror)
                   ! if(1/=0) then
                   !    call mpi_finalize(ierror)
                   !    stop
                   ! endif

     call mpi_allreduce(m_cvars2dw,m_cvars2d,n2d,mpi_integer4,mpi_max,mpi_comm_world,ierror)
     call mpi_allreduce(m_cvars3dw,m_cvars3d,n3d,mpi_integer4,mpi_max,mpi_comm_world,ierror)
     deallocate(m_cvars2dw,m_cvars3dw)

! scatter to subdomains:

      allocate(en_loc(ibsm:ibemz,jbsm:jbemz,kbsm:kbemz,mbsm:mbemz))

          en_loc=zero
          call genex(s_a2b,en_full,en_loc)

      deallocate(en_full)
      call genex_destroy_info(s_a2b)  ! check on actual routine name

! transfer en_loc to en_loc3

      do n=1,n_ens
         do k=1,n2d+n3d*nsig
            jj=0
            do j=jbsm,jbem
               jj=jj+1
               ii=0
               do i=ibsm,ibem
                  ii=ii+1
                  en_loc3(ii,jj,k,n,ntindex)=en_loc(i,j,k,n)
               enddo
            enddo
         enddo
      enddo
      deallocate(istart,jstart,en_loc)

end subroutine get_user_ens_fastread_

subroutine move2bundle_(n,ntindex,atm_bundle,iret)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    move2bundle  transfer 1 ensemble member to bundle
!   prgmmr: mahajan          org: emc/ncep            date: 2016-06-30
!
! abstract: transfer one ensemble member to bundle
!
! program history log:
!   2016-06-30  parrish -- copy and adapt get_user_ens_member_ to transfer 1
!                            ensemble member
!
!   input argument list:
!     n          - ensemble member number
!     atm_bundle - empty atm bundle
!
!   output argument list:
!     atm_bundle - atm bundle w/ fields for ensemble member
!     iret       - return code
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

    use mpimod, only: mype,ierror
    use kinds, only: i_kind,r_kind,r_single
    use constants, only: zero,one,two,fv
    use gridmod, only: regional,use_gfs_nemsio
    use general_sub2grid_mod, only: sub2grid_info,general_sub2grid_destroy_info
    use hybrid_ensemble_parameters, only: uv_hyb_ens
    use hybrid_ensemble_parameters, only: sp_ens
    use hybrid_ensemble_isotropic, only: en_perts
    use hybrid_ensemble_parameters,only: n_ens,ntlevs_ens
    use gsi_bundlemod, only: gsi_bundle
    use gsi_bundlemod, only: gsi_bundlegetpointer,gsi_bundleputvar
    use gsi_bundlemod, only : assignment(=)
!   use control_vectors, only: cvars2d,cvars3d,nc2d,nc3d
    use ensemble_fastread_commvar,only: m_cvars2d,m_cvars3d,clons,slons,en_loc3
    use ensemble_fastread_commvar,only: n2d,n3d,vars2d,vars3d
    use ensemble_fastread_commvar,only: lat2_gfs,lon2_gfs,nlat_gfs,nlon_gfs,nsig_gfs,num_fields_gfs

    implicit none

    ! Declare passed variables
    integer(i_kind)     ,intent(in   ) :: n,ntindex
    type(gsi_bundle)    ,intent(inout) :: atm_bundle
    integer(i_kind),     intent(  out) :: iret

    ! Declare internal variables
    character(len=*),parameter :: myname='move2bundle_'
    character(len=70) :: filename

    integer(i_kind) :: ierr
    integer(i_kind) :: m,i,j,k
    real(r_kind),pointer,dimension(:,:)   :: ps
    real(r_kind),pointer,dimension(:,:)   :: z
  ! real(r_kind),pointer,dimension(:,:)   :: sst
    real(r_kind),pointer,dimension(:,:,:) :: u
    real(r_kind),pointer,dimension(:,:,:) :: v
    real(r_kind),pointer,dimension(:,:,:) :: tv
    real(r_kind),pointer,dimension(:,:,:) :: q
    type(sub2grid_info) :: grd2d,grd3d
    real(r_kind),parameter:: r0_001 = 0.001_r_kind

!    initialize atm_bundle to zero

       atm_bundle=zero

       call gsi_bundlegetpointer(atm_bundle,'ps',ps,  ierr); iret = ierr
       call gsi_bundlegetpointer(atm_bundle,'z',z,  ierr); iret = ierr + iret
     ! call gsi_bundlegetpointer(atm_bundle,'sst',sst, ierr); iret = ierr
       call gsi_bundlegetpointer(atm_bundle,'u',u ,  ierr); iret = ierr + iret
       call gsi_bundlegetpointer(atm_bundle,'v',v ,  ierr); iret = ierr + iret
       call gsi_bundlegetpointer(atm_bundle,'tv' ,tv,  ierr); iret = ierr + iret
       call gsi_bundlegetpointer(atm_bundle,'q' ,q ,  ierr); iret = ierr + iret
       if ( iret /= 0 ) then
          if ( mype == 0 ) then
             write(6,'(A)') trim(myname) // ': ERROR!'
             write(6,'(A)') trim(myname) // ': special, NAM requires : ps,u,v,(sf,vp)tv,q'
             write(6,'(A)') trim(myname) // ': but some have not been found. Aborting ... '
          endif
          goto 100
       endif

       do m=1,n2d
          if(trim(vars2d(m))=='ps') then
             do j=1,lon2_gfs
                do i=1,lat2_gfs
                   ps(i,j)=en_loc3(i,j,m_cvars2d(m),n,ntindex)
                enddo
             enddo
          endif
          if(trim(vars2d(m))=='z') then
             do j=1,lon2_gfs
                do i=1,lat2_gfs
                   z(i,j)=en_loc3(i,j,m_cvars2d(m),n,ntindex)
                enddo
             enddo
          endif
       enddo
       do m=1,n3d
          if(trim(vars3d(m))=='u') then
             do k=1,nsig_gfs
                do j=1,lon2_gfs
                   do i=1,lat2_gfs
                      u(i,j,k)=en_loc3(i,j,m_cvars3d(m)+k-1,n,ntindex)
                   enddo
                enddo
             enddo
          endif
          if(trim(vars3d(m))=='v') then
             do k=1,nsig_gfs
                do j=1,lon2_gfs
                   do i=1,lat2_gfs
                      v(i,j,k)=en_loc3(i,j,m_cvars3d(m)+k-1,n,ntindex)
                   enddo
                enddo
             enddo
          endif
          if(trim(vars3d(m))=='tv') then
             do k=1,nsig_gfs
                do j=1,lon2_gfs
                   do i=1,lat2_gfs
                      tv(i,j,k)=en_loc3(i,j,m_cvars3d(m)+k-1,n,ntindex)
                   enddo
                enddo
             enddo
          endif
          if(trim(vars3d(m))=='q') then
             do k=1,nsig_gfs
                do j=1,lon2_gfs
                   do i=1,lat2_gfs
                      q(i,j,k)=en_loc3(i,j,m_cvars3d(m)+k-1,n,ntindex)
                   enddo
                enddo
             enddo
          endif
       enddo

!   convert ps from Pa to cb
     ps=r0_001*ps
!   convert t to virtual temperature
     tv=tv*(one+fv*q)

!--- now update pole values of atm_bundle using general_sub2grid (so halos also
!       automatically updated.

       call create_grd23d(grd2d,nlat_gfs,nlon_gfs,1)
       call create_grd23d(grd3d,nlat_gfs,nlon_gfs,nsig_gfs)

       call update_scalar_poles(grd2d,z)
       call update_scalar_poles(grd2d,ps)
       call update_vector_poles(grd3d,u,v,clons,slons)
       call update_scalar_poles(grd3d,tv)
       call update_scalar_poles(grd3d,q)

       call gsi_bundleputvar(atm_bundle,'ps',ps,  ierr); iret = ierr
       call gsi_bundleputvar(atm_bundle,'z',z,  ierr); iret = ierr + iret
     ! call gsi_bundleputvar(atm_bundle,'sst',sst,ierr); iret = ierr + iret  ! no sst for now
       call gsi_bundleputvar(atm_bundle,'u',u ,  ierr); iret = ierr + iret
       call gsi_bundleputvar(atm_bundle,'v',v ,  ierr); iret = ierr + iret
       call gsi_bundleputvar(atm_bundle,'tv' ,tv,  ierr); iret = ierr + iret
       call gsi_bundleputvar(atm_bundle,'q' ,q ,  ierr); iret = ierr + iret
                  write(6,'(" after upd poles, for q,  n,ierr=",2i5)') n,ierr
                       write(6,'(" after upd poles, n,iret=",2i8)') n,iret

       call destroy_grd23d(grd2d)
                       write(6,'(" after destroy(grd2d),n=",i3)') n
       call destroy_grd23d(grd3d)
                       write(6,'(" after destroy(grd3d),n,iret=",2i5)') n,iret

       if ( iret /= 0 ) then
          if ( mype == 0 ) then
             write(6,'(A)') trim(myname) // ': ERROR!'
             write(6,'(A)') trim(myname) // ': For NAM, need following: ps,u,v,(sf,vp)tv,q'
             write(6,'(A)') trim(myname) // ': but some have not been found. Aborting ... '
          endif
          goto 100
       endif

100 continue
                       write(6,'(" after 100 continue,iret,n=",2i5)') iret,n

    if ( iret /= 0 ) then
       if ( mype == 0 ) then
          write(6,'(A)') trim(myname) // ': WARNING!'
          write(6,'(3A,I5)') trim(myname) // ': Trouble reading ensemble file : ', trim(filename), ', IRET = ', iret
       endif
    endif

!   if last ensemble member and last time level, remove temporary space en_loc3 used to store raw ensembles
                       write(6,'(" stop before call ensemble_postread,iret,n=",2i5)') iret,n
                     
                               ! if(n==2) then
                               !    call mpi_finalize(ierror)
                               !    stop
                               ! endif

    if(n==n_ens.and.ntindex==ntlevs_ens) call ensemble_postread
    if(n==n_ens.and.ntindex==ntlevs_ens) then
                                     write(6,*)' GOT AS FAR AS ENSEMBLE_POSTREAD '
                               ! if(1/=0) then
                               !    call mpi_finalize(ierror)
                               !    stop
                               ! endif
    endif
                  write(6,'(" exit move2bundle_ iret,n=",2i5)') iret,n

    return

end subroutine move2bundle_

subroutine create_grd23d(grd23d,nlat,nlon,nvert)

   use kinds, only: i_kind
   use general_sub2grid_mod, only: sub2grid_info,general_sub2grid_create_info

   implicit none

   ! Declare local parameters

   ! Declare passed variables
   type(sub2grid_info)               ,intent(inout) :: grd23d
   integer(i_kind)                   ,intent(in   ) :: nlat,nlon,nvert

   ! Declare local variables
   integer(i_kind) inner_vars
   logical regional

   inner_vars=1
   regional=.false.
   call general_sub2grid_create_info(grd23d,inner_vars,nlat,nlon,nvert,nvert,regional)

end subroutine create_grd23d

subroutine destroy_grd23d(grd23d)

   use kinds, only: i_kind
   use general_sub2grid_mod, only: sub2grid_info,general_sub2grid_destroy_info

   implicit none

   ! Declare local parameters

   ! Declare passed variables
   type(sub2grid_info)               ,intent(inout) :: grd23d

   ! Declare local variables

   call general_sub2grid_destroy_info(grd23d)

end subroutine destroy_grd23d

subroutine update_scalar_poles(grd,s)

   use kinds, only: i_kind,r_kind
   use general_sub2grid_mod, only: sub2grid_info,general_sub2grid,general_grid2sub

   implicit none

   ! Declare local parameters

   ! Declare passed variables
   type(sub2grid_info)               ,intent(in   ) :: grd
   real(r_kind)                      ,intent(inout) :: s(grd%lat2,grd%lon2,grd%num_fields)

   ! Declare local variables
   integer(i_kind) inner_vars,lat2,lon2,nlat,nlon,nvert,kbegin_loc,kend_loc,kend_alloc
   integer(i_kind) ii,i,j,k
   real(r_kind),allocatable:: sloc(:),work(:,:,:,:)

   lat2=grd%lat2
   lon2=grd%lon2
   nlat=grd%nlat
   nlon=grd%nlon
   nvert=grd%num_fields
   inner_vars=grd%inner_vars
   kbegin_loc=grd%kbegin_loc
   kend_loc=grd%kend_loc
   kend_alloc=grd%kend_alloc
   allocate(sloc(lat2*lon2*nvert))
   allocate(work(inner_vars,nlat,nlon,kbegin_loc:kend_alloc))
   ii=0
   do k=1,nvert
      do j=1,lon2
         do i=1,lat2
            ii=ii+1
            sloc(ii)=s(i,j,k)
         enddo
      enddo
   enddo
   call general_sub2grid(grd,sloc,work)

   do k=kbegin_loc,kend_loc
      call fillpoles_s(work(1,:,:,k),nlon,nlat)
   enddo
   call general_grid2sub(grd,work,sloc)
   ii=0
   do k=1,nvert
      do j=1,lon2
         do i=1,lat2
            ii=ii+1
            s(i,j,k)=sloc(ii)
         enddo
      enddo
   enddo

   deallocate(sloc,work)

end subroutine update_scalar_poles

subroutine update_vector_poles(grd,u,v,clons,slons)

   use kinds, only: i_kind,r_kind
   use constants, only: zero
   use general_sub2grid_mod, only: sub2grid_info,general_sub2grid,general_grid2sub

   implicit none

   ! Declare local parameters

   ! Declare passed variables
   type(sub2grid_info)               ,intent(in   ) :: grd
   real(r_kind)                      ,intent(inout) :: u(grd%lat2,grd%lon2,grd%num_fields)
   real(r_kind)                      ,intent(inout) :: v(grd%lat2,grd%lon2,grd%num_fields)
   real(r_kind)                      ,intent(in   ) :: clons(grd%nlon),slons(grd%nlon)

   ! Declare local variables
   integer(i_kind) inner_vars,lat2,lon2,nlat,nlon,nvert,kbegin_loc,kend_loc,kend_alloc
   integer(i_kind) ii,i,j,k
   real(r_kind),allocatable:: uloc(:),uwork(:,:,:,:)
   real(r_kind),allocatable:: vloc(:),vwork(:,:,:,:)
   real(r_kind),allocatable:: tempu(:,:),tempv(:,:)

   lat2=grd%lat2
   lon2=grd%lon2
   nlat=grd%nlat
   nlon=grd%nlon
   nvert=grd%num_fields
   inner_vars=grd%inner_vars
   kbegin_loc=grd%kbegin_loc
   kend_loc=grd%kend_loc
   kend_alloc=grd%kend_alloc
   allocate(uloc(lat2*lon2*nvert))
   allocate(vloc(lat2*lon2*nvert))
   allocate(uwork(inner_vars,nlat,nlon,kbegin_loc:kend_alloc))
   allocate(vwork(inner_vars,nlat,nlon,kbegin_loc:kend_alloc))
   allocate(tempu(nlat,nlon),tempv(nlat,nlon))
   uwork=zero ; vwork=zero ; uloc=zero ; vloc=zero
   ii=0
   do k=1,nvert
      do j=1,lon2
         do i=1,lat2
            ii=ii+1
            uloc(ii)=u(i,j,k)
            vloc(ii)=v(i,j,k)
         enddo
      enddo
   enddo
   call general_sub2grid(grd,uloc,uwork)
   call general_sub2grid(grd,vloc,vwork)

   do k=kbegin_loc,kend_loc
      do j=1,nlon
         do i=1,nlat
            tempu(i,j)=uwork(1,i,j,k)
            tempv(i,j)=vwork(1,i,j,k)
         enddo
      enddo
      call fillpoles_v(tempu,tempv,nlon,nlat,clons,slons)
      do j=1,nlon
         do i=1,nlat
            uwork(1,i,j,k)=tempu(i,j)
            vwork(1,i,j,k)=tempv(i,j)
         enddo
      enddo
   enddo
   call general_grid2sub(grd,uwork,uloc)
   call general_grid2sub(grd,vwork,vloc)
   ii=0
   do k=1,nvert
      do j=1,lon2
         do i=1,lat2
            ii=ii+1
            u(i,j,k)=uloc(ii)
            v(i,j,k)=vloc(ii)
         enddo
      enddo
   enddo

   deallocate(uloc,uwork,tempu)
   deallocate(vloc,vwork,tempv)

end subroutine update_vector_poles

subroutine ens_io_partition(n_ens,io_pe,n_io_pe_s,n_io_pe_e,n_io_pe_em,npe,i_ens,mype)

!     do computation on all processors, then assign final local processor
!     values.

      use kinds, only: r_kind,i_kind
      use constants, only: half
      implicit none

!     Declare passed variables
      integer(i_kind),intent(in   ) :: n_ens,npe,mype
      integer(i_kind),intent(  out) :: io_pe,n_io_pe_s,n_io_pe_e,n_io_pe_em,i_ens

!     Declare local variables
      integer(i_kind) :: io_pe0(n_ens)
      integer(i_kind) iskip,jskip,nextra,ipe,n
          integer(i_kind) nsig

      i_ens=-1
      nsig=1
      iskip=npe/n_ens
      nextra=npe-iskip*n_ens
      jskip=iskip
      io_pe=-1
      io_pe0=-1
      n_io_pe_s=1
      n_io_pe_e=0

      ipe=0
      do n=1,n_ens
         io_pe0(n)=ipe
         if(n <= nextra) then
            jskip=iskip+1
         else
            jskip=iskip
         endif
         ipe=ipe+jskip
      enddo
    ! do n=1,n_ens
    !    if(mype==0) write(6,'(" n, io_pe0(n)=",2i8)') n,io_pe0(n)
    ! enddo

      do n=1,n_ens
         if(mype==io_pe0(n)) then
            i_ens=n
            io_pe=mype
            n_io_pe_s=(n-1)*nsig+1
            n_io_pe_e=n*nsig
         endif
      enddo
      n_io_pe_em=max(n_io_pe_s,n_io_pe_e)

end subroutine ens_io_partition

subroutine parallel_read_nemsio_state(en_full,m_cvars2d,m_cvars3d, &
                                        ias ,jas ,mas ,mae  , &
                                        iasm,iaemz,jasm,jaemz,kasm,kaemz,masm,maemz, &
                                        filename,init_head)

   use kinds, only: i_kind,r_kind,r_single
   use constants, only: r60,r3600,zero,one,fv,half,pi,deg2rad
   use mpimod, only: mype
   use nemsio_module, only: nemsio_init,nemsio_open,nemsio_close
   use ncepnems_io, only: error_msg
   use nemsio_module, only: nemsio_gfile,nemsio_getfilehead,nemsio_readrecv
   use nemsio_module, only: nemsio_getrechead
    use general_sub2grid_mod, only: sub2grid_info
    use ensemble_fastread_commvar, only: vars2d,vars3d,n2d,n3d
    use ensemble_fastread_commvar, only: nlon_gfs,nlat_gfs,nsig_gfs,num_fields_gfs

   implicit none

   ! Declare local parameters

   ! Declare passed variables
   integer(i_kind)                   ,intent(in   ) :: ias ,jas ,mas ,mae  
   integer(i_kind)                   ,intent(in   ) :: iasm,iaemz,jasm,jaemz,kasm,kaemz,masm,maemz
   integer(i_kind)                   ,intent(inout) :: m_cvars2d(n2d),m_cvars3d(n3d)
   real(r_single)                    ,intent(inout) :: en_full(iasm:iaemz,jasm:jaemz,kasm:kaemz,masm:maemz)
   character(*)                      ,intent(in   ) :: filename
   logical                           ,intent(in   ) :: init_head

   ! Declare local variables
   integer(i_kind) nlon,nlat,nsig
   integer(i_kind) i,ii,j,jj,k,lonb,latb,levs
   integer(i_kind) k2,k2ps,k2z,k3,k3u,k3v,k3t,k3q,kf
   integer(i_kind) iret,istop
   integer(i_kind),dimension(7):: idate
   integer(i_kind),dimension(4):: odate
   integer(i_kind) nframe,nfhour,nfminute,nfsecondn,nfsecondd
            integer(i_kind) nrec
   character(len=120) :: my_name = 'parallel_read_nemsio_state'
   character(len=1)   :: null = ' '
   real(r_single) :: work(nlon_gfs*(nlat_gfs-2))

   real(r_single) :: work2(nlat_gfs,nlon_gfs)
   real(r_kind) :: fhour
   type(nemsio_gfile) :: gfile

   nlon=nlon_gfs
   nlat=nlat_gfs
   nsig=nsig_gfs

   if ( init_head)call nemsio_init(iret=iret)
   if (iret /= 0) call error_msg(mype,trim(my_name),trim(filename),null,'init',istop,iret)

   call nemsio_open(gfile,filename,'READ',iret=iret)
   if (iret /= 0) call error_msg(mype,trim(my_name),trim(filename),null,'open',istop+1,iret)

   call nemsio_getfilehead(gfile,iret=iret, nframe=nframe, &
        nfhour=nfhour, nfminute=nfminute, nfsecondn=nfsecondn, nfsecondd=nfsecondd, &
        idate=idate, dimx=lonb, dimy=latb,dimz=levs,nrec=nrec)

   if (  nframe /= 0 ) then
      if ( mype == 0 ) &
         write(6,*)trim(my_name),': ***ERROR***  nframe /= 0 for global model read, nframe = ', nframe
      call stop2(101)
   endif

!  check nlat, nlon against latb, lonb

   if(nlat /= latb+2 .or. nlon /= lonb) then
      if ( mype == 0 ) & 
         write(6,*)trim(my_name),': ***ERROR*** incorrect resolution, nlat,nlon=',nlat,nlon, &
                               ', latb+2,lonb=',latb+2,lonb
      call stop2(101)
   endif

!    order of vars on ensemble nems file:

!                       clwmr   -- not used
!                       hgt     -- used      (z)     vars2d(1)
!                       o3mr    -- not used
!                       pres    -- used      (ps)    vars2d(2)
!                       spfh    -- used      (q)     vars3d(4)
!                       tmp     -- used      (tv)    vars3d(3)
!                       ugrd    -- used      (u)     vars3d(1)
!                       vgrd    -- used      (v)     vars3d(2)

   fhour = float(nfhour) + float(nfminute)/r60 + float(nfsecondn)/float(nfsecondd)/r3600
   odate(1) = idate(4)  !hour
   odate(2) = idate(2)  !month
   odate(3) = idate(3)  !day
   odate(4) = idate(1)  !year

!??????????????????????BEGIN NEW CODE????????????????????
!    set up m_cvars3d, m_cvars2d:
   kf=0
   do k3=1,n3d
      m_cvars3d(k3)=kf+1
      do k=1,nsig
         kf=kf+1
      enddo
   enddo
   do k2=1,n2d
      m_cvars2d(k2)=kf+1
      kf=kf+1
   enddo

 ! vars2d(1)='z   '
!   read z (hgt)
   call nemsio_readrecv(gfile,'hgt','sfc',1,work,iret=iret)
   if (iret /= 0) call error_msg(mype,trim(my_name),trim(filename),'hgt','read',istop+9,iret)
   call move1(work,work2,nlon,nlat)
   kf=m_cvars2d(1)       !  z location
   jj=jas-1
   do j=1,nlon
      jj=jj+1
      ii=ias-1
      do i=1,nlat
         ii=ii+1
         en_full(ii,jj,kf,mas)=work2(i,j)
      enddo
   enddo

 ! vars2d(2)='ps  '
!   read ps (pres)
   call nemsio_readrecv(gfile,'pres','sfc',1,work,iret=iret)
   if (iret /= 0) call error_msg(mype,trim(my_name),trim(filename),'pres','read',istop+8,iret)
   call move1(work,work2,nlon,nlat)
   kf=m_cvars2d(2)       !  ps location
   jj=jas-1
   do j=1,nlon
      jj=jj+1
      ii=ias-1
      do i=1,nlat
         ii=ii+1
         en_full(ii,jj,kf,mas)=work2(i,j)
      enddo
   enddo

 ! vars3d(4)='q   '
   kf=m_cvars3d(4)-1       !  q location
   do k=1,nsig
      call nemsio_readrecv(gfile,'spfh','mid layer',k,work,iret=iret)
      if (iret /= 0) call error_msg(mype,trim(my_name),trim(filename),'spfh','read',istop+2,iret)
      call move1(work,work2,nlon,nlat)
      kf=kf+1
      jj=jas-1
      do j=1,nlon
         jj=jj+1
         ii=ias-1
         do i=1,nlat
            ii=ii+1
            en_full(ii,jj,kf,mas)=work2(i,j)
         enddo
      enddo
   enddo

 ! vars3d(3)='tv  '
   kf=m_cvars3d(3)-1       !  tv location
   do k=1,nsig
      call nemsio_readrecv(gfile,'tmp','mid layer',k,work,iret=iret)
      if (iret /= 0) call error_msg(mype,trim(my_name),trim(filename),'tmp','read',istop+3,iret)
      call move1(work,work2,nlon,nlat)
      kf=kf+1
      jj=jas-1
      do j=1,nlon
         jj=jj+1
         ii=ias-1
         do i=1,nlat
            ii=ii+1
            en_full(ii,jj,kf,mas)=work2(i,j)
         enddo
      enddo
   enddo

 ! vars3d(1)='u   '
   kf=m_cvars3d(1)-1       !  u location
   do k=1,nsig
      call nemsio_readrecv(gfile,'ugrd','mid layer',k,work,iret=iret)
      if (iret /= 0) call error_msg(mype,trim(my_name),trim(filename),'ugrd','read',istop+4,iret)
      call move1(work,work2,nlon,nlat)
      kf=kf+1
      jj=jas-1
      do j=1,nlon
         jj=jj+1
         ii=ias-1
         do i=1,nlat
            ii=ii+1
            en_full(ii,jj,kf,mas)=work2(i,j)
         enddo
      enddo
   enddo

 ! vars3d(2)='v   '
   kf=m_cvars3d(2)-1       !  v location
   do k=1,nsig
      call nemsio_readrecv(gfile,'vgrd','mid layer',k,work,iret=iret)
      if (iret /= 0) call error_msg(mype,trim(my_name),trim(filename),'vgrd','read',istop+5,iret)
      call move1(work,work2,nlon,nlat)
      kf=kf+1
      jj=jas-1
      do j=1,nlon
         jj=jj+1
         ii=ias-1
         do i=1,nlat
            ii=ii+1
            en_full(ii,jj,kf,mas)=work2(i,j)
         enddo
      enddo
   enddo

!??????????????????????END NEW CODE????????????????????

     if(mas==1.and.mae==1) write(6,'(" exiting parallel_read_nemsio_state, mype,mas,mae=",3i6)') mype,mas,mae

end subroutine parallel_read_nemsio_state

subroutine fillpoles_s(temp,nlon,nlat)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    fillpoles_s  make pole points average of nearest pole row
!   prgmmr: parrish          org: emc/ncep            date: 2016-10-14
!
! abstract:  make pole points average of nearest pole row.
!
! program history log:
!   2016-10-14  parrish  - initial code
!
!   input argument list:
!     temp     - 2-d input array containing gsi global horizontal field 
!     nlon     - number of gsi/gfs longitudes
!     nlat     - number of gsi latitudes (nlat-2 is gfs--no pole points)
!
!   output argument list:
!     temp     - 2-d output array containing gsi global horizontal field with 
!                    pole values set equal to average of adjacent pole rows
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

   use kinds, only: i_kind,r_kind
   use constants, only: zero,one

   integer(i_kind),intent(in   ) :: nlon,nlat
   real(r_kind), intent(inout) :: temp(nlat,nlon)

   integer(i_kind) nlatm1,i
   real(r_kind) sumn,sums,rnlon

!  Compute mean along southern and northern latitudes
   sumn=zero
   sums=zero
   nlatm1=nlat-1
   do i=1,nlon
      sumn=sumn+temp(nlatm1,i)
      sums=sums+temp(2,i)
   end do
   rnlon=one/float(nlon)
   sumn=sumn*rnlon
   sums=sums*rnlon

!  Load means into local work array
   do i=1,nlon
      temp(1,i)   =sums
      temp(nlat,i)=sumn
   end do

end subroutine fillpoles_s

subroutine fillpoles_v(tempu,tempv,nlon,nlat,clons,slons)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    fillpoles_v  create vector values at pole from nearest pole row
!   prgmmr: parrish          org: emc/ncep            date: 2016-10-14
!
! abstract:  create vector values at pole from nearest pole row.
!
! program history log:
!   2016-10-14  parrish  - initial code
!
!   input argument list:
!     tempu    - 2-d input array containing gsi global horizontal westerly vector component
!     tempv    - 2-d input array containing gsi global horizontal easterly vector component
!     nlon     - number of gsi/gfs longitudes
!     nlat     - number of gsi latitudes (nlat-2 is gfs--no pole points)
!
!   output argument list:
!     temp     - 2-d output array containing gsi global horizontal field with 
!                    pole values set equal to average of adjacent pole rows
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

   use kinds, only: i_kind,r_kind
   use constants, only: zero
   use hybrid_ensemble_parameters, only: sp_ens
         use mpimod, only: mype

   integer(i_kind),intent(in   ) :: nlon,nlat
   real(r_kind), intent(inout) :: tempu(nlat,nlon),tempv(nlat,nlon)
   real(r_kind)                      ,intent(in   ) :: clons(nlon),slons(nlon)

   integer(i_kind) i
   real(r_kind) polnu,polnv,polsu,polsv

!  Compute mean along southern and northern latitudes
   polnu=zero
   polnv=zero
   polsu=zero
   polsv=zero
   do i=1,nlon
      polnu=polnu+tempu(nlat-1,i)*clons(i)-tempv(nlat-1,i)*slons(i)
      polnv=polnv+tempu(nlat-1,i)*slons(i)+tempv(nlat-1,i)*clons(i)
      polsu=polsu+tempu(2,i     )*clons(i)+tempv(2,i     )*slons(i)
      polsv=polsv+tempu(2,i     )*slons(i)-tempv(2,i     )*clons(i)
   end do
   polnu=polnu/float(nlon)
   polnv=polnv/float(nlon)
   polsu=polsu/float(nlon)
   polsv=polsv/float(nlon)
   do i=1,nlon
      tempu(nlat,i)= polnu*clons(i)+polnv*slons(i)
      tempv(nlat,i)=-polnu*slons(i)+polnv*clons(i)
      tempu(1,i   )= polsu*clons(i)+polsv*slons(i)
      tempv(1,i   )= polsu*slons(i)-polsv*clons(i)
   end do

end subroutine fillpoles_v

subroutine move1(work,temp,nlon,nlat)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    move1   move gfs lon lat array to gsi lat lon array
!   prgmmr: parrish          org: emc/ncep            date: 2016-10-14
!
! abstract: move gfs lon lat array to gsi lat lon array.
!
! program history log:
!   2016-10-14  parrish  - initial code
!
!   input argument list:
!     work     - 1-d input array containing gfs horizontal field
!     nlon     - number of gsi/gfs longitudes
!     nlat     - number of gsi latitudes (nlat-2 is gfs--no pole points)
!
!   output argument list:
!     temp     - 2-d output array containing gsi global horizontal field
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

    use kinds, only: i_kind,r_kind,r_single
    use constants, only: zero

    implicit none

    integer(i_kind),intent(in   ) :: nlon,nlat
    real(r_single), intent(in   ) :: work(nlon*(nlat-2))
    real(r_single), intent(  out)   :: temp(nlat,nlon)

    integer(i_kind) ii,i,j

    ii=0
    temp(1,:)   =zero
    temp(nlat,:)=zero
    do i=nlat-1,2,-1
       do j=1,nlon
          ii=ii+1
          temp(i,j)=work(ii)
       enddo
    enddo

end subroutine move1
subroutine ensemble_preread
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    ensemble_preread  initial fast parallel ensemble read
!   prgmmr: parrish          org: np22                date: 2017-01-07
!
! abstract: Do fast parallel read of nems format global ensemble members and
!           redistribute to subdomains in temporary array en_loc3.
!
! program history log:
!   2017-01-07  parrish, initial documentation
!   input argument list:

!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

   use mpimod, only: mpi_comm_world,ierror,mype,npe,mpi_integer4,mpi_real8
   use kinds, only: i_kind,r_kind,r_single
   use constants, only: half,pi,deg2rad
   use constants, only: max_varname_length,i_missing
   use hybrid_ensemble_parameters,only: ensemble_path,n_ens,ntlevs_ens,uv_hyb_ens
   use nemsio_module, only: nemsio_init,nemsio_open,nemsio_close
   use ncepnems_io, only: error_msg
   use nemsio_module, only: nemsio_gfile,nemsio_getfilehead,nemsio_readrecv
   use general_sub2grid_mod, only: sub2grid_info,general_sub2grid_create_info,general_sub2grid_destroy_info
   use guess_grids, only: ges_tsen,ifilesig,hrdifsig,ntguessig
   use obsmod, only: iadate
   use ensemble_fastread_commvar, only: en_loc3,clons,slons,m_cvars2d,m_cvars3d
   use ensemble_fastread_commvar, only: n2d,n3d,vars2d,vars3d,filenames
   use ensemble_fastread_commvar, only: lat1_gfs,lon1_gfs,lat2_gfs,lon2_gfs
   use ensemble_fastread_commvar, only: istart_gfs,jstart_gfs
   use ensemble_fastread_commvar, only: nlon_gfs,nlat_gfs,nsig_gfs,num_fields_gfs

   implicit none

   ! Declare passed variables

   ! Declare internal variables
   integer(i_kind) lonb,latb,levs,j,m,n,ierr
   integer(i_kind) ic2,ic3,k
   integer(i_kind) mas,inner_vars
   integer(i_kind) iret,istop
   integer(i_kind),dimension(7):: idate
   integer(i_kind) nframe,nfhour,nfminute,nfsecondn,nfsecondd
   integer(i_kind) idsl
   integer(i_kind) nrec
   character(255) filename
   logical init_head
   logical regional
   logical,allocatable::vector(:)
   character(len=120) :: my_name = 'ensemble_preread'
   character(len=1)   :: null = ' '
   type(nemsio_gfile) :: gfile_atm
   type(sub2grid_info) grd_ens_temp
   real(r_kind),allocatable :: rlats(:),rlons(:)
   real(r_single),allocatable :: r4lats(:),r4lons(:)

!   set vars2d, vars3d
   n2d=2
   n3d=4
   if(allocated(vars2d)) deallocate(vars2d)
   if(allocated(vars3d)) deallocate(vars3d)
   allocate(vars2d(n2d),vars3d(n3d))
   vars2d(1)='z   '
   vars2d(2)='ps  '
   vars3d(1)='u   '
   vars3d(2)='v   '
   vars3d(3)='tv  '
   vars3d(4)='q   '

!   set initially to zero to avoid debug compile flag.

   lonb=0 ; latb=0
!  HARDWIRE TO ntlevs_ens=1 -- 3denvar only:

          if(ntlevs_ens /= 1) then
             write(6,'(" FAILURE IN ENSEMBLE_PREREAD, HARDWIRED FOR 3DENVAR--PROGRAM STOPS ")')
             if(1/=0) then
                call mpi_finalize(ierror)
                stop
             endif
          endif

       !   iadate doesn't exist yet !!
         !   write(6,*) ' iadate(y,m,d,hr,min)=',iadate

   filename='filelist03'
   if(allocated(filenames)) deallocate(filenames)
   allocate(filenames(n_ens))
   open(10,file=trim(filename),form='formatted',iostat=ierr)
          if(ierr/=0) then
             write(6,'(" FAILURE IN ENSEMBLE_PREREAD, file filelist03 missing--PROGRAM STOPS ")')
             if(1/=0) then
                call mpi_finalize(ierror)
                stop
             endif
          endif
   rewind (10)
   do n=1,n_ens
      read(10,'(a)',iostat=ierr) filenames(n)
          if(ierr/=0) then
             write(6,'(" FAILURE IN ENSEMBLE_PREREAD, problem with filelist03--PROGRAM STOPS ")')
             if(1/=0) then
                call mpi_finalize(ierror)
                stop
             endif
          endif
   enddo

!  obtain nlat, nlon, levs

   if(mype==0) then
      mas=1
      init_head=.true.
      if ( init_head)call nemsio_init(iret=iret)
                write(6,*)' at 3 in ensemble_preread, filenames(1)=',trim(filenames(1))

      call nemsio_open(gfile_atm,filenames(1),'READ',iret=iret)
           idate         = i_missing
           nfhour        = i_missing; nfminute      = i_missing
           nfsecondn     = i_missing; nfsecondd     = i_missing
           idsl  = i_missing
      if (iret /= 0) call error_msg(mype,trim(my_name),trim(filenames(1)),null,'open',istop+1,iret)

      call nemsio_getfilehead(gfile_atm,iret=iret, nframe=nframe, &
           nfhour=nfhour, nfminute=nfminute, nfsecondn=nfsecondn, nfsecondd=nfsecondd, &
           idate=idate, dimx=lonb, dimy=latb,dimz=levs,nrec=nrec)
           if ( nfhour == i_missing .or. nfminute == i_missing .or. &
                nfsecondn == i_missing .or. nfsecondd == i_missing ) then
              write(6,*)'READ_FILES: ***ERROR*** some forecast hour info ', &
                 'are not defined in ', trim(filename)
              write(6,*)'READ_FILES: nfhour = ', nfhour
              call stop2(80)
           endif

      if (  nframe /= 0 ) then
         if ( mype == 0 ) &
            write(6,*)trim(my_name),': ***ERROR***  nframe /= 0 for global model read, nframe = ', nframe
         call stop2(101)
      endif
      nlat_gfs=latb+2 ; nlon_gfs=lonb ; nsig_gfs=levs
   endif
   call mpi_bcast(nlat_gfs,1,mpi_integer4,0,mpi_comm_world,ierror)
   call mpi_bcast(nlon_gfs,1,mpi_integer4,0,mpi_comm_world,ierror)
   call mpi_bcast(nsig_gfs,1,mpi_integer4,0,mpi_comm_world,ierror)
   allocate(rlats(nlat_gfs),rlons(nlon_gfs),r4lats(lonb*latb),r4lons(lonb*latb))

!  obtain r4lats,r4lons,rlats,rlons,clons,slons exactly as computed in
!  general_read_gfsatm_nems:

   if(allocated(clons)) deallocate(clons)
   if(allocated(slons)) deallocate(slons)
   allocate(clons(nlon_gfs),slons(nlon_gfs))
   if(mype==0) then
      call nemsio_getfilehead(gfile_atm,lat=r4lats,iret=iret)
      call nemsio_getfilehead(gfile_atm,lon=r4lons,iret=iret)
      do j=1,latb
         rlats(latb+2-j)=deg2rad*r4lats(lonb/2+(j-1)*lonb)
      enddo
      do j=1,lonb
         rlons(j)=deg2rad*r4lons(j)
      enddo
      rlats(1)=-half*pi
      rlats(latb+2)=half*pi
      do j=1,lonb
         clons(j)=cos(rlons(j))
         slons(j)=sin(rlons(j))
      enddo
   endif
   deallocate(r4lats,r4lons)
   call mpi_bcast(clons,nlon_gfs,mpi_real8,0,mpi_comm_world,ierror)
   call mpi_bcast(slons,nlon_gfs,mpi_real8,0,mpi_comm_world,ierror)
   num_fields_gfs=max(0,n3d)*nsig_gfs+max(0,n2d)

   allocate(vector(num_fields_gfs))
   vector=.false.
   do ic3=1,n3d
      if(trim(vars3d(ic3))=='sf'.or.trim(vars3d(ic3))=='vp' .or. &
         trim(vars3d(ic3))=='u'.or.trim(vars3d(ic3))=='v') then
         do k=1,nsig_gfs
            vector((ic3-1)*nsig_gfs+k)=uv_hyb_ens
         end do
      end if
   end do

   regional=.false.
   inner_vars=1
   call general_sub2grid_create_info(grd_ens_temp,inner_vars,nlat_gfs,nlon_gfs,nsig_gfs, &
                                     num_fields_gfs,regional,vector)

      if(allocated(istart_gfs)) deallocate(istart_gfs)
      if(allocated(jstart_gfs)) deallocate(jstart_gfs)
      allocate(istart_gfs(npe),jstart_gfs(npe))
      lat1_gfs=grd_ens_temp%lat1
      lon1_gfs=grd_ens_temp%lon1
      lat2_gfs=grd_ens_temp%lat2
      lon2_gfs=grd_ens_temp%lon2
      nsig_gfs=grd_ens_temp%nsig
      istart_gfs=grd_ens_temp%istart
      jstart_gfs=grd_ens_temp%jstart
      if(allocated(en_loc3)) deallocate(en_loc3)
      allocate(en_loc3(lat2_gfs,lon2_gfs,num_fields_gfs,n_ens,ntlevs_ens))
      if(allocated(m_cvars2d)) deallocate(m_cvars2d)
      if(allocated(m_cvars3d)) deallocate(m_cvars3d)
      allocate(m_cvars2d(n2d))
      allocate(m_cvars3d(n3d))
      do m=1,ntlevs_ens
         call get_user_ens_fastread_(m,iret)
      enddo
      if(mype==0) write(6,'(" END ENSEMBLE_PREREAD")')

   end subroutine ensemble_preread

   subroutine ensemble_postread
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    ensemble_postread  transfer ensemble to bundle
!   prgmmr: parrish          org: np22                date: 2017-01-07
!
! abstract: Complete fast parallel read by transferring ensemble members from
!           temporary array en_loc3 to atm_bundle(?spelling?).  Deallocate
!           en_loc3.
!
! program history log:
!   2017-01-07  parrish, initial documentation
!   input argument list:

!   output argument list:
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

   use mpimod, only: mype
   use ensemble_fastread_commvar, only: en_loc3,clons,slons,m_cvars2d,m_cvars3d
   use ensemble_fastread_commvar, only: n2d,n3d,vars2d,vars3d,filenames
   use ensemble_fastread_commvar, only: lat1_gfs,lon1_gfs,lat2_gfs,lon2_gfs
   use ensemble_fastread_commvar, only: istart_gfs,jstart_gfs
   use ensemble_fastread_commvar, only: nlon_gfs,nlat_gfs,nsig_gfs,num_fields_gfs

   implicit none

   if(mype==0) write(6,'(" ensemble_postread: deallocate temporary storage in ensemble_fast_commvar")')
   if(allocated(en_loc3)) deallocate(en_loc3)
   if(allocated(clons)) deallocate(clons)
   if(allocated(slons)) deallocate(slons)
   if(allocated(m_cvars2d)) deallocate(m_cvars2d)
   if(allocated(m_cvars3d)) deallocate(m_cvars3d)
   if(allocated(vars2d)) deallocate(vars2d)
   if(allocated(vars3d)) deallocate(vars3d)
   if(allocated(filenames)) deallocate(filenames)
   if(allocated(istart_gfs)) deallocate(istart_gfs)
   if(allocated(jstart_gfs)) deallocate(jstart_gfs)
   if(allocated(m_cvars3d)) deallocate(m_cvars3d)

   end subroutine ensemble_postread
