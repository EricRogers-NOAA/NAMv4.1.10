program nems_merge_nest_4x_step1

  use nemsio_module
  use interp_coef_bgrid

  implicit none

  integer :: iret
  integer :: i,j,l
  type(nemsio_gfile) :: gfile
  character(len=256) :: infile,outfile
  character(len=100) :: inout

  real, allocatable, dimension(:) :: vdata

  integer :: im,jm,lm
  integer :: ie,je,le,ler
  real :: tlm0d,tph0d,dlmd,dphd,wbd,sbd
  real :: pdtop,pt
  real, allocatable, dimension(:) :: dsg1,dsg2
  real, allocatable, dimension(:,:) :: pd,sm,z0,fis
  real, allocatable, dimension(:,:) :: glat,glon
  real, allocatable, dimension(:,:) :: hlat,hlon,vlat,vlon
  real, allocatable, dimension(:,:,:) :: t,q,u,v,u1,v1
  real, allocatable, dimension(:,:,:) :: pint,zint
  real, allocatable, dimension(:,:,:) :: temp, temp1
  real, allocatable, dimension(:,:) :: fi1,fi2
  real, allocatable, dimension(:,:) :: u10,v10,sped10
  real :: ratio

  real :: start_time,stop_time

!  infile = "nmm_b_history.01_nemsio.000h_00m_00.00s"
!  outfile='mmb_d01'
!
! initialize nemsio 
!
  READ(5,*) infile,inout,outfile 
 
  print*,'infile=',infile
  print*,'inout=',inout
  print*,'outfile=',outfile

 call nemsio_init(iret=iret)
  if ( iret /= 0 ) then
    write(0,*) 'nemsio_init, iret=',iret
    stop
  endif

! open and read nemsio file
!

  IF(trim(inout).eq.'relocate') THEN 

  call nemsio_open(gfile,trim(infile),'READ',iret=iret)
  if ( iret /= 0 ) then
    write(0,*) 'nemsio_open, iret=',iret, trim(infile)
    stop
  endif
    write(0,*) 'nemsio_open relocate, iret=',iret, trim(infile)

  call nemsio_getfilehead(gfile,dimx=im,dimy=jm,dimz=lm,iret=iret); call check_err(iret,"dimx,...")

  allocate (dsg1(lm))
  allocate (dsg2(lm))

  call nemsio_getheadvar(gfile,"TLM0D",tlm0d,iret=iret); call check_err(iret,"TLM0D")
  call nemsio_getheadvar(gfile,"TPH0D",tph0d,iret=iret); call check_err(iret,"TPH0D")
  call nemsio_getheadvar(gfile,"DLMD",dlmd,iret=iret); call check_err(iret,"DLMD")
  call nemsio_getheadvar(gfile,"DPHD",dphd,iret=iret); call check_err(iret,"DPHD")
  call nemsio_getheadvar(gfile,"PT",pt,iret=iret); call check_err(iret,"PT")
  call nemsio_getheadvar(gfile,"PDTOP",pdtop,iret=iret); call check_err(iret,"PDTOP")
  call nemsio_getheadvar(gfile,"DSG1",dsg1,iret=iret); call check_err(iret,"DSG1")
  call nemsio_getheadvar(gfile,"DSG2",dsg2,iret=iret); call check_err(iret,"DSG2")

  allocate (vdata(im*jm))

  allocate (hlat(im,jm))
  allocate (hlon(im,jm))
  allocate (vlat(im,jm))
  allocate (vlon(im,jm))
  allocate (pd(im,jm))
  allocate (sm(im,jm))
  allocate (z0(im,jm))
  allocate (fis(im,jm))
  allocate (t(im,jm,lm))
  allocate (q(im,jm,lm))
  allocate (u(im,jm,lm))
  allocate (v(im,jm,lm))
  allocate (pint(im,jm,lm+1))
  allocate (zint(im,jm,lm+1))
  allocate (u10(im,jm))
  allocate (v10(im,jm))
  allocate (sped10(im,jm))

  call nemsio_readrecv(gfile,'glat','sfc',1,vdata,iret=iret); call check_err(iret,'glat')
  hlat(:,:) = reshape(vdata, (/ im,jm /))
      print*,'nems_bin_io_degree.exe hlat(1,1)= ', hlat(1,1)
 print*,'nems_bin_io_degree.exe hlat(1,1)=',hlat(1,1),hlat(1,2),hlat(2,1),hlat(2,2)
  hlat = hlat*180.0/3.14159265359
 print*,'nems_bin_io_degree.exe hlat(1,1)=',hlat(1,1),hlat(1,2),hlat(2,1),hlat(2,2)

  call nemsio_readrecv(gfile,'glon','sfc',1,vdata,iret=iret); call check_err(iret,'glon')
  hlon(:,:) = reshape(vdata, (/ im,jm /))
  hlon = hlon*180.0/3.14159265359

  call nemsio_readrecv(gfile,'vlat','sfc',1,vdata,iret=iret); call check_err(iret,'vlat')
  vlat(:,:) = reshape(vdata, (/ im,jm /))
  vlat = vlat*180.0/3.14159265359

  call nemsio_readrecv(gfile,'vlon','sfc',1,vdata,iret=iret); call check_err(iret,'vlon')
  vlon(:,:) = reshape(vdata, (/ im,jm /))
  vlon = vlon*180.0/3.14159265359

  call nemsio_readrecv(gfile,'dpres','hybrid sig lev',1,vdata,iret=iret); call check_err(iret,'dpres')
  pd(:,:) = reshape(vdata, (/ im,jm /))

  call nemsio_readrecv(gfile,'sm','sfc',1,vdata,iret=iret); call check_err(iret,'sm')
  sm(:,:) = reshape(vdata, (/ im,jm /))

  call nemsio_readrecv(gfile,'zorl','sfc',1,vdata,iret=iret); call check_err(iret,'zorl')
  z0(:,:) = reshape(vdata, (/ im,jm /))

  call nemsio_readrecv(gfile,'hgt','sfc',1,vdata,iret=iret); call check_err(iret,'fis')
  fis(:,:) = reshape(vdata, (/ im,jm /))

  do l=1,lm

    call nemsio_readrecv(gfile,'tmp','mid layer',l,vdata,iret=iret); call check_err(iret,'tmp')
    t(:,:,l) = reshape(vdata, (/ im,jm /))

    call nemsio_readrecv(gfile,'spfh','mid layer',l,vdata,iret=iret); call check_err(iret,'spfh')
    q(:,:,l) = reshape(vdata, (/ im,jm /))

    call nemsio_readrecv(gfile,'ugrd','mid layer',l,vdata,iret=iret); call check_err(iret,'ugrd')
    u(:,:,l) = reshape(vdata, (/ im,jm /))

    call nemsio_readrecv(gfile,'vgrd','mid layer',l,vdata,iret=iret); call check_err(iret,'vgrd')
    v(:,:,l) = reshape(vdata, (/ im,jm /))

  end do

  do l=1,lm+1

    call nemsio_readrecv(gfile,'pres','layer',l,vdata,iret=iret); call check_err(iret,'pint')
    pint(:,:,l) = reshape(vdata, (/ im,jm /))

  end do

  call nemsio_readrecv(gfile,'u10','10 m above gnd',1,vdata,iret=iret); call check_err(iret,'u10')
  u10(:,:) = reshape(vdata, (/ im,jm /))

  call nemsio_readrecv(gfile,'v10','10 m above gnd',1,vdata,iret=iret); call check_err(iret,'v10')
  v10(:,:) = reshape(vdata, (/ im,jm /))

  sped10(:,:)=sqrt(u10(:,:)*u10(:,:)+v10(:,:)*v10(:,:))

  call nemsio_close(gfile,iret=iret)
  if ( iret /= 0 ) then
    write(0,*) 'nemsio_close, iret=',iret
    stop
  endif

  call nemsio_finalize()

  allocate (fi1(im,jm))
  allocate (fi2(im,jm))

  fi1(:,:)=fis(:,:)
  zint(:,:,lm+1)=fis(:,:)/9.8
  do l=lm,1,-1
    fi2(:,:)=t(:,:,l)*(q(:,:,l)*0.608+1.0)*287.04*               &
             (alog(pint(:,:,l+1))-ALOG(pint(:,:,l)))+fi1(:,:)
    zint(:,:,l) = fi2(:,:)/9.8
    fi1(:,:) = fi2(:,:)
  end do

!! WRITE DATA FOR 3DVAR
      open(77,file=trim(outfile),form='unformatted')
      WRITE(77)im,jm,lm
      WRITE(77)dlmd,dphd,TLM0D,TPH0D
      WRITE(77)pt,pdtop
      WRITE(77)(((T(i,j,l),i=1,im),j=1,jm),l=lm,1,-1)
      WRITE(77)(((Q(i,j,l),i=1,im),j=1,jm),l=lm,1,-1)
      WRITE(77)(((U(i,j,l),i=1,im),j=1,jm),l=lm,1,-1)
      WRITE(77)(((V(i,j,l),i=1,im),j=1,jm),l=lm,1,-1)
      WRITE(77)(((ZINT(i,j,l),i=1,im),j=1,jm),l=lm+1,1,-1)
      WRITE(77)HLON,HLAT,VLON,VLAT
      WRITE(77)(((PINT(i,j,l),i=1,im),j=1,jm),l=lm+1,1,-1)
      WRITE(77)PD
      WRITE(77)DSG1
      WRITE(77)DSG2
      WRITE(77)SM
      WRITE(77)Z0
      WRITE(77)sped10
      PRINT *,'before close WRITE(77) outfile= ', trim(outfile)
      CLOSE(77)

      do l=1,lm
        print*,'l,dsg1,dsg2=',l,dsg1(l),dsg2(l)
      end do
 
  else if (trim(inout).eq.'update') THEN
  
  call cpu_time(start_time)

  call nemsio_open(gfile,trim(infile),'RDWR',iret=iret)
  if ( iret /= 0 ) then
    write(0,*) 'nemsio_open, iret=',iret, trim(infile)
    stop
  endif
    write(0,*) 'nemsio_open update , iret=',iret, trim(infile)

  call cpu_time(stop_time)
  print *," after nemsio_open ",stop_time-start_time

  call nemsio_getfilehead(gfile,dimx=im,dimy=jm,dimz=lm,iret=iret); call check_err(iret,"dimx,...")

!  allocate (dsg1(lm))
!  allocate (dsg2(lm))

!  call nemsio_getheadvar(gfile,"TLM0D",tlm0d,iret=iret); call check_err(iret,"TLM0D")
!  call nemsio_getheadvar(gfile,"TPH0D",tph0d,iret=iret); call check_err(iret,"TPH0D")
!  call nemsio_getheadvar(gfile,"DLMD",dlmd,iret=iret); call check_err(iret,"DLMD")
!  call nemsio_getheadvar(gfile,"DPHD",dphd,iret=iret); call check_err(iret,"DPHD")
!  call nemsio_getheadvar(gfile,"PT",pt,iret=iret); call check_err(iret,"PT")
!  call nemsio_getheadvar(gfile,"PDTOP",pdtop,iret=iret); call check_err(iret,"PDTOP")
!  call nemsio_getheadvar(gfile,"DSG1",dsg1,iret=iret); call check_err(iret,"DSG1")
!  call nemsio_getheadvar(gfile,"DSG2",dsg2,iret=iret); call check_err(iret,"DSG2")

  allocate (vdata(im*jm))

  allocate (hlat(im,jm))
  allocate (hlon(im,jm))
  allocate (vlat(im,jm))
  allocate (vlon(im,jm))
  allocate (pd(im,jm))
  allocate (sm(im,jm))
  allocate (z0(im,jm))
  allocate (fis(im,jm))
  allocate (t(im,jm,lm))
  allocate (q(im,jm,lm))
  allocate (u(im,jm,lm),v(im,jm,lm))
  allocate (u1(im,jm,lm),v1(im,jm,lm))
  allocate (pint(im,jm,lm+1))
  allocate (zint(im,jm,lm+1))
  allocate (u10(im,jm))
  allocate (v10(im,jm))

  allocate (temp(im,jm,lm),temp1(im,jm,lm+1))

!! WRITE DATA FOR 3DVAR
      call cpu_time(start_time)
      open(77,file=trim(outfile),form='unformatted')
      call cpu_time(stop_time)
      print *," after file 77 open ",stop_time-start_time
      call cpu_time(start_time)
      READ(77)    ! im,jm,lm
      READ(77)    ! dlmd,dphd,TLM0D,TPH0D
      READ(77)    ! pt,pdtop
!     READ(77)(((T(i,j,l),i=1,im),j=1,jm),l=lm,1,-1)
!     READ(77)(((Q(i,j,l),i=1,im),j=1,jm),l=lm,1,-1)
!     READ(77)(((U(i,j,l),i=1,im),j=1,jm),l=lm,1,-1)
!     READ(77)(((V(i,j,l),i=1,im),j=1,jm),l=lm,1,-1)
!     temperature
      READ(77) TEMP
      do le=1,lm 
        ler=lm-le+1
        do je=1,jm
          do ie=1,im
            t(ie,je,ler)=temp(ie,je,le)
          enddo
        enddo
      enddo
!     Q
      READ(77) TEMP
      do le=1,lm
        ler=lm-le+1
        do je=1,jm
          do ie=1,im
            q(ie,je,ler)=temp(ie,je,le)
          enddo
        enddo
      enddo
!     u
      READ(77) TEMP
      do le=1,lm
        ler=lm-le+1
        do je=1,jm
          do ie=1,im
            u(ie,je,ler)=temp(ie,je,le)
          enddo
        enddo
      enddo
!     V
      READ(77) TEMP
      do le=1,lm
        ler=lm-le+1
        do je=1,jm
          do ie=1,im
            v(ie,je,ler)=temp(ie,je,le)
          enddo
        enddo
      enddo
      READ(77)    ! (((ZINT(i,j,l),i=1,im),j=1,jm),l=lm+1,1,-1)
      READ(77)    ! HLON,HLAT,VLON,VLAT
!     READ(77)(((PINT(i,j,l),i=1,im),j=1,jm),l=lm+1,1,-1)
!     PINT
      READ(77) TEMP1
      do le=1,lm+1
        ler=lm+1-le+1
        do je=1,jm
          do ie=1,im
            pint(ie,je,ler)=temp1(ie,je,le)
          enddo
        enddo
      enddo
      READ(77)PD  
      READ(77)    ! DSG1
      READ(77)    ! DSG2
!      READ(77)SM
!      READ(77)Z0
      PRINT *,'before close READ(77) readin from ', trim(outfile)
      CLOSE(77)
   
      deallocate (temp)
      deallocate (temp1)

  call cpu_time(stop_time)
  print *," after read_77 ",stop_time-start_time

  call cpu_time(start_time)

!  call nemsio_readrecv(gfile,'glat','sfc',1,vdata,iret=iret); call check_err(iret,'glat')
!  hlat(:,:) = reshape(vdata, (/ im,jm /))

!  call nemsio_readrecv(gfile,'glon','sfc',1,vdata,iret=iret); call check_err(iret,'glon')
!  hlon(:,:) = reshape(vdata, (/ im,jm /))

!  call nemsio_readrecv(gfile,'vlat','sfc',1,vdata,iret=iret); call check_err(iret,'vlat')
!  vlat(:,:) = reshape(vdata, (/ im,jm /))

!  call nemsio_readrecv(gfile,'vlon','sfc',1,vdata,iret=iret); call check_err(iret,'vlon')
!  vlon(:,:) = reshape(vdata, (/ im,jm /))

  vdata(:) = reshape(pd, (/ im*jm /))
  call nemsio_writerecv(gfile,'dpres','hybrid sig lev',1,vdata,iret=iret); call check_err(iret,'dpres')

!  call nemsio_readrecv(gfile,'sm','sfc',1,vdata,iret=iret); call check_err(iret,'sm')
!  sm(:,:) = reshape(vdata, (/ im,jm /))

!  call nemsio_readrecv(gfile,'zorl','sfc',1,vdata,iret=iret); call check_err(iret,'zorl')
!  z0(:,:) = reshape(vdata, (/ im,jm /))

!  call nemsio_readrecv(gfile,'hgt','sfc',1,vdata,iret=iret); call check_err(iret,'fis')
!  fis(:,:) = reshape(vdata, (/ im,jm /))

  call nemsio_readrecv(gfile,'u10','10 m above gnd',1,vdata,iret=iret); call check_err(iret,'u10')
  u10(:,:) = reshape(vdata, (/ im,jm /))

  call nemsio_readrecv(gfile,'v10','10 m above gnd',1,vdata,iret=iret); call check_err(iret,'v10')
  v10(:,:) = reshape(vdata, (/ im,jm /))

  call cpu_time(stop_time)
  print *," after read2d ",stop_time-start_time

  call cpu_time(start_time)

  do l=1,lm

    call nemsio_readrecv(gfile,'ugrd','mid layer',l,vdata,iret=iret); call check_err(iret,'ugrd')
    u1(:,:,l) = reshape(vdata, (/ im,jm /))

    call nemsio_readrecv(gfile,'vgrd','mid layer',l,vdata,iret=iret); call check_err(iret,'vgrd')
    v1(:,:,l) = reshape(vdata, (/ im,jm /))

  end do

!u10 and v10 do not do relocation because they involve fairly complicated
!bondary layer model calculations. for temporary fix, we put a .95 to lowest
!level ugrd and vgrd. They are pretty close to each other.
! Once NAM does a digital filter or one step integration, it will be removed.
!Guang Ping Lou, Nov-2013
!  do j=1,jm
!  do i=1,im
!    ratio =sqrt((u(i,j,lm)**2+v(i,j,lm)**2)/(u1(i,j,lm)**2+v1(i,j,lm)**2+1.e-10))
!    ratio = 0.95
!    u10(i,j)=u(i,j,lm)*ratio
!    v10(i,j)=v(i,j,lm)*ratio
!  end do
!  end do

  do l=1,lm

    vdata(:) = reshape(t(:,:,l), (/ im*jm /))
    call nemsio_writerecv(gfile,'tmp','mid layer',l,vdata,iret=iret); call check_err(iret,'tmp')

    vdata(:) =  reshape(q(:,:,l), (/ im*jm /))
    call nemsio_writerecv(gfile,'spfh','mid layer',l,vdata,iret=iret); call check_err(iret,'spfh')

    vdata(:) =  reshape(u(:,:,l), (/ im*jm /))
    call nemsio_writerecv(gfile,'ugrd','mid layer',l,vdata,iret=iret); call check_err(iret,'ugrd')

    vdata(:) =  reshape(v(:,:,l), (/ im*jm /))
    call nemsio_writerecv(gfile,'vgrd','mid layer',l,vdata,iret=iret); call check_err(iret,'vgrd')

  end do

  do l=1,lm+1

    vdata(:) =  reshape(pint(:,:,l), (/ im*jm /))
    call nemsio_writerecv(gfile,'pres','layer',l,vdata,iret=iret); call check_err(iret,'pint')

  end do

  call cpu_time(stop_time)
  print *," after read3d ",stop_time-start_time

  vdata(:) = reshape(u10, (/ im*jm /))
  call nemsio_writerecv(gfile,'u10','10 m above gnd',1,vdata,iret=iret); call check_err(iret,'u10')

  vdata(:) = reshape(v10, (/ im*jm /))
  call nemsio_writerecv(gfile,'v10','10 m above gnd',1,vdata,iret=iret); call check_err(iret,'v10')

  call nemsio_close(gfile,iret=iret)
  if ( iret /= 0 ) then
    write(0,*) 'nemsio_close, iret=',iret
    stop
  endif

  call cpu_time(start_time)

  call nemsio_finalize()

  call cpu_time(stop_time)
  print *," after finalize ",stop_time-start_time

  end if

contains

   subroutine check_err(ierr,clog)
   integer, intent(in) :: ierr
   character(len=*), optional, intent(in) :: clog

   if (ierr /= 0) then
     write(0,*)' ierr ',ierr
     if (present(clog)) write(0,*)' clog = ', clog
     stop 1
!   else
!     if (present(clog)) write(0,*)' clog = ', clog, ' OK'
   end if

   end subroutine check_err
end program nems_merge_nest_4x_step1
