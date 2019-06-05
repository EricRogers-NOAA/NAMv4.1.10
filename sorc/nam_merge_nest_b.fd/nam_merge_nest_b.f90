
      PROGRAM HWRF_MERGE_NEST

!??????????????????????????????????????????????????????????
!     
! ABSTRACT: CREATE 3x DEGREE2 HIGH RESOLUTION DOMAIN, the center of the domain
!           is the same as the newly created wrfinput_d01 and wrfanl_d02
!     
! ORIGINAL AUTHOR: QINGFU LIU, NCEP/EMC, 2007
! REVISED  AUTHOR: XUEJIN ZHANG (AOML/HRD), Qingfu Liu, JULY 2011 
!                : Extend the initialization to 3km
! REVISED  AUTHOR: Qingfu Liu, Dusan Jovic, Jan. 2012
!                 : Extend the initialization to 3km, NMMB grid
!
!     DECLARE VARIABLES
!
!      IMPLICIT NONE

      INTEGER I,J,K,NX,NY,NZ,NST,IFLAG_NEST,JX,JY,KZ
!
      PARAMETER (NST=5,MMB=1)
!      PARAMETER (IX=117,IY=225)
!      PARAMETER (NX=215,NY=431,NZ=42,NST=5)
!      PARAMETER (JX=393,JY=735)       ! fixed for 9 km resolution
      PARAMETER (GAMMA=6.5E-3,G=9.8,Rd=287.05,D608=0.608)
      PARAMETER (Cp=1004.)

! Variables on new outer nest hybrid coordinate (d01)

      REAL(4)  DLMD0,DPHD0,WBD0,SBD0,CLON0,CLAT0
      REAL(4), ALLOCATABLE :: HLON0(:,:),HLAT0(:,:),VLON0(:,:),VLAT0(:,:)

! Variables on new inner nest hybrid coordinate (d03)

      REAL(4), ALLOCATABLE :: HLON4(:,:),HLAT4(:,:),VLON4(:,:),VLAT4(:,:)

! Variables for 3x working domain (30x30 degress)

      REAL(4)  CLON_NHC,CLAT_NHC
      REAL(4)  DLMD3,DPHD3,PT3,PDTOP3              ! use the new inner nest data
      REAL(4)  WBD3,SBD3
      REAL(8)  rlon0,rlat0

      REAL(4), ALLOCATABLE :: HLON3(:,:),HLAT3(:,:)
      REAL(4), ALLOCATABLE :: VLON3(:,:),VLAT3(:,:)

      REAL(4), ALLOCATABLE :: T3(:,:,:),Q3(:,:,:)
      REAL(4), ALLOCATABLE :: U3(:,:,:),V3(:,:,:)
      REAL(4), ALLOCATABLE :: Z3(:,:,:),P3(:,:,:)
      REAL(4), ALLOCATABLE :: PD3(:,:),PMID3(:,:,:)
      REAL(4), ALLOCATABLE :: PMV3(:,:,:)
      REAL(4), ALLOCATABLE :: A13(:,:),B13(:,:),C13(:,:)

      REAL(4), ALLOCATABLE :: SLP3(:,:)
      REAL(4), ALLOCATABLE :: ZS3(:,:),TS3(:,:),QS3(:,:)

! working arrays used for outer nest interpolation  (parent)

      integer(4), ALLOCATABLE :: IIH(:,:,:),JJH(:,:,:)
      integer(4), ALLOCATABLE :: IIV(:,:,:),JJV(:,:,:)
      REAL(4),    ALLOCATABLE :: HBWGT(:,:,:),VBWGT(:,:,:)

! working arrays used for inner nest interpolation  (d02)

      integer(4), ALLOCATABLE :: IIH1(:,:,:),JJH1(:,:,:)
      integer(4), ALLOCATABLE :: IIV1(:,:,:),JJV1(:,:,:)
      REAL(4),    ALLOCATABLE :: HBWGT1(:,:,:),VBWGT1(:,:,:)

! working arrays used for inner nest interpolation  (d03)

      integer(4), ALLOCATABLE :: IIH2(:,:,:),JJH2(:,:,:)
      integer(4), ALLOCATABLE :: IIV2(:,:,:),JJV2(:,:,:)
      REAL(4),    ALLOCATABLE :: HBWGT2(:,:,:),VBWGT2(:,:,:)

! working array

      REAL(4), ALLOCATABLE :: T21(:,:,:,:),Q21(:,:,:,:)
      REAL(4), ALLOCATABLE :: U21(:,:,:,:),V21(:,:,:,:)
      REAL(4), ALLOCATABLE :: SLP21(:,:,:)

      integer(4) IH1(4),JH1(4),IV1(4),JV1(4)

! Variables from outer nest (d01)

      REAL(4) DLMD1,DPHD1,PT1,PDTOP1
      REAL(4) WBD1,SBD1,CLON1,CLAT1

      REAL(4), ALLOCATABLE :: T1(:,:,:),Q1(:,:,:)
      REAL(4), ALLOCATABLE :: U1(:,:,:),V1(:,:,:) 
      REAL(4), ALLOCATABLE :: Z1(:,:,:),P1(:,:,:)
      REAL(4), ALLOCATABLE :: HLON1(:,:),HLAT1(:,:)
      REAL(4), ALLOCATABLE :: VLON1(:,:),VLAT1(:,:)
      REAL(4), ALLOCATABLE :: DSG1(:),DSG2(:)
      REAL(4), ALLOCATABLE :: PD1(:,:)
      REAL(4), ALLOCATABLE :: A101(:,:),B101(:,:)

      REAL(4), ALLOCATABLE :: SLP1(:,:)
      REAL(4), ALLOCATABLE :: PMID1(:,:,:),ZMID1(:,:,:)
      REAL(4), ALLOCATABLE :: PMV1(:,:,:)
    
! Variables from inner nest (d02)

      REAL(4) DLMD2,DPHD2,PT2,PDTOP2
      REAL(4) WBD2,SBD2,CLON2,CLAT2

      REAL(4), ALLOCATABLE :: T2(:,:,:),Q2(:,:,:)
      REAL(4), ALLOCATABLE :: U2(:,:,:),V2(:,:,:)
      REAL(4), ALLOCATABLE :: Z2(:,:,:),P2(:,:,:)
      REAL(4), ALLOCATABLE :: HLON2(:,:),HLAT2(:,:)
      REAL(4), ALLOCATABLE :: VLON2(:,:),VLAT2(:,:)
      REAL(4), ALLOCATABLE :: PD2(:,:)
      REAL(4), ALLOCATABLE :: A102(:,:),B102(:,:),C102(:,:)
   
      REAL(4), ALLOCATABLE :: SLP2(:,:)
      REAL(4), ALLOCATABLE :: PMID2(:,:,:),ZMID2(:,:,:)
      REAL(4), ALLOCATABLE :: PMV2(:,:,:)

! Variables from inner nest (d03)

      REAL(4) DLMD5,DPHD5,PT5,PDTOP5
      REAL(4) WBD5,SBD5,CLON5,CLAT5

      REAL(4), ALLOCATABLE :: T5(:,:,:),Q5(:,:,:)
      REAL(4), ALLOCATABLE :: U5(:,:,:),V5(:,:,:)
      REAL(4), ALLOCATABLE :: Z5(:,:,:),P5(:,:,:)
      REAL(4), ALLOCATABLE :: HLON5(:,:),HLAT5(:,:)
      REAL(4), ALLOCATABLE :: VLON5(:,:),VLAT5(:,:)
      REAL(4), ALLOCATABLE :: PD5(:,:)
      REAL(4), ALLOCATABLE :: A105(:,:),B105(:,:),C105(:,:)
   
      REAL(4), ALLOCATABLE :: SLP5(:,:)
      REAL(4), ALLOCATABLE :: PMID5(:,:,:),ZMID5(:,:,:)
      REAL(4), ALLOCATABLE :: PMV5(:,:,:)
      
      CHARACTER*1 SN,EW

!!!!!!!!!!!!!!!!11


      COEF1=Rd/Cp
      COEF3=Rd*GAMMA/G
      COEF2=1./COEF3

      GRD=G/Rd

      pi=4.*atan(1.)
      pi_deg=180./pi
      pi180=1./pi_deg

      DIST1=6.371E3*pi180

      READ(5,*)ITIM,IVOBS,IBGS,IFLAG_NEST,CLAT0,CLON0         ! CLAT0,CLON0 is the new domain center

      print*,'ITIM,IVOBS,IBGS,IFLAG_NEST,CLAT0,CLON0=',    &
              ITIM,IVOBS,IBGS,IFLAG_NEST,CLAT0,CLON0

      if(CLON0.gt.60.)CLON0=CLON0-360.

      IF(MMB.EQ.1)THEN

! READ PARENT DATA (d01)   ! 6 hour forecast data

      IUNIT=20+ITIM

      READ(IUNIT) KX0,KY0,NZ0

      print*,'KX0,KY0,NZ0=',KX0,KY0,NZ0

      ALLOCATE ( HLON1(KX0,KY0),HLAT1(KX0,KY0) )

      READ(IUNIT) DLMD0,DPHD0,CLON1,CLAT1
      READ(IUNIT)
      READ(IUNIT)
      READ(IUNIT)
      READ(IUNIT)
      READ(IUNIT)
      READ(IUNIT)
      READ(IUNIT) HLON1,HLAT1
      READ(IUNIT)
      READ(IUNIT)

      rewind(IUNIT)

      IF(CLON1.GT.60.)CLON1=CLON1-360.
      DO J=1,KY0
      DO I=1,KX0
        IF(HLON1(I,J).GT.60.)HLON1(I,J)=HLON1(I,J)-360.
      END DO
      END DO

! check domain center

      dc_dist=abs(CLON1-CLON0)+abs(CLAT1-CLAT0)
      if(dc_dist.gt.0.01)then
        print*,'WARNING! domain center is different, check input data'
        CLON0=CLON1
        CLAT0=CLAT1
      end if

      print*,'DLMD0,DPHD0,CLON0,CLAT0=',DLMD0,DPHD0,CLON0,CLAT0

      WBD0=-((KX0-1.)/2.)*DLMD0                  ! PARENT wbd
      SBD0=-((KY0-1.)/2.)*DPHD0                  ! PARENT SBD

      print*,'WBD0,SBD0=',WBD0,SBD0

      ALLOCATE ( HLON0(KX0,KY0),HLAT0(KX0,KY0),VLON0(KX0,KY0),VLAT0(KX0,KY0) )

      CALL EARTH_LATLON_BGRID( HLAT0,HLON0,VLAT0,VLON0,     &  !Earth lat,lon at H and V points
                               DLMD0,DPHD0,WBD0,SBD0,       &  !input res,west & south boundaries,
                               CLAT0,CLON0,                 &  ! central lat,lon, all in degrees
                               KX0,KY0 )

      print*,'this is check:'
      DO J=1,KY0
      DO I=1,KX0
         DIF1=abs(HLON1(I,J)-HLON0(I,J))+abs(HLAT1(I,J)-HLAT0(I,J))
         IF(DIF1.GT.0.001)THEN
            print*,'need to check hlat1 and hlon1'
            print*,'DIF1=',DIF1,I,J,HLON1(I,J),HLON0(I,J),HLAT1(I,J),HLAT0(I,J)
            stop
         END IF
      END DO
      END DO

      ELSE
     
      print*,'This version is for B-grid, it does not work with E-grid'
      stop

      END IF

! READ INNER NEST DATA (d02)   ! 6 hour forecast data

      IUNIT=30+ITIM

      READ(IUNIT) KX,KY,KZ              ! KZ==NZ

      print*,'KX,KY,KZ=',KX,KY,KZ

      ALLOCATE ( HLON2(KX,KY),HLAT2(KX,KY) )

      READ(IUNIT) DLMD2,DPHD2,CLON2,CLAT2
      READ(IUNIT) ! PT2,PDTOP2
      READ(IUNIT) ! T2
      READ(IUNIT) ! Q2
      READ(IUNIT) ! U2
      READ(IUNIT) ! V2
      READ(IUNIT) ! Z2
      READ(IUNIT) HLON2,HLAT2
      READ(IUNIT) ! P2
      READ(IUNIT) ! PD2
      READ(IUNIT) ! DSG1
      READ(IUNIT) ! DSG2
!        READ(IUNIT) A102
!        READ(IUNIT) B102
!        READ(IUNIT) C102
      REWIND(IUNIT)

      IF(CLON2.GT.60.)CLON2=CLON2-360.
      DO J=1,KY
      DO I=1,KX
        IF(HLON2(I,J).GT.60.)HLON2(I,J)=HLON2(I,J)-360.
      END DO
      END DO

      ERR=1.e20
      DO J=1,KY0
      DO I=1,KX0
        DIF1=abs(HLON2(1,1)-HLON1(I,J))+abs(HLAT2(1,1)-HLAT1(I,J))
        IF(DIF1.LT.ERR)THEN
          ILOC=I
          JLOC=J
          ERR=DIF1
        END IF
      END DO
      END DO

!      WBD2= WBD0 + (ILOC -1)*DLMD0
!      SBD2= SBD0 + (JLOC -1)*DPHD0
      WBD2=-((KX-1.)/2.)*DLMD2                  ! PARENT wbd
      SBD2=-((KY-1.)/2.)*DPHD2                  ! PARENT SBD

      print *, 'ILOC,JLOC,ERR',ILOC,JLOC,ERR
      print *, 'WBD2,SBD2',WBD2,SBD2

      IF(IFLAG_NEST.eq.2)THEN

! READ INNER NEST DATA (d03)   ! 6 hour forecast data

      IUNIT=40+ITIM

      READ(IUNIT) KX3,KY3,KZ3              ! KZ==NZ

      print*,'KX3,KY3,KZ3=',KX3,KY3,KZ3

      ALLOCATE ( HLON5(KX3,KY3),HLAT5(KX3,KY3) )

      READ(IUNIT) DLMD3,DPHD3,CLON3,CLAT3
      READ(IUNIT) 
      READ(IUNIT) 
      READ(IUNIT) 
      READ(IUNIT) 
      READ(IUNIT)
      READ(IUNIT) 
      READ(IUNIT) HLON5,HLAT5
      READ(IUNIT) 
      READ(IUNIT) 
      READ(IUNIT) 
      READ(IUNIT) 
      REWIND(IUNIT)

      IF(CLON3.GT.60.)CLON3=CLON3-360.
      DO J=1,KY3
      DO I=1,KX3
        IF(HLON5(I,J).GT.60.)HLON5(I,J)=HLON5(I,J)-360.
      END DO
      END DO

      ERR=1.e20
      DO J=1,KY
      DO I=1,KX
        DIF1=abs(HLON5(1,1)-HLON2(I,J))+abs(HLAT5(1,1)-HLAT2(I,J))
        IF(DIF1.LT.ERR)THEN
          ILOC=I
          JLOC=J
          ERR=DIF1
        END IF
      END DO
      END DO

      dc_dist=abs(CLON3-CLON0)+abs(CLAT3-CLAT0)
      if(dc_dist.gt.0.01)then
        print*,'WARNING! domain center is different, check input data'
        CLON3=CLON0
        CLAT3=CLAT0
      end if

!      WBD5= WBD2 + (ILOC -1)*DLMD2
!      SBD5= SBD2 + (JLOC -1)*DPHD2
      WBD2=-((KX-1.)/2.)*DLMD2                  ! PARENT wbd
      SBD2=-((KY-1.)/2.)*DPHD2                  ! PARENT SBD

      print*,'WBD5,SBD5=',WBD5,SBD5


! LON & LAT at T,U,V
        
      ALLOCATE ( HLON4(KX3,KY3),HLAT4(KX3,KY3),VLON4(KX3,KY3),VLAT4(KX3,KY3) )

       CALL EARTH_LATLON_BGRID( HLAT4,HLON4,VLAT4,VLON4,     &  !Earth lat,lon at H and V points
                                DLMD3,DPHD3,WBD5,SBD5,       &  !input res,west & south boundaries,
                                CLAT0,CLON0,                 &  ! central lat,lon, all in degrees
                                KX3,KY3 )

       print*,'HLAT4,HLON4,VLAT4,VLON4,1,1=',                  &
               HLAT4(1,1),HLON4(1,1),VLAT4(1,1),VLON4(1,1)
       print*,'HLAT4,HLON4,VLAT4,VLON4,KX3,KY3=',                  &
               HLAT4(KX3,KY3),HLON4(KX3,KY3),VLAT4(KX3,KY3),VLON4(KX3,KY3)
       print*,'HLAT5,HLON5,1,1=',                  &
               HLAT5(1,1),HLON5(1,1)
       print*,'HLAT5,HLON5,KX3,KY3=',                  &
               HLAT5(KX3,KY3),HLON5(KX3,KY3)
        
       DO J=1,KY3
       DO I=1,KX3
          DIF1=abs(HLON4(I,J)-HLON5(I,J))+abs(HLAT4(I,J)-HLAT5(I,J))
          IF(DIF1.GT.0.001)THEN
             print*,'need to check hlat41 and hlon41'
             print*,'DIF1=',DIF1
             stop
          END IF
       END DO
       END DO

       END IF     ! IFLAG_NEST

! read in storm center

       read(11,11)ICLAT,SN,ICLON,EW
  11   format(33x,I3,A1,I5,A1)
       CLAT_NHC=ICLAT*0.1
       CLON_NHC=ICLON*0.1

       IF(SN.eq.'S')CLAT_NHC=-CLAT_NHC
       IF(EW.eq.'W')CLON_NHC=-CLON_NHC
!       IF(EW.eq.'E')CLON_NHC=CLON_NHC-360.
       IF(CLON_NHC.GT.60.)CLON_NHC=CLON_NHC-360.
!       if(CLON0.gt.60.)CLON_NHC=CLON_NHC-360.

      print *, 'CLON_NHC,CLAT_NHC=',CLON_NHC,CLAT_NHC

! create 3x high resolution domain

!      DLMD3=DLMD3
!      DPHD3=DPHD3

      DLMD3=DLMD2
      DPHD3=DPHD2

      JX=INT(32./DLMD3)-1                   ! hard wired for relocation domain and for the bogus storm
      JY=INT(32./DLMD3)-1

!      IF(IFLAG_NEST.EQ.0)THEN
!        DLMD3=DLMD0
!        DPHD3=DPHD0
!        JX=INT(31./DLMD3)+1
!        JY=INT(31./DLMD3)+1
!      END IF

      JX1=JX+1
      JY1=JY+1
      NZ =NZ0
      NZ1=NZ+1

      print*,'DLMD0,DPHD0=',DLMD0,DPHD0,DLMD3,DPHD3
      print*,'JX,JY=',JX,JY

      ERR=1.e20 
      DO J=1,KY
      DO I=1,KX
        DIF1=abs(CLON_NHC-HLON2(I,J))+abs(CLAT_NHC-HLAT2(I,J))
        IF(DIF1.LT.ERR)THEN
          ILOC=I
          JLOC=J
          ERR=DIF1
        END IF
      END DO
      END DO

      print*,'large ERR is OK here'
      print*,'inside merge_nest1_3n,ILOC,JLOC,ERR=',ILOC,JLOC,ERR
!      print*,'inside merge_nest1_3n,HLON1(ILOC,JLOC),HLAT1(ILOC,JLOC)='
!      print*,HLON1(ILOC,JLOC),HLAT1(ILOC,JLOC)
!      print*,'inside merge_nest1_3n,HLON2(ILOC,JLOC),HLAT2(ILOC,JLOC)='
!      print*,HLON2(ILOC,JLOC),HLAT2(ILOC,JLOC)

! make JLOC odd number

!      WBD3= WBD0 + (ILOC -1)*DLMD0
!      SBD3= SBD0 + (JLOC -1)*DPHD0


!!      iloc=-99
!!      jloc=-99
       slon0=HLON2(ILOC,JLOC)
       slat0=HLAT2(ILOC,JLOC)
!!      slon0=CLON_NHC
!!     slat0=CLAT_NHC
! check to see if the hurricane is on outside of the domain edge. 
! if yes, stop here. (use temp variables WBD3 and SBD3) - Guang Ping Lou
!
      SBD3 = slat0-CLAT_NHC
      WBD3 = slon0-CLON_NHC
      print*, 'slon0,slat0, CLAT_NHC,CLON_NHC, dlat, dlon = '
      print*, slon0,slat0,CLAT_NHC,CLON_NHC, SBD3,WBD3
!      if ( abs(SBD3) .gt. 0.4 .or. abs(WBD3) .gt. 0.4 ) then
      if ( ILOC .gt. (KX-8) .or. ILOC .lt. 8 .or. JLOC .lt. 8 ) then
         print*, 'the storm is either on edge or outside '
         write(12,*) 'OUT' 
        stop
        else
         print*, 'the storm is in '
         write(12,*) 'IN' 
      endif
 
      clon05=clon2
      clat05=clat2
      print*, 'slon0,slat0,HLON2(ILOC,JLOC),HLAT2(ILOC,JLOC)= '
      print*, slon0,slat0,HLON2(ILOC,JLOC),HLAT2(ILOC,JLOC)
      call TLL(slon0,slat0,rlon0,rlat0,clat05,clon05)
      print*, 'slon0,slat0,rlon0,rlat0,clat05,clon05'
      print*,  slon0,slat0,rlon0,rlat0,clat05,clon05

      WBD3 = rlon0 - int(JX/2)*dlmd3
      SBD3 = rlat0 - int(JY/2)*dphd3

      write(*,*)'DLMD0,DPHD0=',DLMD0,DPHD0
      write(*,*)'WBD3,SBD3,CLON0,CLAT0=',    &
                 WBD3,SBD3,CLON0,CLAT0
      print*,'ILOC,JLOC,WBD3,SBD3=',ILOC,JLOC,WBD3,SBD3

      ALLOCATE ( HLON3(JX,JY),HLAT3(JX,JY) )
      ALLOCATE ( VLON3(JX,JY),VLAT3(JX,JY) )

      ALLOCATE ( T3(JX,JY,NZ),Q3(JX,JY,NZ) )
      ALLOCATE ( U3(JX,JY,NZ),V3(JX,JY,NZ) )
      ALLOCATE ( Z3(JX,JY,NZ1),P3(JX,JY,NZ1) )

      ALLOCATE ( PD3(JX,JY), PMID3(JX,JY,NZ) )

      ALLOCATE ( SLP3(JX,JY) )
      ALLOCATE ( ZS3(JX,JY),TS3(JX,JY),QS3(JX,JY) )

      ALLOCATE ( IIH(JX,JY,4),JJH(JX,JY,4) )
      ALLOCATE ( IIV(JX,JY,4),JJV(JX,JY,4) )
      ALLOCATE ( HBWGT(JX,JY,4),VBWGT(JX,JY,4) )

      ALLOCATE ( IIH1(JX,JY,4),JJH1(JX,JY,4) )
      ALLOCATE ( IIV1(JX,JY,4),JJV1(JX,JY,4) )
      ALLOCATE ( HBWGT1(JX,JY,4),VBWGT1(JX,JY,4) )

      ALLOCATE ( IIH2(JX,JY,4),JJH2(JX,JY,4) )
      ALLOCATE ( IIV2(JX,JY,4),JJV2(JX,JY,4) )
      ALLOCATE ( HBWGT2(JX,JY,4),VBWGT2(JX,JY,4) )

      ALLOCATE ( A13(JX,JY),B13(JX,JY),C13(JX,JY) )

! working array

      ALLOCATE ( T21(JX,JY,NZ,4),Q21(JX,JY,NZ,4) )
      ALLOCATE ( U21(JX,JY,NZ,4),V21(JX,JY,NZ,4) )
      ALLOCATE ( SLP21(JX,JY,4) )
      ALLOCATE ( PMV3(JX,JY,NZ) )

! LON & LAT at T,U,V

       CALL EARTH_LATLON_BGRID( HLAT3,HLON3,VLAT3,VLON3,     &  !Earth lat,lon at H and V points
                                DLMD3,DPHD3,WBD3,SBD3,       &  !input res,west & south boundaries,
                                CLAT2,CLON2,                 &  ! central lat,lon, all in degrees
                                JX,JY )
    
       print*,'HLAT3,HLON3,VLAT3,VLON3,1,1=',                  &
               HLAT3(1,1),HLON3(1,1),VLAT3(1,1),VLON3(1,1)
       print*,'HLAT3,HLON3,VLAT3,VLON3,JX,JY=',                  &
               HLAT3(JX,JY),HLON3(JX,JY),VLAT3(JX,JY),VLON3(JX,JY)

! 3x vs. parent

      ERR=1.e20
      DO J=1,KY
      DO I=1,KX
        DIF1=abs(HLON3(1,1)-HLON2(I,J))+abs(HLAT3(1,1)-HLAT2(I,J))
!        DIF1=abs(HLON3(1,JY)-HLON2(I,J))+abs(HLAT3(1,JY)-HLAT2(I,J))
        IF(DIF1.LT.ERR)THEN
          I_ST4=I
          J_ST4=J
          ERR=DIF1
        END IF
      END DO
      END DO

      print*,'ERR should be ver small'
      PRINT*,'3x vs. parent, I_ST4,J_ST4,ERR 2=',I_ST4,J_ST4,ERR
      PRINT*,'HLON3(1,JY),HLON2(I_ST4,J_ST4),HLAT3(1,JY),HLAT2(I_ST4,J_ST4)= '
      PRINT*,HLON3(1,JY),HLON2(I_ST4,J_ST4),HLAT3(1,JY),HLAT2(I_ST4,J_ST4)

      IF(IFLAG_NEST.EQ.2)THEN  
  
! find the grid index at the lower left corner
! double check

! inner nest vs. parent (40x40)

      ERR=1.e20
      DO J=1,KY0
      DO I=1,KX0
        DIF1=abs(HLON4(1,1)-HLON0(I,J))+abs(HLAT4(1,1)-HLAT0(I,J))
        IF(DIF1.LT.ERR)THEN
          I_ST=I
          J_ST=J
          ERR=DIF1
        END IF
      END DO
      END DO
                 
      print*,'Large ERR should be OK'
      PRINT*,'nest vs. parent, I_ST,J_ST,ERR 1=',I_ST,J_ST,ERR 

! inner nest vs. 3x

      FLON=HLON4(1,1)
      FLAT=HLAT4(1,1)

      print*,'FLON,FLAT=',FLON,FLAT

      ERR=1.e20
      DO J=1,JY
      DO I=1,JX
        DIF1=abs(FLON-HLON3(I,J))+abs(FLAT-HLAT3(I,J))
        IF(DIF1.LT.ERR)THEN
          I_ST14=I
          J_ST14=J
          ERR=DIF1
        END IF
      END DO
      END DO

      print*,'ERR should be ver small'
      PRINT*,'3x vs. nest, I_ST14,J_ST14,ERR=',I_ST14,J_ST14,ERR


      DO J=1,KY3
      DO I=1,KX3
        I1=I+I_ST14-1
        J1=J+J_ST14-1
        DIFFI=abs(HLON4(I,J)-HLON3(I1,J1))+    &
              abs(HLAT4(I,J)-HLAT3(I1,J1))
       IF(DIFFI.GT.0.001)THEN
         print*,'I,J,I1,J1,DIFFI=',I,J,I1,J1,DIFFI
         print*,'HLON4(I,J),HLON3(I1,J1)',HLON4(I,J),HLON3(I1,J1)
         print*,'HLAT4(I,J),HLAT3(I1,J1)',HLAT4(I,J),HLAT3(I1,J1)
         print*,'try change JLOC to JLOC-1'
         STOP 333
       END IF
      END DO
      END DO

      END IF    ! IFLAG_NEST

!     DEALLOCATE ( HLON0,HLAT0,VLON0,VLAT0 )
!     DEALLOCATE ( HLON0,HLAT0,HLON1,HLAT1,HLON2,HLAT2,HLON5,HLAT5 )
     DEALLOCATE ( HLON0,HLAT0,HLON1,HLAT1,HLON2,HLAT2)

! READ PARENT DATA (d01)   ! 6 hour forecast data

      IUNIT=20+ITIM

      READ(IUNIT) NX,NY,NZ

      print*,'NX,NY,NZ=',NX,NY,NZ

      NZ1=NZ+1

      ALLOCATE ( T1(NX,NY,NZ),Q1(NX,NY,NZ) )
      ALLOCATE ( U1(NX,NY,NZ),V1(NX,NY,NZ) )
      ALLOCATE ( Z1(NX,NY,NZ1),P1(NX,NY,NZ1) )
      ALLOCATE ( HLON1(NX,NY),HLAT1(NX,NY),VLON1(NX,NY),VLAT1(NX,NY) )
      ALLOCATE ( PD1(NX,NY),A101(NX,NY),B101(NX,NY) )
      ALLOCATE ( DSG1(NZ),DSG2(NZ) )

      READ(IUNIT) DLMD1,DPHD1,CLON1,CLAT1
      READ(IUNIT) PT1,PDTOP1
      READ(IUNIT) T1
      READ(IUNIT) Q1
      READ(IUNIT) U1
      READ(IUNIT) V1
      READ(IUNIT) Z1
      READ(IUNIT) HLON1,HLAT1
      READ(IUNIT) P1
      READ(IUNIT) PD1
      READ(IUNIT) DSG1
      READ(IUNIT) DSG2

!      IF(IBGS.NE.2)THEN

        print*,'reading A101'
        READ(IUNIT) A101                        ! A101 = land sea mask, B101 = ZNT
        READ(IUNIT) B101
        print*,'finishing reading A101'
      print*,'d01 PT1,PDTOP1=', PT1,PDTOP1
        DO I=1,NZ
         print*,'NZ, DSG1,DSG2=', I, DSG1(I),DSG2(I)
        ENDDO


!      END IF
 
      CLOSE(IUNIT)

      if(CLON1.gt.60.)CLON1=CLON1-360.
      DO J=1,NY
      DO I=1,NX
        IF(HLON1(I,J).GT.60.)HLON1(I,J)=HLON1(I,J)-360.
      END DO
      END DO

      WBD1=-((NX-1.)/2.)*DLMD1
      SBD1=-((NY-1.)/2.)*DPHD1

!       IF(IFLAG_NEST.NE.1.and.KX0.EQ.NX)THEN

         VLON1=VLON0
         VLAT1=VLAT0
   
         IUNIT=21+ITIM

         WRITE(IUNIT) NX,NY,NZ
         WRITE(IUNIT) DLMD1,DPHD1,CLON1,CLAT1
         WRITE(IUNIT) PT1,PDTOP1
         WRITE(IUNIT) T1
         WRITE(IUNIT) Q1
         WRITE(IUNIT) U1
         WRITE(IUNIT) V1
         WRITE(IUNIT) Z1
         WRITE(IUNIT) HLON1,HLAT1,VLON1,VLAT1
         WRITE(IUNIT) P1
         WRITE(IUNIT) PD1
         WRITE(IUNIT) DSG1
         WRITE(IUNIT) DSG2
         WRITE(IUNIT) A101            ! A101 = land sea mask, B101 = ZNT
         WRITE(IUNIT) B101
!       END IF

       print*,'read in d01, HLON1,HLAT1,VLON1,VLAT1=',     &
              HLON1(1,1),HLAT1(1,1),VLON1(1,1),VLAT1(1,1)

! READ INNER NEST DATA (d02)   ! 6 hour forecast data

      IUNIT=30+ITIM

      READ(IUNIT) IX,IY,IZ              ! IZ==NZ

      print*,'IX,IY,IZ=',IX,IY,IZ

      IZ1=IZ+1

      ALLOCATE ( T2(IX,IY,IZ),Q2(IX,IY,IZ) )
      ALLOCATE ( U2(IX,IY,IZ),V2(IX,IY,IZ) )
      ALLOCATE ( Z2(IX,IY,IZ1),P2(IX,IY,IZ1) )
      print*,'check111'
      ALLOCATE ( HLON2(IX,IY),HLAT2(IX,IY) )
      print*,'check112'
      ALLOCATE ( VLON2(IX,IY),VLAT2(IX,IY) )
      print*,'check113'
      ALLOCATE ( A102(IX,IY),B102(IX,IY),C102(IX,IY) )
      ALLOCATE ( PD2(IX,IY) )

!       A102=1.                 ! all ocean
!       B102=0.002
!       C102=0.88
      DEALLOCATE ( DSG1,DSG2 )
      ALLOCATE ( DSG1(IZ),DSG2(IZ) )
 
      READ(IUNIT) DLMD2,DPHD2,CLON2,CLAT2
      READ(IUNIT) PT2,PDTOP2
      READ(IUNIT) T2
      READ(IUNIT) Q2
      READ(IUNIT) U2
      READ(IUNIT) V2
      READ(IUNIT) Z2
      READ(IUNIT) HLON2,HLAT2,VLON2,VLAT2
      READ(IUNIT) P2
      READ(IUNIT) PD2
      READ(IUNIT) DSG1
      READ(IUNIT) DSG2

!      IF(IBGS.NE.2)THEN
        READ(IUNIT) A102
        READ(IUNIT) B102
        READ(IUNIT) C102
        DO J=1,IY
        DO I=1,IX
          wind_s=sqrt(U2(I,J,1)**2+V2(I,J,1)**2)+1.E-10
          C102(I,J)=min(1.0,C102(I,J)/wind_s)
        END DO
        END DO
!      END IF 

      CLOSE(IUNIT)

      print*,'check1'
      print*,'d02 PT2,PDTOP2=', PT2,PDTOP2
        DO I=1,IZ
         print*,'IZ, DSG1,DSG2=', I, DSG1(I),DSG2(I)
        ENDDO

      IF(CLON2.GT.60.)CLON2=CLON2-360.
      DO J=1,IY
      DO I=1,IX
        IF(HLON2(I,J).GT.60.)HLON2(I,J)=HLON2(I,J)-360.
        IF(VLON2(I,J).GT.60.)VLON2(I,J)=VLON2(I,J)-360.
      ENDDO
      ENDDO

      print*,'read in d02, HLON2,HLAT2,VLON2,VLAT2=',     &
              HLON2(1,1),HLAT2(1,1),VLON2(1,1),VLAT2(1,1)


      ERR=1.e20
!      DO J=1,NY
!      DO I=1,NX
!        DIF1=abs(HLON2(1,1)-HLON1(I,J))+abs(HLAT2(1,1)-HLAT1(I,J))
!        IF(DIF1.LT.ERR)THEN
!          ILOC=I
!          JLOC=J
!          ERR=DIF1
!        END IF
!      END DO
!      END DO

      print*,'new test11111'
      do k=1,NZ1
        print*,k,p1(nx/2,ny/2,k),p2(ix/2,iy/2,k)
      end do


      WBD1=-((NX-1.)/2.)*DLMD1
      SBD1=-((NY-1.)/2.)*DPHD1

!      WBD2= WBD1 + (ILOC -1)*DLMD1
!      SBD2= SBD1 + (JLOC -1)*DPHD1

!      print*,'ERR should be ver small'
      PRINT*,'ILOC,JLOC,ERR=',ILOC,JLOC,ERR

       print*,'IX,IY,NZ=',IX,IY,NZ
       write(*,*)'inside merge K,T2,Q2,U2,V2,Z2,P2='
       do k=1,NZ
         write(*,32)K,T2(9,9,K),            &
           Q2(9,9,K),U2(9,9,K),V2(9,9,K),Z2(9,9,K),P2(9,9,K)
       end do

       IF(IFLAG_NEST.EQ.2)THEN

! READ 6h FORECAST data (d03)

      IUNIT=40+ITIM

      READ(IUNIT)NX5,NY5,NZ5
      print*,'NX5,NY5,NZ5=',NX5,NY5,NZ5

      ALLOCATE ( T5(NX5,NY5,NZ5),Q5(NX5,NY5,NZ5) )
      ALLOCATE ( U5(NX5,NY5,NZ5),V5(NX5,NY5,NZ5) )
      ALLOCATE ( Z5(NX5,NY5,NZ1),P5(NX5,NY5,NZ1) )
      ALLOCATE ( HLON5(NX5,NY5),HLAT5(NX5,NY5) )
      ALLOCATE ( VLON5(NX5,NY5),VLAT5(NX5,NY5) )
      ALLOCATE ( A105(NX5,NY5),B105(NX5,NY5),C105(NX5,NY5) )
      ALLOCATE ( PD5(NX5,NY5) )

!       A105=1.                 ! all ocean
!       B105=0.002
!       C105=0.88

      READ(IUNIT) DLMD5,DPHD5,CLON5,CLAT5
      READ(IUNIT) DT5,PDTOP5
      READ(IUNIT) T5
      READ(IUNIT) Q5
      READ(IUNIT) U5
      READ(IUNIT) V5
      READ(IUNIT) Z5
      READ(IUNIT) HLON5,HLAT5,VLON5,VLAT5
      READ(IUNIT) P5
      READ(IUNIT) PD5
      READ(IUNIT) DSG1
      READ(IUNIT) DSG2

!      IF(IBGS.NE.2)THEN
        READ(IUNIT) A105
        READ(IUNIT) B105
        READ(IUNIT) C105
        DO J=1,NY5
        DO I=1,NX5
          wind_s=sqrt(U5(I,J,1)**2+V5(I,J,1)**2)+1.E-10
          C105(I,J)=min(1.0,C105(I,J)/wind_s)
        END DO
        END DO
!      END IF 

      CLOSE(IUNIT)

      IF(CLON5.GT.60.)CLON5=CLON5-360.
      DO J=1,NY5
      DO I=1,NX5
        IF(HLON5(I,J).GT.60.)HLON5(I,J)=HLON5(I,J)-360.
        IF(VLON5(I,J).GT.60.)VLON5(I,J)=VLON5(I,J)-360.
      ENDDO
      ENDDO

       print*,'read in d03, HLON5,HLAT5,VLON5,VLAT5=',     &
               HLON5(1,1),HLAT5(1,1),VLON5(1,1),VLAT5(1,1)
      ERR=1.e20
      DO J=1,NY5
      DO I=1,NX5
        DIF1=abs(HLON5(1,1)-HLON2(I,J))+abs(HLAT5(1,1)-HLAT2(I,J))
        IF(DIF1.LT.ERR)THEN
          ILOC=I
          JLOC=J
          ERR=DIF1
        END IF
      END DO
      END DO

      WBD5= WBD2 + (ILOC -1)*DLMD2
      SBD5= SBD2 + (JLOC -1)*DPHD2

      IF(IBGS.EQ.1)THEN

         print*,'check if the nest domain created correctly'

         DO J=1,NY5
         DO I=1,NX5
           DIFFI=abs(HLON4(I,J)-HLON5(I,J))+    &
                 abs(HLAT4(I,J)-HLAT5(I,J))
           IF(DIFFI.GT.0.001)THEN
             print*,'I,J,DIFFI=',I,J,DIFFI
             STOP 333
           END IF
         END DO
         END DO

         print*,'finished checking nest domain created correctly'

      END IF

       print*,'NX5,NY5,NZ=',NX5,NY5,NZ
       write(*,*)'inside merge K,T5,Q5,U5,V5,Z5,P5='
       do k=1,NZ
         write(*,32)K,T5(9,9,K),            &
           Q5(9,9,K),U5(9,9,K),V5(9,9,K),Z5(9,9,K),P5(9,9,K)
       end do

      END IF     ! IFLAG_NEST

!!!!!!!!!!!!!!!!!!!!

       print*,'NX,NY,NZ=',NX,NY,NZ
       write(*,*)'inside merge K,T1,Q1,U1,V1,Z1,P1='
       do k=1,NZ
         write(*,32)K,T1(9,9,K),            &
           Q1(9,9,K),U1(9,9,K),V1(9,9,K),Z1(9,9,K),P1(9,9,K)
       end do


      DO J=1,NY
      DO I=1,NX
        IF(HLON1(I,J).gt.60.)HLON1(I,J)=HLON1(I,J)-360.
        IF(VLON1(I,J).gt.60.)VLON1(I,J)=VLON1(I,J)-360.
      END DO
      END DO

      PT3=PT1
      PDTOP3=PDTOP1

      print*,'DLMD1,DPHD1,PT1,PDTOP1=',DLMD1,DPHD1,PT1,PDTOP1
      print*,'HLAT1,HLON1=',HLAT1(1,1),HLON1(1,1)


      ALLOCATE ( SLP1(NX,NY) )
      ALLOCATE ( PMID1(NX,NY,NZ),ZMID1(NX,NY,NZ) )
      ALLOCATE ( PMV1(NX,NY,NZ) )

       DO K=1,NZ
       DO J=1,NY
       DO I=1,NX
         PMID1(I,J,K)=EXP((ALOG(P1(I,J,K))+ALOG(P1(I,J,K+1)))*0.5)
         ZMID1(I,J,K)=0.5*(Z1(I,J,K)+Z1(I,J,K+1))
       ENDDO
       ENDDO
       ENDDO

!C        COMPUTE SEA LEVEL PRESSURE.
!C
       DO J=1,NY
       DO I=1,NX
         ZSF1 = ZMID1(I,J,1)
         PSF1 = PMID1(I,J,1)
         TV1 = T1(I,J,1)*(1.+D608*Q1(I,J,1))
         A = (GAMMA * ZSF1) / TV1
         SLP1(I,J) = PSF1*(1+A)**COEF2
      ENDDO
      ENDDO

      print *,'fort.63','NX=',NX,'NY=',NY

!      WRITE(63)((SLP1(I,J),I=1,NX),J=1,NY)
!      DO K=1,NZ+1
!        WRITE(63)((Z1(I,J,K),I=1,NX),J=1,NY)
!      END DO
!      DO K=1,NZ+1
!       WRITE(63)((P1(I,J,K),I=1,NX),J=1,NY)
!      END DO
!      DO K=1,NZ
!        WRITE(63)((T1(I,J,K),I=1,NX),J=1,NY)
!      END DO
!      DO K=1,NZ
!        WRITE(63)((Q1(I,J,K),I=1,NX),J=1,NY)
!      END DO
!      DO K=1,NZ
!        WRITE(63)((U1(I,J,K),I=1,NX),J=1,NY)
!      END DO
!      DO K=1,NZ
!        WRITE(63)((V1(I,J,K),I=1,NX),J=1,NY)
!      END DO

!  parameters for the outer nest

       WBD1=-((NX-1.)/2.)*DLMD1
       SBD1=-((NY-1.)/2.)*DPHD1

       print*,'CLON1,CLAT1=',CLON1,CLAT1
       print*,'WBD1,SBD1,NX,NY=',WBD1,SBD1,NX,NY

        CALL G2T2H_BGRID( IIH,JJH,             & ! output grid index and weights
                    HBWGT,                         & ! output weights in terms of parent grid
                    HLAT3,HLON3,                   & ! target (nest) input lat lon in degrees
                    DLMD1,DPHD1,WBD1,SBD1,         & ! parent res, western and south boundaries
                    CLAT1,CLON1,                   & ! parent central lat,lon, all in degrees
                    NX,NY,                         & ! parent imax and jmax
                    JX,JY)                           ! target (nest) grid dimensions

        CALL G2T2V_BGRID( IIV,JJV,             & ! output grid index and weights
                    VBWGT,                         & ! output weights in terms of parent grid
                    VLAT3,VLON3,                   & ! target (nest) input lat lon in degrees
                    DLMD1,DPHD1,WBD1,SBD1,         & ! parent res, western and south boundaries
                    CLAT1,CLON1,                   & ! parent central lat,lon, all in degrees
                    NX,NY,                         & ! parent imax and jmax
                    JX,JY)                           ! target (nest) grid dimensions

       DO J=1,JY
       DO I=1,JX
         IF(IIH(I,J,4).GT.NX.or.IIH(I,J,1).LT.1)print*,I,J,IIH(I,J,1)
         IF(JJH(I,J,4).GT.NY.or.JJH(I,J,1).LT.1)print*,I,J,JJH(I,J,1)
       END DO
       END DO

       print*,'check coef'

       DO J=1,JY
       DO I=1,JX
         CSUM=0.
         DO N1=1,4
           CSUM=CSUM+abs(HBWGT(I,J,N1))
         END DO
         IF(CSUM.GT.1.01.or.CSUM.lt.0.99)THEN
           print*,'test failed, I,J,CSUM=',I,J,CSUM
         END IF
       END DO
       END DO

       ZS3=0.
       DO J=1,JY
       DO I=1,JX
         DO N1=1,4
           IH1(N1)=IIH(I,J,N1)
           JH1(N1)=JJH(I,J,N1)
         END DO
         ZS3(I,J) = HBWGT(I,J,1)*Z1(IH1(1),JH1(1),1)             &
                  + HBWGT(I,J,2)*Z1(IH1(2),JH1(2),1)             &
                  + HBWGT(I,J,3)*Z1(IH1(3),JH1(3),1)             &
                  + HBWGT(I,J,4)*Z1(IH1(4),JH1(4),1)
       ENDDO
       ENDDO

       print*,'test1='

       DO J=1,JY
       DO I=1,JX
         Z3(I,J,1)=ZS3(I,J)
       END DO
       END DO

       print*,'test2='

       DO K=1,IZ
         v_maxk=0.
         DO J=1,IY
         DO I=1,IX
  	   v_max1=U2(I,J,K)*U2(I,J,K)+V2(I,J,K)*V2(I,J,K)
	   if(v_maxk.lt.v_max1)v_maxk=v_max1
	 END DO
	 END DO
	 print*,'native grid data K,v_maxk=',K,sqrt(v_maxk)
       END DO

! variables from d02

      DO J=1,IY
      DO I=1,IX
  	IF(HLON2(I,J).gt.60.)HLON2(I,J)=HLON2(I,J)-360.
  	IF(VLON2(I,J).gt.60.)VLON2(I,J)=VLON2(I,J)-360.
      END DO
      END DO

      ALLOCATE ( SLP2(IX,IY) )
      ALLOCATE ( PMID2(IX,IY,IZ),ZMID2(IX,IY,IZ) )
      ALLOCATE ( PMV2(IX,IY,IZ) )

       DO K=1,IZ
       DO J=1,IY
       DO I=1,IX
         PMID2(I,J,K)=EXP((ALOG(P2(I,J,K))+ALOG(P2(I,J,K+1)))*0.5)
         ZMID2(I,J,K)=0.5*(Z2(I,J,K)+Z2(I,J,K+1))
       ENDDO
       ENDDO
       ENDDO

       DO J=1,IY
       DO I=1,IX
         ZSF2 = ZMID2(I,J,1)
         PSF2 = PMID2(I,J,1)
         TV2 = T2(I,J,1)*(1.+D608*Q2(I,J,1))
!C
!C        COMPUTE SEA LEVEL PRESSURE.
         A = (GAMMA * ZSF2) / TV2
         SLP2(I,J) = PSF2*(1.+A)**COEF2
      ENDDO
      ENDDO

       print *,'fort.62','IX=',IX,'IY=',IY,'IZ=',IZ
       print *,'fort.62','DLMD1,DPHD1= ', DLMD1,DPHD1

!      WRITE(62)((SLP2(I,J),I=1,IX),J=1,IY)
!      DO K=1,IZ+1
!        WRITE(62)((Z2(I,J,K),I=1,IX),J=1,IY)
!      END DO
!      DO K=1,IZ+1
!        WRITE(62)((P2(I,J,K),I=1,IX),J=1,IY)
!      END DO
!      DO K=1,IZ
!        WRITE(62)((T2(I,J,K),I=1,IX),J=1,IY)
!      END DO
!      DO K=1,IZ
!        WRITE(62)((Q2(I,J,K),I=1,IX),J=1,IY)
!      END DO
!      DO K=1,IZ
!        WRITE(62)((U2(I,J,K),I=1,IX),J=1,IY)
!      END DO
!      DO K=1,IZ
!        WRITE(62)((V2(I,J,K),I=1,IX),J=1,IY)
!      END DO

       IF(IFLAG_NEST.EQ.2)THEN

! variables from d03

      print*,'test1 qqqqq'

      DO J=1,NY5
      DO I=1,NX5
  	IF(HLON5(I,J) .GT. 60.)HLON5(I,J)=HLON5(I,J)-360.
  	IF(VLON5(I,J) .GT. 60.)VLON5(I,J)=VLON5(I,J)-360.
      END DO
      END DO

      print*,'test2 qqqqq'

      ALLOCATE ( SLP5(NX5,NY5) )
      ALLOCATE ( PMID5(NX5,NY5,NZ5),ZMID5(NX5,NY5,NZ5) )
      ALLOCATE ( PMV5(NX5,NY5,NZ5) )

      print*,'test21 qqqqq'
       DO K=1,NZ5
       DO J=1,NY5
       DO I=1,NX5
         PMID5(I,J,K)=EXP((ALOG(P5(I,J,K))+ALOG(P5(I,J,K+1)))*0.5)
         ZMID5(I,J,K)=0.5*(Z5(I,J,K)+Z5(I,J,K+1))
       ENDDO
       ENDDO
       ENDDO

      print*,'test3 qqqqq'
       DO J=1,NY5
       DO I=1,NX5
         ZSF5 = ZMID5(I,J,1)
         PSF5 = PMID5(I,J,1)
         TV5 = T5(I,J,1)*(1.+D608*Q5(I,J,1))
!C
!C        COMPUTE SEA LEVEL PRESSURE.
         A = (GAMMA * ZSF5) / TV5
         SLP5(I,J) = PSF5*(1.+A)**COEF2
      ENDDO
      ENDDO

      print*,'test4 qqqqq'
      NX22=NX5/2
      NY22=NY5/2
      print*,'SLP5,T5,Z5,P5=',SLP5(NX22,NY22),T5(NX22,NY22,1),     &
              Q5(NX22,NY22,1),Z5(NX22,NY22,1),P5(NX22,NY22,1)

      WRITE(62)((SLP2(I,J),I=1,IX),J=1,IY,2)
      DO K=1,IZ+1
        WRITE(62)((Z2(I,J,K),I=1,IX),J=1,IY,2)
      END DO
      DO K=1,IZ+1
        WRITE(62)((P2(I,J,K),I=1,IX),J=1,IY,2)
      END DO
      DO K=1,IZ
        WRITE(62)((T2(I,J,K),I=1,IX),J=1,IY,2)
      END DO
      DO K=1,IZ
        WRITE(62)((Q2(I,J,K),I=1,IX),J=1,IY,2)
      END DO
      DO K=1,IZ
        WRITE(62)((U2(I,J,K),I=1,IX),J=1,IY,2)
      END DO
      DO K=1,IZ
        WRITE(62)((V2(I,J,K),I=1,IX),J=1,IY,2)
      END DO

      END IF     ! IFLAG_NEST

! compute interpolation for outer  nest

       NX_1=NX-1
       NY_1=NY-1              

! First, compute variables at surface level

       DO J=1,JY
       DO I=1,JX
         IF(IIH(I,J,1).LE.2.or.IIH(I,J,1).GE.NX_1)GO TO 40
         IF(JJH(I,J,1).LE.2.or.JJH(I,J,1).GE.NY_1)GO TO 40
         DO N1=1,4
           IH1(N1)=IIH(I,J,N1)
           JH1(N1)=JJH(I,J,N1)
         END DO
         DO N1=1,4
           SLP21(I,J,N1)=SLP1(IH1(N1),JH1(N1))
         END DO
         K=1             ! surface
         DO N1=1,4
           IF(Z3(I,J,K).LT.ZMID1(IH1(N1),JH1(N1),1))THEN
             DZ1=ZMID1(IH1(N1),JH1(N1),1)-Z3(I,J,K)
             T21(I,J,K,N1)=T1(IH1(N1),JH1(N1),1)+GAMMA*DZ1
             Q21(I,J,K,N1)=Q1(IH1(N1),JH1(N1),1) 
           ELSE IF(Z3(I,J,K).GE.ZMID1(IH1(N1),JH1(N1),NZ))THEN  ! never occur for K=1
             DZ1=ZMID1(IH1(N1),JH1(N1),NZ)-Z3(I,J,K)
             T21(I,J,K,N1)=T1(IH1(N1),JH1(N1),NZ)+GAMMA*DZ1
             Q21(I,J,K,N1)=Q1(IH1(N1),JH1(N1),NZ)
           ELSE
             DO K1=2,NZ
               ZDIF1=ZMID1(IH1(N1),JH1(N1),K1)-Z3(I,J,K)
               IF(ZDIF1.GE.0.)THEN
                 FACT1=ZDIF1/(ZMID1(IH1(N1),JH1(N1),K1)-   &
                              ZMID1(IH1(N1),JH1(N1),K1-1))
                 T21(I,J,K,N1)=T1(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +T1(IH1(N1),JH1(N1),K1-1)*FACT1
                 Q21(I,J,K,N1)=Q1(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +Q1(IH1(N1),JH1(N1),K1-1)*FACT1
                 GO TO 38 
               END IF
             END DO
           END IF
 38        CONTINUE
         END DO
 40      CONTINUE
       ENDDO
       ENDDO

       print*,'test3='

       DO J=1,JY
       DO I=1,JX
         IF(IIH(I,J,1).LE.2.or.IIH(I,J,1).GE.NX_1)GO TO 44
         IF(JJH(I,J,1).LE.2.or.JJH(I,J,1).GE.NY_1)GO TO 44
         SLP3(I,J) = HBWGT(I,J,1)*SLP21(I,J,1)             &
                   + HBWGT(I,J,2)*SLP21(I,J,2)             &
                   + HBWGT(I,J,3)*SLP21(I,J,3)             &
                   + HBWGT(I,J,4)*SLP21(I,J,4)
         TS3(I,J)  = HBWGT(I,J,1)*T21(I,J,1,1)             &
                   + HBWGT(I,J,2)*T21(I,J,1,2)             &
                   + HBWGT(I,J,3)*T21(I,J,1,3)             &
                   + HBWGT(I,J,4)*T21(I,J,1,4)
         QS3(I,J)  = HBWGT(I,J,1)*Q21(I,J,1,1)             &
                   + HBWGT(I,J,2)*Q21(I,J,1,2)             &
                   + HBWGT(I,J,3)*Q21(I,J,1,3)             &
                   + HBWGT(I,J,4)*Q21(I,J,1,4)
 44      CONTINUE
       ENDDO
       ENDDO


      SP_max=-1.e20
      SP_min=1.e20
      DO J=1,NY
      DO I=1,NX
        IF(SP_max.lt.SLP1(i,j))then
          SP_max=SLP1(i,j)
        end if
        if(SP_min.gt.SLP1(i,j))then
           SP_min=SLP1(i,j)
        end if
      END DO
      END DO
      print*,'SP1_max,SP1_min=',SP_max,SP_min

      SP_max=-1.e20
      SP_min=1.e20
      DO J=1,JY
      DO I=1,JX
        IF(SP_max.lt.SLP3(i,j))then
          SP_max=SLP3(i,j)
        end if
        if(SP_min.gt.SLP3(i,j))then
           SP_min=SLP3(i,j)
        end if
      END DO
      END DO
      print*,'SP_max,SP_min=',SP_max,SP_min

!      DO J=1,JY
!      DO I=1,JX
!        IF(I.EQ.J)THEN
!          PRINT*,'I,J,SLP3,TS3,QS3=',I,J,SLP3(I,J),TS3(I,J),QS3(I,J)
!        END IF
!      END DO
!      END DO

! compute interpolation for d02

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  computer interpolation index for inner nest

      ERR=1.e20
      DO J=1,NY
      DO I=1,NX
        DIF1=abs(HLON2(1,1)-HLON1(I,J))+abs(HLAT2(1,1)-HLAT1(I,J))
        IF(DIF1.LT.ERR)THEN
          ILOC=I
          JLOC=J
          ERR=DIF1
        END IF
      END DO
      END DO

!      WBD2= WBD1 + (ILOC -1)*DLMD1
!      SBD2= SBD1 + (JLOC -1)*DPHD1
      WBD2=-((KX-1.)/2.)*DLMD2                  ! PARENT wbd
      SBD2=-((KY-1.)/2.)*DPHD2                  ! PARENT SBD

      PRINT*,'ILOC,JLOC,ERR 1=',ILOC,JLOC,ERR

      print*,'CLON2,CLAT2=',CLON2,CLAT2


!      CLON2=HLON2(1+(IX-1)/2,1+(IY-1)/2)
!      CLAT2=HLAT2(1+(IX-1)/2,1+(IY-1)/2)

!      print*,'news CLON2,CLAT2=',CLON2,CLAT2

!      WBD2=-(IX-1)*DLMD2
!      SBD2=-((IY-1)/2)*DPHD2

       print*,'DLMD2,DPHD2,WBD2,SBD2=',DLMD2,DPHD2,WBD2,SBD2
       print*,'HLON2(1,1),HLAT2(1,1)=',HLON2(1,1),HLAT2(1,1)

!      NDX=(ILOC+20)*2
!      NDY=(JLOC+40)*2
      NDX=0
      NDY=0

      IX_1=IX-1
      IY_1=IY-1

      IXDX=IX+NDX
      IYDY=IY+NDY

      IIH1=999999
      JJH1=999999
      IIV1=999999
      JJV1=999999
      HBWGT1=0.
      VBWGT1=0.

!      WBD2=WBD2-NDX*DLMD2/2                ! WBD2=WBD2-DLMD2*2*NDX/2
!      SBD2=SBD2-NDY*DPHD2/2
      WBD2=-((KX-1.)/2.)*DLMD2                  ! PARENT wbd
      SBD2=-((KY-1.)/2.)*DPHD2                  ! PARENT SBD

       CALL G2T2H_BGRID( IIH1,JJH1,            & ! output grid index and weights
                   HBWGT1,                         & ! output weights in terms of parent grid
                   HLAT3,HLON3,                    & ! target (nest) input lat lon in degrees
                   DLMD2,DPHD2,WBD2,SBD2,          & ! parent res, western and south boundaries
                   CLAT2,CLON2,                    & ! parent central lat,lon, all in degrees
                   IXDX,IYDY,                      & ! parent imax and jmax
                   JX,JY)                            ! target (nest) grid dimensions


       CALL G2T2V_BGRID( IIV1,JJV1,            & ! output grid index and weights
                   VBWGT1,                         & ! output weights in terms of parent grid
                   VLAT3,VLON3,                    & ! target (nest) input lat lon in degrees
                   DLMD2,DPHD2,WBD2,SBD2,          & ! parent res, western and south boundaries
                   CLAT2,CLON2,                    & ! parent central lat,lon, all in degrees
                   IXDX,IYDY,                      & ! parent imax and jmax
                   JX,JY)                            ! target (nest) grid dimensions

       DO J=1,JY
       DO I=1,JX
         DO N1=1,4
           IIH1(I,J,N1)=IIH1(I,J,N1)-NDX/2
           JJH1(I,J,N1)=JJH1(I,J,N1)-NDY/2
           IIV1(I,J,N1)=IIV1(I,J,N1)-NDX/2
           JJV1(I,J,N1)=JJV1(I,J,N1)-NDY/2
         END DO
       END DO
       END DO


       print*,'IIH1(1,1,1),JJH1(1,1,1)=',IIH1(1,1,1),JJH1(1,1,1)
       print*,'IIH1(JX,JY,1),JJH1(JX,JY,1)=',IIH1(JX,JY,1),JJH1(JX,JY,1)
       print*,'IIV1(1,1,1),JJV1(1,1,1)=',IIV1(1,1,1),JJV1(1,1,1)
       print*,'IIV1(JX,JY,1),JJV1(JX,JY,1)=',IIV1(JX,JY,1),JJV1(JX,JY,1)


       print*,'check coef'

       DO J=1,JY
       DO I=1,JX
         IF(IIH1(I,J,1).LE.1.or.IIH1(I,J,1).GE.IX_1)GO TO 41
         IF(JJH1(I,J,1).LE.1.or.JJH1(I,J,1).GE.IY_1)GO TO 41
         CSUM=0.
         DO N1=1,4
           CSUM=CSUM+abs(HBWGT1(I,J,N1))
         END DO
         IF(CSUM.GT.1.01.or.CSUM.lt.0.99)THEN
           print*,'test failed, I,J,CSUM=',I,J,CSUM
         END IF
 41      CONTINUE
       END DO
       END DO

       print*,'test66'
    
! replace the topgraphy data for the inner nest

       A13=1.                 ! all ocean
       B13=0.002
       C13=0.88

         DO J=1,JY
         DO I=1,JX
           IF(IIH1(I,J,1).LE.1.or.IIH1(I,J,1).GE.IX_1)GO TO 34
           IF(JJH1(I,J,1).LE.1.or.JJH1(I,J,1).GE.IY_1)GO TO 34
           DO N1=1,4
             IH1(N1)=IIH1(I,J,N1)
             JH1(N1)=JJH1(I,J,N1)
           END DO
           ZS3(I,J) = HBWGT1(I,J,1)*Z2(IH1(1),JH1(1),1)             &
                    + HBWGT1(I,J,2)*Z2(IH1(2),JH1(2),1)             &
                    + HBWGT1(I,J,3)*Z2(IH1(3),JH1(3),1)             &
                    + HBWGT1(I,J,4)*Z2(IH1(4),JH1(4),1)
 34      CONTINUE
         END DO
         END DO

       DO J=1,JY
       DO I=1,JX
         Z3(I,J,1)=ZS3(I,J)
       END DO
       END DO
             
!!!!!!!!!!!!!

       DO J=1,JY
       DO I=1,JX
         IF(IIH1(I,J,1).LE.1.or.IIH1(I,J,1).GE.IX_1)GO TO 45
         IF(JJH1(I,J,1).LE.1.or.JJH1(I,J,1).GE.IY_1)GO TO 45
            DO N1=1,4
              IH1(N1)=IIH1(I,J,N1)
              JH1(N1)=JJH1(I,J,N1)
            END DO
         DO N1=1,4
           SLP21(I,J,N1)=SLP2(IH1(N1),JH1(N1))
         END DO
         K=1             ! surface
         DO N1=1,4
           IF(Z3(I,J,K).LT.ZMID2(IH1(N1),JH1(N1),1))THEN
             DZ1=ZMID2(IH1(N1),JH1(N1),1)-Z3(I,J,K)
             T21(I,J,K,N1)=T2(IH1(N1),JH1(N1),1)+GAMMA*DZ1
             Q21(I,J,K,N1)=Q2(IH1(N1),JH1(N1),1)
           ELSE IF(Z3(I,J,K).GE.ZMID2(IH1(N1),JH1(N1),NZ))THEN  ! never occur for K=1
             DZ1=ZMID2(IH1(N1),JH1(N1),NZ)-Z3(I,J,K)
             T21(I,J,K,N1)=T2(IH1(N1),JH1(N1),NZ)+GAMMA*DZ1
             Q21(I,J,K,N1)=Q2(IH1(N1),JH1(N1),NZ)
           ELSE
             DO K1=2,IZ
               ZDIF1=ZMID2(IH1(N1),JH1(N1),K1)-Z3(I,J,K)
               IF(ZDIF1.GE.0.)THEN
                 FACT1=ZDIF1/(ZMID2(IH1(N1),JH1(N1),K1)-   &
                              ZMID2(IH1(N1),JH1(N1),K1-1))
                 T21(I,J,K,N1)=T2(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +T2(IH1(N1),JH1(N1),K1-1)*FACT1
                 Q21(I,J,K,N1)=Q2(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +Q2(IH1(N1),JH1(N1),K1-1)*FACT1
                 GO TO 39
               END IF
             END DO
           END IF
 39        CONTINUE
         END DO

 45      CONTINUE
       ENDDO
       ENDDO
       
       print*,'test67'

       DO J=1,JY
       DO I=1,JX
         IF(IIH1(I,J,1).LE.1.or.IIH1(I,J,1).GE.IX_1)GO TO 47
         IF(JJH1(I,J,1).LE.1.or.JJH1(I,J,1).GE.IY_1)GO TO 47
         SLP3(I,J) = HBWGT1(I,J,1)*SLP21(I,J,1)             &
                   + HBWGT1(I,J,2)*SLP21(I,J,2)             &
                   + HBWGT1(I,J,3)*SLP21(I,J,3)             &
                   + HBWGT1(I,J,4)*SLP21(I,J,4)
         TS3(I,J)  = HBWGT1(I,J,1)*T21(I,J,1,1)             &
                   + HBWGT1(I,J,2)*T21(I,J,1,2)             &
                   + HBWGT1(I,J,3)*T21(I,J,1,3)             &
                   + HBWGT1(I,J,4)*T21(I,J,1,4)
         QS3(I,J)  = HBWGT1(I,J,1)*Q21(I,J,1,1)             &
                   + HBWGT1(I,J,2)*Q21(I,J,1,2)             &
                   + HBWGT1(I,J,3)*Q21(I,J,1,3)             &
                   + HBWGT1(I,J,4)*Q21(I,J,1,4)
 47      CONTINUE
       ENDDO
       ENDDO

!??????????????????

      IF(IFLAG_NEST.EQ.2)THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  interpolate from d03
!
!      NDY=(JLOC+40)*2

      NDX=0
      NDY=0
      NXDX=NX5+NDX
      NYDY=NY5+NDY

      WBD5=WBD5-NDX*DLMD5/2                ! WBD2=WBD2-DLMD2*2*NDX/2
      SBD5=SBD5-NDY*DPHD5/2

      IIH2=999999
      JJH2=999999
      IIV2=999999
      JJV2=999999
      HBWGT2=0.
      VBWGT2=0.

       CALL G2T2H_BGRID( IIH2,JJH2,            & ! output grid index and weights
                   HBWGT2,                         & ! output weights in terms of parent grid
                   HLAT3,HLON3,                    & ! target (nest) input lat lon in degrees
                   DLMD5,DPHD5,WBD5,SBD5,          & ! parent res, western and south boundaries
                   CLAT5,CLON5,                    & ! parent central lat,lon, all in degrees
                   NXDX,NYDY,                      & ! parent imax and jmax
                   JX,JY)                            ! target (nest) grid dimensions


       CALL G2T2V_BGRID( IIV2,JJV2,            & ! output grid index and weights
                   VBWGT2,                         & ! output weights in terms of parent grid
                   VLAT3,VLON3,                    & ! target (nest) input lat lon in degrees
                   DLMD5,DPHD5,WBD5,SBD5,          & ! parent res, western and south boundaries
                   CLAT5,CLON5,                    & ! parent central lat,lon, all in degrees
                   NXDX,NYDY,                      & ! parent imax and jmax
                   JX,JY)                            ! target (nest) grid dimensions

       DO J=1,JY
       DO I=1,JX
         DO N1=1,4
           IIH2(I,J,N1)=IIH2(I,J,N1)-NDX/2
           JJH2(I,J,N1)=JJH2(I,J,N1)-NDY/2
           IIV2(I,J,N1)=IIV2(I,J,N1)-NDX/2
           JJV2(I,J,N1)=JJV2(I,J,N1)-NDY/2
         END DO
       END DO
       END DO


       NX5_1=NX5-1
       NY5_1=NY5-1

       print*,'check coef'

       DO J=1,JY
       DO I=1,JX
         IF(IIH2(I,J,1).LE.1.or.IIH2(I,J,1).GE.NX5_1)GO TO 441
         IF(JJH2(I,J,1).LE.1.or.JJH2(I,J,1).GE.NY5_1)GO TO 441
         CSUM=0.
         DO N1=1,4
           CSUM=CSUM+abs(HBWGT2(I,J,N1))
         END DO
         IF(CSUM.GT.1.01.or.CSUM.lt.0.99)THEN
           print*,'test failed, I,J,CSUM=',I,J,CSUM
         END IF
 441      CONTINUE
       END DO
       END DO

       print*,'test669'
    
! replace the topgraphy data for the inner nest

       A13=1.                 ! all ocean
       B13=0.002
       C13=0.88

       IF(IBGS.EQ.1)THEN

         DO J=1,NY5
         DO I=1,NX5
           J1=J+J_ST14-1
           I1=I+I_ST14-1
           ZS3(I1,J1)=Z5(I,J,1)
           A13(I1,J1)=A105(I,J)
           B13(I1,J1)=B105(I,J)
           C13(I1,J1)=C105(I,J)
         END DO
         END DO

       ELSE IF(IBGS.EQ.0)THEN

         DO J=1,JY
         DO I=1,JX
           IF(IIH2(I,J,1).LE.1.or.IIH2(I,J,1).GE.NX5_1)GO TO 344
           IF(JJH2(I,J,1).LE.1.or.JJH2(I,J,1).GE.NY5_1)GO TO 344
           DO N1=1,4
             IH1(N1)=IIH2(I,J,N1)
             JH1(N1)=JJH2(I,J,N1)
           END DO
           ZS3(I,J) = HBWGT2(I,J,1)*Z5(IH1(1),JH1(1),1)             &
                    + HBWGT2(I,J,2)*Z5(IH1(2),JH1(2),1)             &
                    + HBWGT2(I,J,3)*Z5(IH1(3),JH1(3),1)             &
                    + HBWGT2(I,J,4)*Z5(IH1(4),JH1(4),1)
 344      CONTINUE
         END DO
         END DO

       END IF
 
         DO J=1,JY
         DO I=1,JX
           Z3(I,J,1)=ZS3(I,J)
         END DO
         END DO
            
!!!!!!!!!!!!!

       DO J=1,JY
       DO I=1,JX
         IF(IIH2(I,J,1).LE.1.or.IIH2(I,J,1).GE.NX5_1)GO TO 445
         IF(JJH2(I,J,1).LE.1.or.JJH2(I,J,1).GE.NY5_1)GO TO 445
            DO N1=1,4
              IH1(N1)=IIH2(I,J,N1)
              JH1(N1)=JJH2(I,J,N1)
            END DO
         DO N1=1,4
           SLP21(I,J,N1)=SLP5(IH1(N1),JH1(N1))
         END DO
         K=1             ! surface
         DO N1=1,4
           IF(Z3(I,J,K).LT.ZMID5(IH1(N1),JH1(N1),1))THEN
             DZ1=ZMID5(IH1(N1),JH1(N1),1)-Z3(I,J,K)
             T21(I,J,K,N1)=T5(IH1(N1),JH1(N1),1)+GAMMA*DZ1
             Q21(I,J,K,N1)=Q5(IH1(N1),JH1(N1),1)
           ELSE IF(Z3(I,J,K).GE.ZMID5(IH1(N1),JH1(N1),NZ))THEN  ! never occur for K=1
             DZ1=ZMID5(IH1(N1),JH1(N1),NZ)-Z3(I,J,K)
             T21(I,J,K,N1)=T5(IH1(N1),JH1(N1),NZ)+GAMMA*DZ1
             Q21(I,J,K,N1)=Q5(IH1(N1),JH1(N1),NZ)
           ELSE
             DO K1=2,IZ
               ZDIF1=ZMID5(IH1(N1),JH1(N1),K1)-Z3(I,J,K)
               IF(ZDIF1.GE.0.)THEN
                 FACT1=ZDIF1/(ZMID5(IH1(N1),JH1(N1),K1)-   &
                              ZMID5(IH1(N1),JH1(N1),K1-1))
                 T21(I,J,K,N1)=T5(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +T5(IH1(N1),JH1(N1),K1-1)*FACT1
                 Q21(I,J,K,N1)=Q5(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +Q5(IH1(N1),JH1(N1),K1-1)*FACT1
                 GO TO 339
               END IF
             END DO
           END IF
 339        CONTINUE
         END DO

 445      CONTINUE
       ENDDO
       ENDDO
       
       print*,'test67'

       DO J=1,JY
       DO I=1,JX
         IF(IIH2(I,J,1).LE.1.or.IIH2(I,J,1).GE.NX5_1)GO TO 447
         IF(JJH2(I,J,1).LE.1.or.JJH2(I,J,1).GE.NY5_1)GO TO 447
         SLP3(I,J) = HBWGT2(I,J,1)*SLP21(I,J,1)             &
                   + HBWGT2(I,J,2)*SLP21(I,J,2)             &
                   + HBWGT2(I,J,3)*SLP21(I,J,3)             &
                   + HBWGT2(I,J,4)*SLP21(I,J,4)
         TS3(I,J)  = HBWGT2(I,J,1)*T21(I,J,1,1)             &
                   + HBWGT2(I,J,2)*T21(I,J,1,2)             &
                   + HBWGT2(I,J,3)*T21(I,J,1,3)             &
                   + HBWGT2(I,J,4)*T21(I,J,1,4)
         QS3(I,J)  = HBWGT2(I,J,1)*Q21(I,J,1,1)             &
                   + HBWGT2(I,J,2)*Q21(I,J,1,2)             &
                   + HBWGT2(I,J,3)*Q21(I,J,1,3)             &
                   + HBWGT2(I,J,4)*Q21(I,J,1,4)
 447      CONTINUE
       ENDDO
       ENDDO


      END IF     ! IFLAG_NEST

!??????????????????

       print*,'test68'

! 
! Construct 3D pressure grid

       DO J=1,JY
       DO I=1,JX
         ZSFC = ZS3(I,J)
         TSFC = TS3(I,J)*(1.+D608*QS3(I,J))
         A = (GAMMA * ZSFC) / TSFC
         P3(I,J,1) = SLP3(I,J)/(1+A)**COEF2
!         PD3(I,J)=P3(I,J,1)-PDTOP3-PT3
         PD3(I,J)=P3(I,J,1)-PT3     ! for MMB
       ENDDO
       ENDDO
                                                                                                                                      
! PD(I,J)=P1(I,J,1)-PDTOP-PT=PSFC(I,J)-PDTOP-PT
       P3(:,:,NZ1)=PT3
       DO K=1,NZ
       DO J=1,JY
       DO I=1,JX
!         P3(I,J,NZ1-K)=PT3+PDTOP3*DSG1(K)+PD3(I,J)*DSG2(K)     ! PD(I,J) changed
         P3(I,J,NZ1-K)=P3(I,J,NZ1-K+1)+PDTOP3*DSG1(K)+PD3(I,J)*DSG2(K)     ! PD(I,J) changed
       ENDDO
       ENDDO
       ENDDO

      print*,'check pressure calculation'

      print*,'P3(9,9,NZ1)=',P3(9,9,NZ1)
      do k=1,NZ
        print*,K,P3(9,9,K),DSG1(K),DSG2(K)
      end do

       DO K=1,NZ
       DO J=1,JY
       DO I=1,JX
         PMID3(I,J,K)=EXP((ALOG(P3(I,J,K))+ALOG(P3(I,J,K+1)))*0.5)
       ENDDO
       ENDDO
       ENDDO

! interpolate vertically to 3D-P level in new coordinate  (H Points)
! from outer nest data

       DO J=1,JY
       DO I=1,JX
         IF(IIH(I,J,1).LE.2.or.IIH(I,J,1).GE.NX_1)GO TO 70
         IF(JJH(I,J,1).LE.2.or.JJH(I,J,1).GE.NY_1)GO TO 70
         DO N1=1,4
           IH1(N1)=IIH(I,J,N1)
           JH1(N1)=JJH(I,J,N1)
         END DO
!         DO N1=1,4
!           SLP21(I,J,N1)=SLP2(IH1(N1),JH1(N1))
!         END DO
         DO K=1,NZ
         DO N1=1,4
           IF(PMID3(I,J,K).GT.PMID1(IH1(N1),JH1(N1),1))THEN
             DZ1=T1(IH1(N1),JH1(N1),1)      &
                 /GAMMA*(1.-(PMID3(I,J,K)/PMID1(IH1(N1),JH1(N1),1))**COEF3)
             T21(I,J,K,N1)=T1(IH1(N1),JH1(N1),1)-GAMMA*DZ1
             Q21(I,J,K,N1)=Q1(IH1(N1),JH1(N1),1)
           ELSE IF(PMID3(I,J,K).LE.PMID1(IH1(N1),JH1(N1),NZ))THEN
             DZ1=T1(IH1(N1),JH1(N1),NZ)      &
                 /GAMMA*(1.-(PMID3(I,J,K)/PMID1(IH1(N1),JH1(N1),NZ))**COEF3)
             T21(I,J,K,N1)=T1(IH1(N1),JH1(N1),NZ)-GAMMA*DZ1
             Q21(I,J,K,N1)=Q1(IH1(N1),JH1(N1),NZ)
           ELSE
             DO K1=2,NZ
               PDIF1=ALOG(PMID3(I,J,K))-ALOG(PMID1(IH1(N1),JH1(N1),K1))
               IF(PDIF1.GE.0.)THEN
                 FACT1=PDIF1/(ALOG(PMID1(IH1(N1),JH1(N1),K1-1))-   &
                               ALOG(PMID1(IH1(N1),JH1(N1),K1)))
                 T21(I,J,K,N1)=T1(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +T1(IH1(N1),JH1(N1),K1-1)*FACT1
                 Q21(I,J,K,N1)=Q1(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +Q1(IH1(N1),JH1(N1),K1-1)*FACT1
                 GO TO 50
               END IF
             END DO
           END IF
 50        CONTINUE
         END DO
         END DO
 70      CONTINUE
       ENDDO
       ENDDO

       print*,'test69'
!
       DO J=1,JY
       DO I=1,JX
         IF(IIH(I,J,1).LE.2.or.IIH(I,J,1).GE.NX_1)GO TO 120
         IF(JJH(I,J,1).LE.2.or.JJH(I,J,1).GE.NY_1)GO TO 120
         DO K=1,NZ
            T3(I,J,K) =                               &
                HBWGT(I,J,1)*T21(I,J,K,1)             &
              + HBWGT(I,J,2)*T21(I,J,K,2)             &
              + HBWGT(I,J,3)*T21(I,J,K,3)             &
              + HBWGT(I,J,4)*T21(I,J,K,4)
            Q3(I,J,K) =                               &
                HBWGT(I,J,1)*Q21(I,J,K,1)             &
              + HBWGT(I,J,2)*Q21(I,J,K,2)             &
              + HBWGT(I,J,3)*Q21(I,J,K,3)             &
              + HBWGT(I,J,4)*Q21(I,J,K,4)
         ENDDO
 120     CONTINUE
       ENDDO
       ENDDO

       print*,'test70'

! interpolate vertically to 3D-P level in new coordinate  (H Points)
! from inner nest data

       DO J=1,JY
       DO I=1,JX
         IF(IIH1(I,J,1).LE.1.or.IIH1(I,J,1).GE.IX_1)GO TO 73
         IF(JJH1(I,J,1).LE.1.or.JJH1(I,J,1).GE.IY_1)GO TO 73
         DO N1=1,4
           IH1(N1)=IIH1(I,J,N1)
           JH1(N1)=JJH1(I,J,N1)
         END DO
!         DO N1=1,4
!           SLP21(I,J,N1)=SLP2(IH1(N1),JH1(N1))
!         END DO
         DO K=1,NZ
         DO N1=1,4
           IF(PMID3(I,J,K).GT.PMID2(IH1(N1),JH1(N1),1))THEN
             DZ1=T2(IH1(N1),JH1(N1),1)      &
                 /GAMMA*(1.-(PMID3(I,J,K)/PMID2(IH1(N1),JH1(N1),1))**COEF3)
             T21(I,J,K,N1)=T2(IH1(N1),JH1(N1),1)-GAMMA*DZ1
             Q21(I,J,K,N1)=Q2(IH1(N1),JH1(N1),1)
           ELSE IF(PMID3(I,J,K).LE.PMID2(IH1(N1),JH1(N1),NZ))THEN
             DZ1=T2(IH1(N1),JH1(N1),NZ)      &
                 /GAMMA*(1.-(PMID3(I,J,K)/PMID2(IH1(N1),JH1(N1),NZ))**COEF3)
             T21(I,J,K,N1)=T2(IH1(N1),JH1(N1),NZ)-GAMMA*DZ1
             Q21(I,J,K,N1)=Q2(IH1(N1),JH1(N1),NZ)
           ELSE
             DO K1=2,NZ
               PDIF1=ALOG(PMID3(I,J,K))-ALOG(PMID2(IH1(N1),JH1(N1),K1))
               IF(PDIF1.GE.0.)THEN
                 FACT1=PDIF1/(ALOG(PMID2(IH1(N1),JH1(N1),K1-1))-   &
                      ALOG(PMID2(IH1(N1),JH1(N1),K1)))
                 T21(I,J,K,N1)=T2(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +T2(IH1(N1),JH1(N1),K1-1)*FACT1
                 Q21(I,J,K,N1)=Q2(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +Q2(IH1(N1),JH1(N1),K1-1)*FACT1
                 GO TO 53
               END IF
             END DO
           END IF
 53        CONTINUE
         END DO
         END DO
 73      CONTINUE
       ENDDO
       ENDDO

       DO J=1,JY
       DO I=1,JX
         IF(IIH1(I,J,1).LE.1.or.IIH1(I,J,1).GE.IX_1)GO TO 123
         IF(JJH1(I,J,1).LE.1.or.JJH1(I,J,1).GE.IY_1)GO TO 123
         DO K=1,NZ
            T3(I,J,K) =                               &
                HBWGT1(I,J,1)*T21(I,J,K,1)             &
              + HBWGT1(I,J,2)*T21(I,J,K,2)             &
              + HBWGT1(I,J,3)*T21(I,J,K,3)             &
              + HBWGT1(I,J,4)*T21(I,J,K,4)             
            Q3(I,J,K) =                               &
                HBWGT1(I,J,1)*Q21(I,J,K,1)             &
              + HBWGT1(I,J,2)*Q21(I,J,K,2)             &
              + HBWGT1(I,J,3)*Q21(I,J,K,3)             &
              + HBWGT1(I,J,4)*Q21(I,J,K,4)             
         ENDDO
 123     CONTINUE
       ENDDO
       ENDDO

       IF(IFLAG_NEST.EQ.2)THEN

! interpolate from d03

       DO J=1,JY
       DO I=1,JX
         IF(IIH2(I,J,1).LE.1.or.IIH2(I,J,1).GE.NX5_1)GO TO 773
         IF(JJH2(I,J,1).LE.1.or.JJH2(I,J,1).GE.NY5_1)GO TO 773
         DO N1=1,4
           IH1(N1)=IIH2(I,J,N1)
           JH1(N1)=JJH2(I,J,N1)
         END DO
!         DO N1=1,4
!           SLP21(I,J,N1)=SLP5(IH1(N1),JH1(N1))
!         END DO
         DO K=1,NZ
         DO N1=1,4
           IF(PMID3(I,J,K).GT.PMID5(IH1(N1),JH1(N1),1))THEN
             DZ1=T5(IH1(N1),JH1(N1),1)      &
                 /GAMMA*(1.-(PMID3(I,J,K)/PMID5(IH1(N1),JH1(N1),1))**COEF3)
             T21(I,J,K,N1)=T5(IH1(N1),JH1(N1),1)-GAMMA*DZ1
             Q21(I,J,K,N1)=Q5(IH1(N1),JH1(N1),1)
           ELSE IF(PMID3(I,J,K).LE.PMID5(IH1(N1),JH1(N1),NZ))THEN
             DZ1=T5(IH1(N1),JH1(N1),NZ)      &
                 /GAMMA*(1.-(PMID3(I,J,K)/PMID5(IH1(N1),JH1(N1),NZ))**COEF3)
             T21(I,J,K,N1)=T5(IH1(N1),JH1(N1),NZ)-GAMMA*DZ1
             Q21(I,J,K,N1)=Q5(IH1(N1),JH1(N1),NZ)
           ELSE
             DO K1=2,NZ
               PDIF1=ALOG(PMID3(I,J,K))-ALOG(PMID5(IH1(N1),JH1(N1),K1))
               IF(PDIF1.GE.0.)THEN
                 FACT1=PDIF1/(ALOG(PMID5(IH1(N1),JH1(N1),K1-1))-   &
                      ALOG(PMID5(IH1(N1),JH1(N1),K1)))
                 T21(I,J,K,N1)=T5(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +T5(IH1(N1),JH1(N1),K1-1)*FACT1
                 Q21(I,J,K,N1)=Q5(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                              +Q5(IH1(N1),JH1(N1),K1-1)*FACT1
                 GO TO 553
               END IF
             END DO
           END IF
 553        CONTINUE
         END DO
         END DO
 773      CONTINUE
       ENDDO
       ENDDO

       DO J=1,JY
       DO I=1,JX
         IF(IIH2(I,J,1).LE.1.or.IIH2(I,J,1).GE.NX5_1)GO TO 223
         IF(JJH2(I,J,1).LE.1.or.JJH2(I,J,1).GE.NY5_1)GO TO 223
         DO K=1,NZ
            T3(I,J,K) =                               &
                HBWGT1(I,J,1)*T21(I,J,K,1)             &
              + HBWGT1(I,J,2)*T21(I,J,K,2)             &
              + HBWGT1(I,J,3)*T21(I,J,K,3)             &
              + HBWGT1(I,J,4)*T21(I,J,K,4)             
            Q3(I,J,K) =                               &
                HBWGT1(I,J,1)*Q21(I,J,K,1)             &
              + HBWGT1(I,J,2)*Q21(I,J,K,2)             &
              + HBWGT1(I,J,3)*Q21(I,J,K,3)             &
              + HBWGT1(I,J,4)*Q21(I,J,K,4)             
         ENDDO
 223     CONTINUE
       ENDDO
       ENDDO

      END IF     ! IFLAG_NEST

       print*,'test71'
 
! Compute Geopotentital height, INTEGRATE HEIGHT HYDROSTATICLY
            
       DO J=1,JY
       DO I=1,JX
         ZSFC = ZS3(I,J)
         TSFC = T3(I,J,1)*(1.+D608*Q3(I,J,1))
         A = (GAMMA * ZSFC) / TSFC
         P3(I,J,1) = SLP3(I,J)/(1+A)**COEF2
!         PD3(I,J)=P3(I,J,1)-PDTOP3-PT3
         PD3(I,J)=P3(I,J,1)-PT3
       ENDDO
       ENDDO

! PD(I,J)=P1(I,J,1)-PDTOP-PT=PSFC(I,J)-PDTOP-PT
       P3(:,:,NZ1)=PT3
       DO K=1,NZ
       DO J=1,JY
       DO I=1,JX
!         P3(I,J,NZ1-K)=PT3+PDTOP3*DSG1(K)+PD3(I,J)*DSG2(K)     ! PD(I,J) changed
         P3(I,J,NZ1-K)=P3(I,J,NZ1-K+1)+PDTOP3*DSG1(K)+PD3(I,J)*DSG2(K)     ! PD(I,J) changed
       ENDDO
       ENDDO
       ENDDO

      do j = 1,JY
      do i = 1,JX
        Z3(I,J,1)=ZS3(I,J)
        DO L=2,NZ+1
          Z3(I,J,L)=Z3(I,J,L-1)+T3(I,J,L-1)*          &
              (Q3(I,J,L-1)*0.608+1.0)*287.04*         &
              (ALOG(P3(I,J,L-1))-ALOG(P3(I,J,L)))/G
        ENDDO
       ENDDO
      END DO

       DO K=1,NZ
       DO J=1,JY
       DO I=1,JX
         PMID3(I,J,K)=EXP((ALOG(P3(I,J,K))+ALOG(P3(I,J,K+1)))*0.5)
       ENDDO
       ENDDO
       ENDDO

       PMV3=PMID3

       print*,'test711'

! interpolate vertically to P level in new coordinate  (V Points)
       DO J=2,JY-1
       DO K=1,NZ
       DO I=2,JX-1
         PMV3(I,J,K)=0.25*(PMID3(I,J,K)+PMID3(I+1,J,K)+            &
                           PMID3(I,J+1,K)+PMID3(I+1,J+1,K))
       END DO
       END DO
       END DO

       PMV1=PMID1     

       DO J=2,NY-1
       DO I=2,NX-1
       DO K=1,NZ
         PMV1(I,J,K)=0.25*(PMID1(I,J,K)+PMID1(I+1,J,K)+            &
                           PMID1(I,J+1,K)+PMID1(I+1,J+1,K))
       END DO
       END DO
       END DO
 
       print*,'test712'

       DO J=1,JY
       DO I=1,JX
         IF(IIV(I,J,1).LE.2.or.IIV(I,J,1).GE.NX_1)GO TO 80
         IF(JJV(I,J,1).LE.2.or.JJV(I,J,1).GE.NY_1)GO TO 80
         DO N1=1,4
           IV1(N1)=IIV(I,J,N1)
           JV1(N1)=JJV(I,J,N1)
         END DO
!         A13(I,J) = VBWGT(I,J,1)*A101(IV1(1),JV1(1))             &
!                  + VBWGT(I,J,2)*A101(IV1(2),JV1(2))             &
!                  + VBWGT(I,J,3)*A101(IV1(3),JV1(3))             &
!                  + VBWGT(I,J,4)*A101(IV1(4),JV1(4))

         DO K=1,NZ
           DO N1=1,4
             IF(PMV3(I,J,K).GT.PMV1(IV1(N1),JV1(N1),1))THEN
               DP1=PMV3(I,J,K)-PMV1(IV1(N1),JV1(N1),1)
!               U21(I,J,K,N1)=U1(IV1(N1),JV1(N1),1)*(1.-DP1*1.4E-5)
!               V21(I,J,K,N1)=V1(IV1(N1),JV1(N1),1)*(1.-DP1*1.4E-5)
               U21(I,J,K,N1)=U1(IV1(N1),JV1(N1),1)
               V21(I,J,K,N1)=V1(IV1(N1),JV1(N1),1)
             ELSE IF(PMV3(I,J,K).LE.PMV1(IV1(N1),JV1(N1),NZ))THEN
               U21(I,J,K,N1)=U1(IV1(N1),JV1(N1),NZ)
               V21(I,J,K,N1)=V1(IV1(N1),JV1(N1),NZ)
             ELSE
               DO K1=2,NZ
                 PDIF1=ALOG(PMV3(I,J,K))-ALOG(PMV1(IV1(N1),JV1(N1),K1))
                 IF(PDIF1.GE.0.)THEN
                   FACT1=PDIF1/(ALOG(PMV1(IV1(N1),JV1(N1),K1-1))-             &
                                ALOG(PMV1(IV1(N1),JV1(N1),K1)))
                   U21(I,J,K,N1)=U1(IV1(N1),JV1(N1),K1)*(1.-FACT1)     &
                                +U1(IV1(N1),JV1(N1),K1-1)*FACT1
                   V21(I,J,K,N1)=V1(IV1(N1),JV1(N1),K1)*(1.-FACT1)     &
                                +V1(IV1(N1),JV1(N1),K1-1)*FACT1
                   GO TO 60
                 END IF
               END DO
             END IF
 60          CONTINUE
           END DO
         END DO
 80      CONTINUE
       END DO
       END DO

       print*,'test72'

!23456789012345678901234567890123456789012345678901234567890123456789012

       DO J=1,JY
       DO I=1,JX
         IF(IIV(I,J,1).LE.2.or.IIV(I,J,1).GE.NX_1)GO TO 110
         IF(JJV(I,J,1).LE.2.or.JJV(I,J,1).GE.NY_1)GO TO 110
         DO K=1,NZ
            U3(I,J,K) =                       & 
               VBWGT(I,J,1)*U21(I,J,K,1)      &
             + VBWGT(I,J,2)*U21(I,J,K,2)      &
             + VBWGT(I,J,3)*U21(I,J,K,3)      &
             + VBWGT(I,J,4)*U21(I,J,K,4)
            V3(I,J,K) =                       &
               VBWGT(I,J,1)*V21(I,J,K,1)      &
             + VBWGT(I,J,2)*V21(I,J,K,2)      &
             + VBWGT(I,J,3)*V21(I,J,K,3)      &
             + VBWGT(I,J,4)*V21(I,J,K,4)
         END DO
!         if(J.GT.0.95*JY.and.I.gt.0.95*JX)print*,'test72111',I,J
 110     CONTINUE
       ENDDO
       ENDDO

        print*,'JX,JY,NZ=',JX,JY,NZ
        write(*,*)'interpolation after d01='
      do k=1,NZ
        write(*,32)K,T3(9,9,K),            &
          Q3(9,9,K),U3(9,9,K),V3(9,9,K),Z3(9,9,K),P3(9,9,K)
      end do

       print*,'test721'

!      WRITE(64)((SLP3(I,J),I=1,JX),J=1,JY)
!      DO K=1,NZ+1
!        WRITE(64)((Z3(I,J,K),I=1,JX),J=1,JY)
!      END DO
!      DO K=1,NZ+1
!        WRITE(64)((P3(I,J,K),I=1,JX),J=1,JY)
!      END DO
!      DO K=1,NZ
!       WRITE(64)((T3(I,J,K),I=1,JX),J=1,JY)
!     END DO
!     DO K=1,NZ
!       WRITE(64)((Q3(I,J,K),I=1,JX),J=1,JY)
!     END DO
!     DO K=1,NZ
!       WRITE(64)((U3(I,J,K),I=1,JX),J=1,JY)
!     END DO
!     DO K=1,NZ
!       WRITE(64)((V3(I,J,K),I=1,JX),J=1,JY)
!     END DO

       print*,'test73'

       PMV2=PMID2     

       DO J=2,IY-1
       DO I=2,IX-1
       DO K=1,IZ
         PMV2(I,J,K)=0.25*(PMID2(I,J,K)+PMID2(I+1,J,K)+            &
                           PMID2(I,J+1,K)+PMID2(I+1,J+1,K))
       END DO
       END DO
       END DO
 
       IF(IBGS.EQ.1)THEN
       DO J=1,JY
       DO I=1,JX
         IF(IIV1(I,J,1).LE.1.or.IIV1(I,J,1).GE.IX_1)GO TO 37
         IF(JJV1(I,J,1).LE.1.or.JJV1(I,J,1).GE.IY_1)GO TO 37
         DO N1=1,4
           IV1(N1)=IIV1(I,J,N1)
           JV1(N1)=JJV1(I,J,N1)
         END DO
         A13(I,J) = VBWGT1(I,J,1)*A102(IV1(1),JV1(1))             &
                  + VBWGT1(I,J,2)*A102(IV1(2),JV1(2))             &
                  + VBWGT1(I,J,3)*A102(IV1(3),JV1(3))             &
                  + VBWGT1(I,J,4)*A102(IV1(4),JV1(4))
         C13(I,J) = VBWGT1(I,J,1)*C102(IV1(1),JV1(1))             &
                  + VBWGT1(I,J,2)*C102(IV1(2),JV1(2))             &
                  + VBWGT1(I,J,3)*C102(IV1(3),JV1(3))             &
                  + VBWGT1(I,J,4)*C102(IV1(4),JV1(4))

 37      CONTINUE
       END DO
       END DO
       ELSE IF(IBGS.EQ.0)THEN
       DO J=1,JY
       DO I=1,JX
         IF(IIV1(I,J,1).LE.1.or.IIV1(I,J,1).GE.IX_1)GO TO 43
         IF(JJV1(I,J,1).LE.1.or.JJV1(I,J,1).GE.IY_1)GO TO 43
         DO N1=1,4
           IV1(N1)=IIV1(I,J,N1)
           JV1(N1)=JJV1(I,J,N1)
         END DO
         A13(I,J) = VBWGT1(I,J,1)*A102(IV1(1),JV1(1))             &
                  + VBWGT1(I,J,2)*A102(IV1(2),JV1(2))             &
                  + VBWGT1(I,J,3)*A102(IV1(3),JV1(3))             &
                  + VBWGT1(I,J,4)*A102(IV1(4),JV1(4))
         B13(I,J) = VBWGT1(I,J,1)*B102(IV1(1),JV1(1))             &
                  + VBWGT1(I,J,2)*B102(IV1(2),JV1(2))             &
                  + VBWGT1(I,J,3)*B102(IV1(3),JV1(3))             &
                  + VBWGT1(I,J,4)*B102(IV1(4),JV1(4))
         C13(I,J) = VBWGT1(I,J,1)*C102(IV1(1),JV1(1))             &
                  + VBWGT1(I,J,2)*C102(IV1(2),JV1(2))             &
                  + VBWGT1(I,J,3)*C102(IV1(3),JV1(3))             &
                  + VBWGT1(I,J,4)*C102(IV1(4),JV1(4))

 43      CONTINUE
       END DO
       END DO
       END IF

       DO J=1,JY
       DO I=1,JX
         IF(IIV1(I,J,1).LE.1.or.IIV1(I,J,1).GE.IX_1)GO TO 85
         IF(JJV1(I,J,1).LE.1.or.JJV1(I,J,1).GE.IY_1)GO TO 85
         DO N1=1,4
           IV1(N1)=IIV1(I,J,N1)
           JV1(N1)=JJV1(I,J,N1)
         END DO
         DO K=1,NZ
           DO N1=1,4
             IF(PMV3(I,J,K).GT.PMV2(IV1(N1),JV1(N1),1))THEN
               DP1=PMV3(I,J,K)-PMV2(IV1(N1),JV1(N1),1)
!               U21(I,J,K,N1)=U2(IV1(N1),JV1(N1),1)*(1.-DP1*1.4E-5)
!               V21(I,J,K,N1)=V2(IV1(N1),JV1(N1),1)*(1.-DP1*1.4E-5)
               U21(I,J,K,N1)=U2(IV1(N1),JV1(N1),1)
               V21(I,J,K,N1)=V2(IV1(N1),JV1(N1),1)
             ELSE IF(PMV3(I,J,K).LE.PMV2(IV1(N1),JV1(N1),NZ))THEN
               U21(I,J,K,N1)=U2(IV1(N1),JV1(N1),NZ)
               V21(I,J,K,N1)=V2(IV1(N1),JV1(N1),NZ)
             ELSE
               DO K1=2,NZ
                 PDIF1=ALOG(PMV3(I,J,K))-ALOG(PMV2(IV1(N1),JV1(N1),K1))
                 IF(PDIF1.GE.0.)THEN
                   FACT1=PDIF1/(ALOG(PMV2(IV1(N1),JV1(N1),K1-1))-             &
                                ALOG(PMV2(IV1(N1),JV1(N1),K1)))
                   U21(I,J,K,N1)=U2(IV1(N1),JV1(N1),K1)*(1.-FACT1)     &
                                +U2(IV1(N1),JV1(N1),K1-1)*FACT1
                   V21(I,J,K,N1)=V2(IV1(N1),JV1(N1),K1)*(1.-FACT1)     &
                                +V2(IV1(N1),JV1(N1),K1-1)*FACT1
                   GO TO 65
                 END IF
               END DO
             END IF
 65          CONTINUE
           END DO
         END DO
 85      CONTINUE
       END DO
       END DO
!23456789012345678901234567890123456789012345678901234567890123456789012

       DO J=1,JY
       DO I=1,JX
         IF(IIV1(I,J,1).LE.1.or.IIV1(I,J,1).GE.IX_1)GO TO 115
         IF(JJV1(I,J,1).LE.1.or.JJV1(I,J,1).GE.IY_1)GO TO 115
         DO K=1,NZ
            U3(I,J,K) =                       & 
               VBWGT1(I,J,1)*U21(I,J,K,1)      &
             + VBWGT1(I,J,2)*U21(I,J,K,2)      &
             + VBWGT1(I,J,3)*U21(I,J,K,3)      &
             + VBWGT1(I,J,4)*U21(I,J,K,4)
            V3(I,J,K) =                       &
               VBWGT1(I,J,1)*V21(I,J,K,1)      &
             + VBWGT1(I,J,2)*V21(I,J,K,2)      &
             + VBWGT1(I,J,3)*V21(I,J,K,3)      &
             + VBWGT1(I,J,4)*V21(I,J,K,4)
         END DO
 115     CONTINUE
       ENDDO
       ENDDO

        print*,'JX,JY,NZ=',JX,JY,NZ
        write(*,*)'interpolation after d02='
      do k=1,NZ
        write(*,32)K,T3(9,9,K),            &
          Q3(9,9,K),U3(9,9,K),V3(9,9,K),Z3(9,9,K),P3(9,9,K)
      end do

      IF(IFLAG_NEST.EQ.2)THEN

! interpolate from d03

       PMV5=PMID5

       DO J=2,NY5-1
       DO I=2,NX5-1
       DO K=1,IZ
         PMV5(I,J,K)=0.25*(PMID5(I,J,K)+PMID5(I+1,J,K)+            &
                           PMID5(I,J+1,K)+PMID5(I+1,J+1,K))
       END DO
       END DO
       END DO

       DO J=1,JY
       DO I=1,JX
         IF(IIV2(I,J,1).LE.1.or.IIV2(I,J,1).GE.NX5_1)GO TO 885
         IF(JJV2(I,J,1).LE.1.or.JJV2(I,J,1).GE.NY5_1)GO TO 885
         DO N1=1,4
           IV1(N1)=IIV2(I,J,N1)
           JV1(N1)=JJV2(I,J,N1)
         END DO
         DO K=1,KZ
           DO N1=1,4
             IF(PMV3(I,J,K).GT.PMV5(IV1(N1),JV1(N1),1))THEN
               DP1=PMV3(I,J,K)-PMV5(IV1(N1),JV1(N1),1)
               U21(I,J,K,N1)=U5(IV1(N1),JV1(N1),1)
               V21(I,J,K,N1)=V5(IV1(N1),JV1(N1),1)
             ELSE IF(PMV3(I,J,K).LE.PMV5(IV1(N1),JV1(N1),NZ))THEN
               U21(I,J,K,N1)=U5(IV1(N1),JV1(N1),NZ)
               V21(I,J,K,N1)=V5(IV1(N1),JV1(N1),NZ)
             ELSE
               DO K1=2,KZ
                 PDIF1=ALOG(PMV3(I,J,K))-ALOG(PMV5(IV1(N1),JV1(N1),K1))
                 IF(PDIF1.GE.0.)THEN
                   FACT1=PDIF1/(ALOG(PMV5(IV1(N1),JV1(N1),K1-1))-             &
                                ALOG(PMV5(IV1(N1),JV1(N1),K1)))
                   U21(I,J,K,N1)=U5(IV1(N1),JV1(N1),K1)*(1.-FACT1)     &
                                +U5(IV1(N1),JV1(N1),K1-1)*FACT1
                   V21(I,J,K,N1)=V5(IV1(N1),JV1(N1),K1)*(1.-FACT1)     &
                                +V5(IV1(N1),JV1(N1),K1-1)*FACT1
                   GO TO 665
                 END IF
               END DO
             END IF
 665          CONTINUE
           END DO
         END DO
 885      CONTINUE
       END DO
       END DO
!23456789012345678901234567890123456789012345678901234567890123456789012

       DO J=1,JY
       DO I=1,JX
         IF(IIV2(I,J,1).LE.1.or.IIV2(I,J,1).GE.NX5_1)GO TO 215
         IF(JJV2(I,J,1).LE.1.or.JJV2(I,J,1).GE.NY5_1)GO TO 215
         DO K=1,KZ
            U3(I,J,K) =                       &
               VBWGT2(I,J,1)*U21(I,J,K,1)      &
             + VBWGT2(I,J,2)*U21(I,J,K,2)      &
             + VBWGT2(I,J,3)*U21(I,J,K,3)      &
             + VBWGT2(I,J,4)*U21(I,J,K,4)
            V3(I,J,K) =                       &
               VBWGT2(I,J,1)*V21(I,J,K,1)      &
             + VBWGT2(I,J,2)*V21(I,J,K,2)      &
             + VBWGT2(I,J,3)*V21(I,J,K,3)      &
             + VBWGT2(I,J,4)*V21(I,J,K,4)
         END DO
 215     CONTINUE
       ENDDO
       ENDDO

      END IF     ! IFLAG_NEST

        print*,'test74'

! save 4x data

      IUNIT=50+ITIM

      WRITE(IUNIT) JX,JY,NZ
!      WRITE(IUNIT) DLMD3,DPHD3,CLON0,CLAT0
      WRITE(IUNIT) DLMD3,DPHD3,CLON2,CLAT2
      WRITE(IUNIT) PT3,PDTOP3,WBD3,SBD3
      WRITE(IUNIT) T3
      WRITE(IUNIT) Q3
      WRITE(IUNIT) U3
      WRITE(IUNIT) V3
      WRITE(IUNIT) Z3
      WRITE(IUNIT) HLON3,HLAT3,VLON3,VLAT3             ! in degree
      WRITE(IUNIT) P3
      WRITE(IUNIT) PD3
      WRITE(IUNIT) DSG1
      WRITE(IUNIT) DSG2

      CLOSE(IUNIT)

! save A13 

      IF(IFLAG_NEST.EQ.2)THEN

      IF(IBGS.EQ.1)THEN
         DO J=1,NY5
         DO I=1,NX5
           J1=J+J_ST14-1
           I1=I+I_ST14-1
           A13(I1,J1)=A105(I,J)
           B13(I1,J1)=B105(I,J)
           C13(I1,J1)=C105(I,J)
         END DO
         END DO
       END IF

         IUNIT=60+ITIM

         WRITE(IUNIT) JX,JY
         WRITE(IUNIT) HLON3,HLAT3,VLON3,VLAT3
         WRITE(IUNIT) A13
         WRITE(IUNIT) B13
         WRITE(IUNIT) C13

         CLOSE(IUNIT) 

        END IF    ! IFLAG_NEST

        print*,'JX,JY,NZ=',JX,JY,NZ
        write(*,*)'inside merge K,T1,Q1,U1,V1,Z1,P1='
      do k=1,NZ
        write(*,32)K,T3(9,9,K),            &
          Q3(9,9,K),U3(9,9,K),V3(9,9,K),Z3(9,9,K),P3(9,9,K)
      end do
 32   format(I3,6F12.2)

      print *,'fort.61','JX=',JX,'JY=',JY

!      WRITE(61)((SLP3(I,J),I=1,JX),J=1,JY)
!      DO K=1,NZ+1
!        WRITE(61)((Z3(I,J,K),I=1,JX),J=1,JY)
!      END DO
!      DO K=1,NZ+1
!        WRITE(61)((P3(I,J,K),I=1,JX),J=1,JY)
!      END DO
!      DO K=1,NZ
!        WRITE(61)((T3(I,J,K),I=1,JX),J=1,JY)
!      END DO
!      DO K=1,NZ
!        WRITE(61)((Q3(I,J,K),I=1,JX),J=1,JY)
!      END DO
!      DO K=1,NZ
!        WRITE(61)((U3(I,J,K),I=1,JX),J=1,JY)
!      END DO
!      DO K=1,NZ
!        WRITE(61)((V3(I,J,K),I=1,JX),J=1,JY)
!      END DO


       END

!

