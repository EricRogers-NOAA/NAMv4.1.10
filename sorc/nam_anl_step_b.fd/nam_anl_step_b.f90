!******************************************************************************
! 
! PROGRAM HWRF RELOCATE 
! 
! SUBPROGRAM    
! PRGRMMR 
! ABSTRACT
! 
!______________________________________________________________________________ 

! DECLARE VARIABLES

      INTEGER I,J,K,NX,NY,NZ,IFLAG

!     PARAMETER (NX=215,NY=431,NZ=42,NST=5)
      PARAMETER (NST=5)
      PARAMETER (GAMMA=6.5E-3,G=9.8,Rd=287.05,D608=0.608)
      PARAMETER (Cp=1004.)

!     PARAMETER (KMX=2*NZ+1)
     
! Variables on 4x hybrid coordinate (ENV)

      REAL(4) DLMD,DPHD,PT,PDTOP
      REAL(4) WBD,SBD,CENTRAL_LON,CENTRAL_LAT

      REAL(4), ALLOCATABLE :: T1(:,:,:),Q1(:,:,:)
      REAL(4), ALLOCATABLE :: U1(:,:,:),V1(:,:,:)
      REAL(4), ALLOCATABLE :: Z1(:,:,:),P1(:,:,:)
      REAL(4), ALLOCATABLE :: GLON(:,:),GLAT(:,:)
      REAL(4), ALLOCATABLE :: PD1(:,:),ETA1(:),ETA2(:)

      REAL(4), ALLOCATABLE :: USCM(:,:),VSCM(:,:) ! Env. wind at new grids

      REAL(4), ALLOCATABLE :: HLON(:,:),HLAT(:,:)
      REAL(4), ALLOCATABLE :: VLON(:,:),VLAT(:,:)

! variables for hurricane component
 
      REAL(4), ALLOCATABLE :: HLON2(:,:),HLAT2(:,:)
      REAL(4), ALLOCATABLE :: VLON2(:,:),VLAT2(:,:)

      REAL(4), ALLOCATABLE :: SLPE(:,:),SLP_1(:,:),TENV(:,:,:)
!     REAL(4), ALLOCATABLE :: QENV(:,:,:)
      REAL(4), ALLOCATABLE :: T_1(:,:,:),Q_1(:,:,:)
      REAL(4), ALLOCATABLE :: U_1(:,:,:),V_1(:,:,:)

      REAL(4), ALLOCATABLE :: U_850(:,:),V_850(:,:)
!     REAL(4), ALLOCATABLE :: SLP_2(:,:)

! working array

      REAL(4), ALLOCATABLE :: SLP1(:,:),ANG2(:,:)
      REAL(4), ALLOCATABLE :: PMID1(:,:,:),ZMID1(:,:,:)
      REAL(4), ALLOCATABLE :: ZS1(:,:),TS1(:,:),QS1(:,:)
 
      REAL(4), ALLOCATABLE :: U_S(:,:),U_A(:,:)
      REAL(4), ALLOCATABLE :: V_S(:,:),V_A(:,:)
       
      REAL(4), ALLOCATABLE :: HLON3(:,:),HLAT3(:,:)
      REAL(4), ALLOCATABLE :: VLON3(:,:),VLAT3(:,:)

      REAL(4) CLON1,CLAT1

      REAL(4), ALLOCATABLE :: USC_1(:,:),VSC_1(:,:) ! hurr comp wind at level 1
      REAL(4), ALLOCATABLE :: USC1(:,:),VSC1(:,:) ! Hurricane wind at new grids
      REAL(4), ALLOCATABLE :: USC2(:,:),VSC2(:,:) ! Hurricane wind at new grids
      REAL(4), ALLOCATABLE :: SLPV(:,:)
 
      REAL(4), ALLOCATABLE :: HLON1(:,:),HLAT1(:,:)
      REAL(4), ALLOCATABLE :: VLON1(:,:),VLAT1(:,:)
      REAL(4), ALLOCATABLE :: T21(:,:,:,:),Q21(:,:,:,:)
      REAL(4), ALLOCATABLE :: U21(:,:,:,:),V21(:,:,:,:)
      REAL(4), ALLOCATABLE :: SLP21(:,:,:)
      REAL(4), ALLOCATABLE :: PMV1(:,:,:)

      REAL(4), ALLOCATABLE :: A101(:,:),B101(:,:),C101(:,:)

      REAL(4), ALLOCATABLE :: U_2(:,:,:),V_2(:,:,:)
      REAL(4), ALLOCATABLE :: U_3(:,:,:),V_3(:,:,:)
      REAL(4), ALLOCATABLE :: T_4(:,:,:),Q_4(:,:,:)

      REAL(8), ALLOCATABLE :: WRK1(:),WRK2(:),WRK3(:),WRK4(:)

      REAL(4), ALLOCATABLE :: PCST(:),HP(:,:,:)

      REAL(4), ALLOCATABLE :: PW(:)

      REAL(4), ALLOCATABLE ::    HBWGT1(:,:,:),VBWGT1(:,:,:)
      integer(4), ALLOCATABLE :: IIH1(:,:,:),JJH1(:,:,:)
      integer(4), ALLOCATABLE :: IIV1(:,:,:),JJV1(:,:,:)

      REAL(4), ALLOCATABLE :: T_X(:,:,:),Q_X(:,:,:),SLP_X(:,:)

      REAL(4), ALLOCATABLE :: dist(:,:)

      integer(4) IH1(4),JH1(4),IV1(4),JV1(4)

      CHARACTER ST_NAME(NST)*3,SN*1,EW*1,DEPTH*1
      REAL(4) CLON_NEW,CLAT_NEW,CLON_NHC,CLAT_NHC
      REAL(4) CLON_NEW1,CLAT_NEW1

      DIMENSION TWM(101),RWM(101),TH1(200),RP(200) 

      REAL(4) zmax 
     
      integer Ir_v4(4),   IR34_mod(4)
      REAL    R34_obs(4), R34_mod(4)

      CHARACTER PART1*2,basin*2,NUM*2
      LOGICAL TSTH

      COEF1=Rd/Cp
      COEF3=Rd*GAMMA/G
      COEF2=1./COEF3

      GRD=G/Rd

      pi=4.*atan(1.)
      pi_deg=180./pi     !* rad -> deg 
      pi180=pi/180.      !* deg -> rad 

      DST1=6.371E6*pi180 !* deg ->  m 
      arad=DST1          !* deg ->  m 

      eps1 = 1.E-1 
      eps2 = 1.E-2 
      eps3 = 1.E-3 
      eps4 = 1.E-4 
      eps5 = 1.E-5 
      eps6 = 1.E-6 
      v34kt= 34./1.944        !* m/s 
      v50kt= 50./1.944        !* m/s 
      v64kt= 64./1.944        !* m/s 
      rad2deg= 180./pi        !* rad -> deg 
      deg2rad= pi/180.        !* deg -> rad 
      deg2km= deg2rad*6.371E3 !* deg -> km 
      deg2m = deg2rad*6.371E6 !* deg ->  m 

      READ(5,*)ITIM,IGFS_FLAG,IFLAG_NEST
               
      print*,'ITIM,IGFS_FLAG,IFLAG_NEST=',ITIM,IGFS_FLAG,IFLAG_NEST
                
! Read TC vitals ... 
                               
      READ(11,11) id_storm,ICLAT,SN,ICLON,EW,Ipsfc,Ipcls,  &
                  Irmax,ivobs,Ir_vobs,(Ir_v4(I),I=1,4),DEPTH

 11   format(5x,I2,26x,I3,A1,I5,A1,9x,I4,1x,I4,1x,I4,I3,I4,4I5,1x,A1)

      CLAT_NHC=ICLAT*0.1
      CLON_NHC=ICLON*0.1

      IF(SN.eq.'S')CLAT_NHC=-CLAT_NHC
      IF(EW.eq.'W')CLON_NHC=-CLON_NHC
!     IF(EW.eq.'E')CLON_NHC=CLON_NHC-360.
      IF(CLON_NHC.GT.60.)CLON_NHC=CLON_NHC-360.

      psfc_obs = Ipsfc*100. !* pmin (Pa) 
      psfc_cls = Ipcls*100. !* pout (Pa)
      delt_p=psfc_obs-psfc_cls 
      PRMAX = Irmax*1.      !* ROCI (km) 

      vobs  = ivobs*1.0     !* Vmax (m/s)
      vobs_o= vobs 
      VRmax = Ir_vobs*1.    !* RMW  (km) 
      vobs  = max(vobs,1.) 
      vrmax = max(vrmax,19.) 

      R34obs = 0.
      R34obsm= 0. 
      acount = 0.
      R34_obs= 0. 
      DO i = 1, 4 
      if ( Ir_v4(i) > 0 ) then
!         if(Ir_v4(i) < 0 ) Ir_v4(i) = 0
         R34_obs(i) = Ir_v4(i) 
         R34obs = R34obs + Ir_v4(i) 
         acount = acount + 1.
         if(R34obsm.lt.R34_obs(i)) R34obsm = R34_obs(i) 
      endif 
      ENDDO 
      IF ( acount > 0.5 ) R34obs = R34obs/acount     !* avg R34 [km] 

      PRINT*, 'Obsereved Vmax:   vobs [m/s], DEPTH =', vobs, DEPTH
      PRINT*, 'Obsereved RMW:  Ir_vobs, VRmax [km] =', Ir_vobs, VRmax
      PRINT*, 'Obsereved R34: NE,SE,SW,NW,AVG [km] =', Ir_v4, R34obs,R34obsm 
      PRINT*, 'Obsereved ROCI:   Irmax, PRmax [km] =', Irmax, PRmax
      PRINT*, 'Obsereved Pressure: pmin, pout [Pa] =', psfc_obs, psfc_cls 

      REWIND 11

! READ 4x area env. data HWRF 

      IUNIT=20+ITIM

      READ(IUNIT) NX,NY,NZ

      print*,'NX,NY,NZ=',NX,NY,NZ

      NZ1=NZ+1

      KMX=2*NZ+1

      ALLOCATE ( T1(NX,NY,NZ),Q1(NX,NY,NZ) )
      ALLOCATE ( U1(NX,NY,NZ),V1(NX,NY,NZ) )
      ALLOCATE ( Z1(NX,NY,NZ+1),P1(NX,NY,NZ+1) )
      ALLOCATE ( GLON(NX,NY),GLAT(NX,NY) )
      ALLOCATE ( PD1(NX,NY),ETA1(NZ),ETA2(NZ) )
      ALLOCATE ( USCM(NX,NY),VSCM(NX,NY) )          ! Env. wind at new grids

      ALLOCATE ( HLON(NX,NY),HLAT(NX,NY) )
      ALLOCATE ( VLON(NX,NY),VLAT(NX,NY) )

      ALLOCATE ( dist(NX,NY) )

      READ(IUNIT) DLMD,DPHD,CENTRAL_LON,CENTRAL_LAT ! Domain res & center (deg)
      READ(IUNIT) PT,PDTOP,WBD,SBD
      READ(IUNIT) T1
      READ(IUNIT) Q1
      READ(IUNIT) U1
      READ(IUNIT) V1
      READ(IUNIT) Z1
!     READ(IUNIT) GLON,GLAT
      READ(IUNIT) HLON,HLAT,VLON,VLAT
      READ(IUNIT) P1
      READ(IUNIT) PD1
      READ(IUNIT) ETA1
      READ(IUNIT) ETA2
      READ(IUNIT) USCM
      READ(IUNIT) VSCM
 
      CLOSE(IUNIT)

      DO J=1,NY
      DO I=1,NX
	IF(HLON(I,J).GT.60.)HLON(I,J)=HLON(I,J)-360.
	IF(VLON(I,J).GT.60.)VLON(I,J)=VLON(I,J)-360.
        GLON(I,J)=HLON(I,J)
        GLAT(I,J)=HLAT(I,J)
      END DO
      END DO

      ALLOCATE ( SLP1(NX,NY),ANG2(NX,NY) )
      ALLOCATE ( PMID1(NX,NY,NZ),ZMID1(NX,NY,NZ) )
      ALLOCATE ( ZS1(NX,NY),TS1(NX,NY),QS1(NX,NY) )

      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
         PMID1(I,J,K)=EXP((ALOG(P1(I,J,K))+ALOG(P1(I,J,K+1)))*0.5)
         ZMID1(I,J,K)=0.5*(Z1(I,J,K)+Z1(I,J,K+1))
      ENDDO
      ENDDO
      ENDDO

! Compute sea level pressure 

      DO J=1,NY
      DO I=1,NX
         ZSF1 = ZMID1(I,J,1)
         PSF1 = PMID1(I,J,1)
         TV1 = T1(I,J,1)*(1.+D608*Q1(I,J,1))
         A = (GAMMA * ZSF1) / TV1
         SLP1(I,J) = PSF1*(1+A)**COEF2
      ENDDO
      ENDDO

      DO J=1,NY
      DO I=1,NX
         ZS1(I,J) = Z1(I,J,1)
         TS1(I,J) = T1(I,J,1)+GAMMA*(Z1(I,J,2)-Z1(I,J,1))*0.5
         QS1(I,J) = Q1(I,J,1)
      ENDDO
      ENDDO

! compute 10m wind
 
      IF(IFLAG_NEST.EQ.2)THEN

      IUNIT=40+ITIM
 
      READ(IUNIT) JX,JY
 
      ALLOCATE ( A101(JX,JY),B101(JX,JY),C101(JX,JY) )

      READ(IUNIT)
      READ(IUNIT) A101
      READ(IUNIT) B101
      READ(IUNIT) C101

      CLOSE(IUNIT)
 
      ELSE
        JX=NX
        JY=NY
        ALLOCATE ( A101(JX,JY),B101(JX,JY),C101(JX,JY) )
        A101=1.0
        B101=0.002
        C101=0.88
      END IF      

      PRINT*,'JX,JY,NX,NY=',JX,JY,NX,NY
 
! finsih compute 10m wind

!      WBD=-(NX-1)*DLMD
!      SBD=-((NY-1)/2)*DPHD

      write(*,*)'DLMD,DPHD,PT,PDTOP=',DLMD,DPHD,PT,PDTOP
      write(*,*)'WBD,SBD,CENTRAL_LON,CENTRAL_LAT=',    &
                 WBD,SBD,CENTRAL_LON,CENTRAL_LAT
      do k=1,nz
        write(*,*)'K,ETA1,ETA2=',K,ETA1(k),ETA2(k)
      end do

      print*,'CLON,CLAT=',GLON(1+(NX-1)/2,1+(NY-1)/2),   &
                          GLAT(1+(NX-1)/2,1+(NY-1)/2)
      print*,'SLON,SLAT=',GLON(1,1),           &
                          GLAT(1,1)

!      CLON1=HLON(1+(NX-1)/2,1+(NY-1)/2)
!      CLAT1=HLAT(1+(NX-1)/2,1+(NY-1)/2)
      CLON1=CENTRAL_LON
      CLAT1=CENTRAL_LAT

      ALLOCATE ( HLON3(NX,NY),HLAT3(NX,NY) )
      ALLOCATE ( VLON3(NX,NY),VLAT3(NX,NY) )

! LON & LAT at U,V

      CALL EARTH_LATLON_BGRID( HLAT3,HLON3,VLAT3,VLON3,  & ! lat, lon at H, V points
                               DLMD,DPHD,WBD,SBD,  &  ! input res, west & south bdry
                               CLAT1,CLON1,      & ! central lat,lon, all in degrees
                               NX,NY )

      print*,'HLAT,HLON,VLAT,VLON=',                  &
              HLAT(1,1),HLON(1,1),VLAT(1,1),VLON(1,1)
      print*,'HLAT3,HLON3,VLAT3,VLON3=',                  &
              HLAT3(1,1),HLON3(1,1),VLAT3(1,1),VLON3(1,1)

!           *     *     *     *     *     *     *     *     *     * 
   
! Read hurricane perturbation 

      ALLOCATE ( HLON2(NX,NY),HLAT2(NX,NY) )
      ALLOCATE ( VLON2(NX,NY),VLAT2(NX,NY) )

      ALLOCATE ( PCST(KMX),HP(NX,NY,KMX) )
      ALLOCATE ( SLPE(NX,NY),SLP_1(NX,NY),TENV(NX,NY,KMX) )
!     ALLOCATE ( QENV(NX,NY,KMX) )
      ALLOCATE ( T_1(NX,NY,KMX),Q_1(NX,NY,KMX) )
      ALLOCATE ( U_1(NX,NY,KMX),V_1(NX,NY,KMX) )
      ALLOCATE ( U_850(NX,NY),V_850(NX,NY) )
!     ALLOCATE ( SLP_2(NX,NY) )

      ALLOCATE ( U_S(NX,NY),U_A(NX,NY) )
      ALLOCATE ( V_S(NX,NY),V_A(NX,NY) )

      ALLOCATE ( USC_1(NX,NY),VSC_1(NX,NY) ) ! hurr comp wind at level 1
      ALLOCATE ( USC1(NX,NY),VSC1(NX,NY) )   ! Hurricane wind at new grids
      ALLOCATE ( USC2(NX,NY),VSC2(NX,NY) )   ! Hurricane wind at new grids
      ALLOCATE ( SLPV(NX,NY) )

      ALLOCATE ( T21(NX,NY,KMX,4),Q21(NX,NY,KMX,4) )
      ALLOCATE ( U21(NX,NY,KMX,4),V21(NX,NY,KMX,4) )
      ALLOCATE ( SLP21(NX,NY,4) )
      ALLOCATE ( PMV1(NX,NY,NZ) )

      ALLOCATE ( T_4(NX,NY,KMX),Q_4(NX,NY,KMX) )

      ALLOCATE ( T_X(NX,NY,KMX),Q_X(NX,NY,KMX),SLP_X(NX,NY) )

      ALLOCATE ( WRK1(KMX),WRK2(KMX),WRK3(KMX),WRK4(KMX) )

      ALLOCATE ( PW(KMX) )
 
      IR1=120

      SLP_1=0.
      T_1=0.
      Q_1=0.
      U_1=0.
      V_1=0. 
                 
      T_X=0.
      Q_X=0.
      SLP_X=0.
 
      NCHT=71
      READ(NCHT)KSTM
      PRINT*,'test1',KSTM

      READ(NCHT)HLAT2,HLON2
      READ(NCHT)VLAT2,VLON2

      print*,'HLAT2,HLON2=',HLAT2(1,1),HLON2(1,1)
      print*,'VLAT2,VLON2=',VLAT2(1,1),VLON2(1,1)

      deltp=1.e20

      READ(NCHT)PCST
      DO K=1,KMX
         PRINT*,'K,PCST=',K,PCST(K)
         deltp1=abs(PCST(K)-85000.)
         IF (deltp1.LT.deltp) THEN
            deltp=deltp1
            k850=k
         END IF
      END DO

      k850=1       ! use the surface wind
 
      READ(NCHT)HP

      KST=1
      READ(NCHT) ST_NAME(KST)
      PRINT*,'ST_NAME=',ST_NAME(KST)
      READ(NCHT) CLON_NEW,CLAT_NEW

      IF (CLON_NEW.gt.60.) CLON_NEW=CLON_NEW-360.
      DO J=1,NY
      DO I=1,NX
         IF(HLON2(I,J).GT.60.) HLON2(I,J)=HLON2(I,J)-360.
         IF(VLON2(I,J).GT.60.) VLON2(I,J)=VLON2(I,J)-360.
      END DO
      END DO

      PRINT*,CLON_NEW,CLAT_NEW

      READ(NCHT) zmax
      print*,'zmax=',zmax

      READ(NCHT)IWMIN1,IWMAX1,JWMIN1,JWMAX1
      PRINT*,IWMIN1,IWMAX1,JWMIN1,JWMAX1
      READ(NCHT)((SLPE(I,J),I=1,NX),J=1,NY)          ! SLP
      READ(NCHT)((SLP_1(I,J),I=IWMIN1,IWMAX1),J=JWMIN1,JWMAX1)    ! pert SLP
!!     WRITE(25)((SLP_1(I,J),I=1,NX),J=1,NY)
      PRINT*,'TEST1'
      DO K=1,KMX
         READ(NCHT)((TENV(I,J,K),I=1,NX),J=1,NY)
         READ(NCHT)((T_1(I,J,K),I=IWMIN1,IWMAX1),J=JWMIN1,JWMAX1)
!        WRITE(25) ((T_1(I,J,K),I=1,NX),J=1,NY)
      END DO
      DO K=1,KMX
         READ(NCHT)((U_1(I,J,K),I=IWMIN1,IWMAX1),J=JWMIN1,JWMAX1)
         READ(NCHT)((V_1(I,J,K),I=IWMIN1,IWMAX1),J=JWMIN1,JWMAX1)
!        WRITE(25)((U_1(I,J,K),I=1,NX),J=1,NY)
      END DO
!     DO K=1,KMX
!        WRITE(25)((V_1(I,J,K),I=1,NX),J=1,NY)
!     END DO

      PRINT*,'TEST31'
      DO K=1,KMX
!!!      READ(NCHT)((QENV(I,J,K),I=1,NX),J=1,NY)
         READ(NCHT)((Q_1(I,J,K),I=IWMIN1,IWMAX1),J=JWMIN1,JWMAX1)
!        WRITE(25) ((Q_1(I,J,K),I=1,NX),J=1,NY)
         DO J=JWMIN1,JWMAX1
         DO I=IWMIN1,IWMAX1
            Q_1(I,J,K)=Q_1(I,J,K)*1.E-3
         END DO
         END DO
         PRINT*,'K,TEST21',K
      END DO

      U_850=0.
      V_850=0.

      IF(IFLAG_NEST.EQ.2)THEN
     
      READ(NCHT)((U_850(I,J),I=IWMIN1,IWMAX1),J=JWMIN1,JWMAX1)
      READ(NCHT)((V_850(I,J),I=IWMIN1,IWMAX1),J=JWMIN1,JWMAX1)

      ELSE

        DO J=1,NY
        DO I=1,NX
           U_850(I,J)=U1(I,J,1)+U_1(I,J,1)
           V_850(I,J)=V1(I,J,1)+V_1(I,J,1)
        END DO
        END DO

      END IF     ! IFLAG_NEST

      CLOSE(NCHT)

!           *     *     *     *     *     *     *     *     *     * 
   
      aaa=1.0
      bbb=0.     ! for NAM relocation

      fact=1.0

      distm=1.E20
      do j=1,ny
      do i=1,nx !* Search for indices nearest to center of storm 
         distt=(HLON2(i,j)-CLON_NEW)**2+(HLAT2(i,j)-CLAT_NEW)**2
         if (distm.GT.distt) then
            distm=distt
            icst=i
            jcst=j
         end if
      end do
      end do

      print*,'qliu test icst,jcst=',icst,jcst
      print*,'HLON2,HLAT2(NX/2,NY/2)=',HLON2(NX/2,NY/2),HLAT2(NX/2,NY/2)
      print*,'CLON_NEW,CLAT_NEW=',CLON_NEW,CLAT_NEW

      PRINT*,'IWMIN1,IWMAX1,JWMIN1,JWMAX1 before=',IWMIN1,IWMAX1,JWMIN1,JWMAX1
 
      ID_INDX=MAX(icst-IWMIN1,IWMAX1-icst) ! delete "*fact", OK 
      JD_INDX=MAX(jcst-JWMIN1,JWMAX1-jcst) ! delete "*fact", OK 
      IWMIN1=MAX(icst-ID_INDX,1)-8
      IWMAX1=MIN(icst+ID_INDX,NX)+8
      JWMIN1=MAX(jcst-JD_INDX,1)-8
      JWMAX1=MIN(jcst+JD_INDX,NY)+8

!      IWMIN1=1
!      IWMAX1=NX
!      JWMIN1=1
!      JWMAX1=NY 
 
      PRINT*,'IWMIN1,IWMAX1,JWMIN1,JWMAX1 after=',IWMIN1,IWMAX1,JWMIN1,JWMAX1
 
      DO J=1,NY
      DO I=1,NX
        USC_1(I,J)=U_1(I,J,1)
        VSC_1(I,J)=V_1(I,J,1)
      END DO
      END DO

!           *     *     *     *     *     *     *     *     *     * 
   
! scale the coordinate

      iparam=2
 
      IF ( iparam == 1 ) THEN !* --------------------------------------- 

      DO J=1,NY
      DO I=1,NX
         HLAT2(I,J)=(HLAT2(I,J)-CLAT_NEW)*fact+CLAT_NHC
         HLON2(I,J)=(HLON2(I,J)-CLON_NEW)*fact+CLON_NHC
         VLAT2(I,J)=(VLAT2(I,J)-CLAT_NEW)*fact+CLAT_NHC
         VLON2(I,J)=(VLON2(I,J)-CLON_NEW)*fact+CLON_NHC
      END DO
      END DO

      ELSEIF ( iparam == 2 ) THEN !* for inverse stretch --------------- 

      DO J=1,NY
      DO I=1,NX
         HLAT2(I,J)=(HLAT2(I,J)-CLAT_NEW)+CLAT_NHC
         HLON2(I,J)=(HLON2(I,J)-CLON_NEW)+CLON_NHC
         VLAT2(I,J)=(VLAT2(I,J)-CLAT_NEW)+CLAT_NHC
         VLON2(I,J)=(VLON2(I,J)-CLON_NEW)+CLON_NHC
      END DO
      END DO

      ENDIF !* --------------------------------------------------------- 

!           *     *     *     *     *     *     *     *     *     * 
      KNHC=1+(NX-1)/2
      MNHC=1+(NY-1)/2
      IC1=KNHC+1
      JC1=MNHC+1
      DDX=(GLON(IC1,MNHC)-GLON(KNHC,MNHC))
      DDY=(GLAT(KNHC,JC1)-GLAT(KNHC,MNHC))
      MDX=((CLON_NHC-CLON_NEW)/DDX)
      MDY=((CLAT_NHC-CLAT_NEW)/DDY)
 
      PRINT*,'DDX,DDY,MDX,MDY=',DDX,DDY,MDX,MDY
      PRINT*,'IWMIN1,IWMAX1,JWMIN1,JWMAX1=',        &
              IWMIN1,IWMAX1,JWMIN1,JWMAX1

      IF ( iparam == 1 ) THEN !* - - - - - - - - - - - - - - - - - - - - 

         DLMD2=DLMD*fact 
         DPHD2=DPHD*fact 
         WBD2=WBD*fact
         SBD2=SBD*fact

      ELSEIF ( iparam == 2 ) THEN !* for inverse stretch - - - - - - - - 

         DLMD2=DLMD
         DPHD2=DPHD
         WBD2=WBD
         SBD2=SBD

      ENDIF !* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

!      WBD2=-(NX-1)*DLMD2 
!      SBD2=-((NY-1)/2)*DPHD2 

!      CENTRAL_LON2=HLON2(1+(NX-1)/2,1+(NY-1)/2) 
!      CENTRAL_LAT2=HLAT2(1+(NX-1)/2,1+(NY-1)/2) 
      CENTRAL_LON2=CENTRAL_LON
      CENTRAL_LAT2=CENTRAL_LAT

      PRINT*,'DLMD2,DPHD2,WBD2,SBD2,CENTRAL_LAT2,CENTRAL_LON2=', &
              DLMD2,DPHD2,WBD2,SBD2,CENTRAL_LAT2,CENTRAL_LON2

      NXT=IWMAX1-IWMIN1+1
      NYT=JWMAX1-JWMIN1+1
      NXT1=min(NXT+1,NX)
      NYT1=min(NYT+1,NY)

      ALLOCATE ( HLON1(NXT,NYT),HLAT1(NXT,NYT) )
      ALLOCATE ( VLON1(NXT,NYT),VLAT1(NXT,NYT) )
      ALLOCATE ( IIH1(NXT,NYT,4), JJH1(NXT,NYT,4)  )
      ALLOCATE ( IIV1(NXT,NYT,4), JJV1(NXT,NYT,4)  )
      ALLOCATE ( HBWGT1(NXT,NYT,4), VBWGT1(NXT,NYT,4) )

      IF(FLAG_NEST.NE.0)THEN

      IF ( iparam == 1 ) THEN !* --------------------------------------- 
 
      DO J=1,NYT
      DO I=1,NXT
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
         HLAT1(I,J)=HLAT3(I1,J1) 
         HLON1(I,J)=HLON3(I1,J1) 
         VLAT1(I,J)=VLAT3(I1,J1) 
         VLON1(I,J)=VLON3(I1,J1) 
      END DO
      END DO
 
      PRINT*, 'Old storm-size correction with 1 parameter.' 
 
      ELSEIF ( iparam == 2 ) THEN !* inverse stretch ------------------- 
 
!      if ( abs(bbb) < eps5 ) then !* - - - - - - - - - - - - - - - - 
      if ( abs(bbb) < eps4 ) then !* - - - - - - - - - - - - - - - - 
 
      DO J=1,NYT
      DO I=1,NXT
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
         HLON1(I,J) = (HLON3(I1,J1)-CLON_NEW)/aaa + CLON_NEW
         HLAT1(I,J) = (HLAT3(I1,J1)-CLAT_NEW)/aaa + CLAT_NEW  
         VLON1(I,J) = (VLON3(I1,J1)-CLON_NEW)/aaa + CLON_NEW  
         VLAT1(I,J) = (VLAT3(I1,J1)-CLAT_NEW)/aaa + CLAT_NEW  
      ENDDO
      ENDDO
 
      PRINT*, '2-parameter storm-size correction reduced to 1 parameter.' 
 
      else  !* - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 
      PRINT*, 'Storm-size correction with 2 parameters.' 
 
      DO J=1,NYT
      DO I=1,NXT
 
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
 
         xxx = HLON3(I1,J1)-CLON_NEW 
         yyy = HLAT3(I1,J1)-CLAT_NEW 
         zzz = sqrt(xxx**2 + yyy**2) 
         IF ( zzz < eps6 ) THEN 
            HLON1(I,J) = CLON_NEW
            HLAT1(I,J) = CLAT_NEW
         ELSE 
            ddd = aaa**2 + 2.*bbb*zzz 
            IF ( ddd < 0. ) THEN 
               HLON1(I,J) = xxx/aaa + CLON_NEW 
               HLAT1(I,J) = yyy/aaa + CLAT_NEW 
!               PRINT*,'Storm center is too far off analysis domain center !' 
!               PRINT*,'Constant inverse stretch for outer storm area.' 
            ELSE !* Variable inverse stretch 
               rrr = ( sqrt(ddd) - aaa )/bbb 
               HLON1(I,J) = rrr*xxx/zzz + CLON_NEW 
               HLAT1(I,J) = rrr*yyy/zzz + CLAT_NEW 
            ENDIF 
         ENDIF 
 
         xxx = VLON3(I1,J1)-CLON_NEW 
         yyy = VLAT3(I1,J1)-CLAT_NEW 
         zzz = sqrt(xxx**2 + yyy**2) 
         IF ( zzz < eps6 ) THEN  
            VLON1(I,J) = CLON_NEW
            VLAT1(I,J) = CLAT_NEW
         ELSE 
            ddd = aaa**2 + 2.*bbb*zzz 
            IF ( ddd < 0. ) THEN 
               VLON1(I,J) = xxx/aaa + CLON_NEW 
               VLAT1(I,J) = yyy/aaa + CLAT_NEW 
!               PRINT*,'Storm center is too far off analysis domain center !' 
!               PRINT*,'Constant inverse stretch for outer storm area.' 
            ELSE !* Variable inverse stretch 
               rrr = ( sqrt(ddd) - aaa )/bbb 
               VLON1(I,J) = rrr*xxx/zzz + CLON_NEW 
               VLAT1(I,J) = rrr*yyy/zzz + CLAT_NEW 
            ENDIF 
         ENDIF 
 
      ENDDO
      ENDDO
 
      endif !* - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 
        DO J=1,NYT
        DO I=1,NXT
           HLAT1(I,J)=HLAT1(I,J)+CLAT_NEW-CLAT_NHC
           HLON1(I,J)=HLON1(I,J)+CLON_NEW-CLON_NHC
           VLAT1(I,J)=VLAT1(I,J)+CLAT_NEW-CLAT_NHC
           VLON1(I,J)=VLON1(I,J)+CLON_NEW-CLON_NHC
        END DO
        END DO

      ENDIF !* --------------------------------------------------------- 

      ELSE      ! FLAG_NEST=0

      DO J=1,NYT
      DO I=1,NXT
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
         HLAT1(I,J)=HLAT3(I1,J1)
         HLON1(I,J)=HLON3(I1,J1)
         VLAT1(I,J)=VLAT3(I1,J1)
         VLON1(I,J)=VLON3(I1,J1)
      END DO
      END DO


        DO J=1,NYT
        DO I=1,NXT
           HLAT1(I,J)=HLAT1(I,J)+CLAT_NEW-CLAT_NHC
           HLON1(I,J)=HLON1(I,J)+CLON_NEW-CLON_NHC
           VLAT1(I,J)=VLAT1(I,J)+CLAT_NEW-CLAT_NHC
           VLON1(I,J)=VLON1(I,J)+CLON_NEW-CLON_NHC
        END DO
        END DO
     
     END IF

     print*,'CLAT_NEW,CLAT_NHC=',CLAT_NEW,CLAT_NHC
     print*,'CLON_NEW,CLON_NHC=',CLON_NEW,CLON_NHC

!      CLAT_NEW1=CLAT_NEW
!      CLON_NEW1=CLON_NEW
      CLAT_NEW=CLAT_NHC
      CLON_NEW=CLON_NHC
 
      print*,'domain 3=',HLON3(1,1),HLON3(NX,NY),HLAT3(1,1),HLAT3(NX,NY)
      print*,'domain 1=',HLON1(1,1),HLON1(NXT,NYT),HLAT1(1,1),HLAT1(NXT,NYT)

!*    Calculate interpolating indices and weights from modified vortex 
!*    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
!*    Input DLMD2,DPHD2,CENTRAL_LAT,CENTRAL_LON in degrees 
!*    Input HLAT1,HLON1,VLAT1,VLON1 in degrees 

      IIH1=999999
      JJH1=999999
      IIV1=999999
      JJV1=999999
      HBWGT1=0.
      VBWGT1=0.


      CALL G2T2H_BGRID( IIH1,JJH1,  & !* interp indices of source H points
           HBWGT1,                      & !* interp weights to target H points
           HLAT1,HLON1,                 & !* lat/lon of H points on target grid
           DLMD2,DPHD2,WBD2,SBD2,       & !* dlat/dlon, W/S bdry of source grid
           CENTRAL_LAT2,CENTRAL_LON2,   & !* central lat/lon of source  grid
           NX,NY,                       & !*  x/y dimensions of source  grid
           NXT,NYT )                      !* dimensions of target grid (nest)

      CALL G2T2V_BGRID( IIV1,JJV1,  & !* interp indices of source V points
           VBWGT1,                      & !* interp weights to target V points
           VLAT1,VLON1,                 & !* lat/lon of V points on target grid
           DLMD2,DPHD2,WBD2,SBD2,       & !* dlat/dlon, W/S bdry of source grid
           CENTRAL_LAT2,CENTRAL_LON2,   & !* central lat/lon of source  grid
           NX,NY,                       & !*  x/y dimensions of source  grid
           NXT,NYT )                      !* dimensions of target grid (nest)


      print*,'finish calculate index and coef'

      DEALLOCATE( HLAT1, HLON1, VLAT1, VLON1 ) 

!           *     *     *     *     *     *     *     *     *     * 
   
! Add hurricane component for surface wind (approximate for pert only)

      U21=0.
      V21=0.
      
      USC1=USC_1
      VSC1=VSC_1

      NX_1=NX-1
      NY_1=NY-1

      DO J=1,NYT
      DO I=1,NXT
         DO N1=1,4
            IV1(N1)=IIV1(I,J,N1)
            JV1(N1)=JJV1(I,J,N1)
         END DO
         DO N1=1,4
            U21(I,J,1,N1)=USC_1(IV1(N1),JV1(N1))
            V21(I,J,1,N1)=VSC_1(IV1(N1),JV1(N1))
         END DO
       ENDDO
       ENDDO

      DO J=1,NYT
      DO I=1,NXT
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
         USC1(I1,J1) =                                     &
                     + VBWGT1(I,J,1)*U21(I,J,1,1)          &
                     + VBWGT1(I,J,2)*U21(I,J,1,2)          &
                     + VBWGT1(I,J,3)*U21(I,J,1,3)          &
                     + VBWGT1(I,J,4)*U21(I,J,1,4)
         VSC1(I1,J1) =                                     &
                     + VBWGT1(I,J,1)*V21(I,J,1,1)          &
                     + VBWGT1(I,J,2)*V21(I,J,1,2)          &
                     + VBWGT1(I,J,3)*V21(I,J,1,3)          &
                     + VBWGT1(I,J,4)*V21(I,J,1,4)
 36      CONTINUE
      END DO
      END DO

      DO J=1,NY
      DO I=1,NX
         USCM(I,J)=U1(I,J,1)
         VSCM(I,J)=V1(I,J,1)
      END DO
      END DO

      vmax1=0.
      d_max=3.5
      IF(vobs.gt.30..and.CLAT_NHC.gt.30.)d_max=4.5

      dist=1.e20

      DO J=1,NYT
         J1=J+JWMIN1-1+MDY
         DO I=1,NXT
            I1=I+IWMIN1-1+MDX
            vmax2=USC1(I1,J1)**2+VSC1(I1,J1)**2
            R05=SQRT((VLAT(I1,J1)-CLAT_NEW)**2+((VLON(I1,J1)-CLON_NEW)*cost)**2)
            dist(I1,J1)=R05
            if (vmax2.gt.vmax1.and.R05.le.d_max) then
               vmax1=vmax2
               imax1=I1
               jmax1=j1
            end if
         END DO
      END DO
 
      print*,'I,J,vmax=',imax1,jmax1,sqrt(vmax1)
                               
      RMX_d=1.2*sqrt(((VLON(imax1,jmax1)-CLON_NEW)*cost)**2+           &
                     (VLAT(imax1,jmax1)-CLAT_NEW)**2)

      PRINT*,'RMX_d1=',RMX_d

      RMX_d=2.*RWMAX*fact !* OK 

      RMX_d=max(RMX_d,2.0)
!      RMX_d=max(RMX_d,3.0)
      IF(RMX_d.gt.3.5)RMX_d=3.5

      PRINT*,'RMX_d2=',RMX_d

      smax1=0.
      DO J=1,NYT
         J1=J+JWMIN1-1+MDY
         DO I=1,NXT
            I1=I+IWMIN1-1+MDX
!            R05=SQRT((VLAT(I1,J1)-CLAT_NEW)**2+((VLON(I1,J1)-CLON_NEW)*cost)**2)
            R05=dist(I1,J1)
            smax2=SQRT((USC1(I1,J1)+USCM(I1,J1))**2+(VSC1(I1,J1)+VSCM(I1,J1))**2)*C101(I1,J1)
            if (smax2.gt.smax1.and.R05.LT.RMX_d) then
               smax1=smax2
               i_max=I1
               j_max=J1
            end if
         END DO
      END DO

      vs_t=smax1
      vmax_s=smax1
      PRINT*, 'Observed initial wind [m/s] at 10m level: vobs =', vobs 
      PRINT*, 'Adjusted prev.6h wind [m/s] at 10m level: vs_t =', vs_t 

      DO J=1,NY
      DO I=1,NX
         USC2(I,J)=USC1(I,J)+USCM(I,J)
         VSC2(I,J)=VSC1(I,J)+VSCM(I,J)
      END DO
      END DO

      IFLAG=0
      PRINT*,'IFLAG0=',IFLAG

      RWMAX1=RWMAX*fact*111.12*cos(CLAT_NEW*pi180) 
      if (vobs.LE.24.) THEN
!        RWMAX1=VRmax
         TWMAX=TWMAX
      end if

      vs_t=1.0        ! for NAM relocation

   IF ( vs_t .LT. vobs) THEN !* --------------------------------------------> 

      IFLAG=1
      PRINT*,'IFLAG1=',IFLAG
!      OPEN(69,file='flag_file',form='formatted')
!         WRITE(69,*)IFLAG
!         write(69,*)K850
!         write(69,*)TWMAX,RWMAX1,fact,psfc_obs 
!      CLOSE(69)

      PRINT*, 'Recycled vortex will be enhanced by bogus vortex.' 

      beta=1.0
      ics=2

      IF(IFLAG_NEST.NE.0)THEN

      CALL CORT_MAT_1(IR1,NX,NY,NZ,KMX,                   &
		      T_X,Q_X,SLP_X,HLON,HLAT,VLON,VLAT,  &
		      CLON_NEW,CLAT_NEW,                  &
                      beta,fact,aaa,bbb,iparam,ics) 

      DO J=1,NY
      DO I=1,NX
         SLP_1(I,J)=SLP_1(I,J)*beta+SLP_X(I,J)
      END DO
      END DO

      pct_m=0.
      DO J=1,NY
      DO I=1,NX
	 IF (pct_m.GT.SLP_1(I,J)) THEN
	    pct_m=SLP_1(I,J)
	 END IF
      END DO
      END DO
!     pct_m=pct_m+50.

      ps_rat2=psfc_obs1/(pct_m-100.)              ! - 1mb 

      IF (ps_rat2.LT.1.0) THEN
     
         print*,'ps_rat2=',ps_rat2,psfc_obs1,pct_m

         DO J=1,NY
         DO I=1,NX
            SLP_1(I,J)=SLP_1(I,J)*ps_rat2
            DO K=1,KMX
	       TEK1=TENV(I,J,K)+T_1(I,J,K)
               T_1(I,J,K)=(T_1(I,J,K)*beta+T_X(I,J,K))*ps_rat2
	       TEK2=TENV(I,J,K)+T_1(I,J,K)
	       ESRR=exp(4302.645*(TEK2-TEK1)/((TEK2-29.66)*(TEK1-29.66)))
	       Q_1(I,J,K)=ESRR*Q_1(I,J,K)
            END DO
         END DO
         END DO

      ELSE

         DO J=1,NY
         DO I=1,NX
         DO K=1,KMX
	    TEK1=TENV(I,J,K)+T_1(I,J,K)
            T_1(I,J,K)=T_1(I,J,K)*beta+T_X(I,J,K)
!           Q_1(I,J,K)=Q_1(I,J,K)*beta+Q_X(I,J,K)
	    TEK2=TENV(I,J,K)+T_1(I,J,K)
	    ESRR=exp(4302.645*(TEK2-TEK1)/((TEK2-29.66)*(TEK1-29.66)))
	    Q_1(I,J,K)=ESRR*Q_1(I,J,K)
         END DO
         END DO
         END DO
   
      END IF

      END IF           !IF(IFLAG_NEST.NE.0)THEN


! Add hurricane component and write the guess field out
! the bogus program will add some of the storm increments

      SLP21=0.
      U21=0.
      V21=0.
      T21=0.
      Q21=0.

      DO J=1,NYT
      DO I=1,NXT
         DO N1=1,4
           IH1(N1)=IIH1(I,J,N1)
           JH1(N1)=JJH1(I,J,N1)
         END DO
         DO N1=1,4
            SLP21(I,J,N1)=SLP_1(IH1(N1),JH1(N1))
            DO K=1,KMX
               T21(I,J,K,N1)=T_1(IH1(N1),JH1(N1),K)
               Q21(I,J,K,N1)=Q_1(IH1(N1),JH1(N1),K)
            END DO
         END DO
      ENDDO
      ENDDO

      print*,'NXT,NYT,IWMIN1,JWMIN1,MDX,MDY=',NXT,NYT,IWMIN1,JWMIN1,MDX,MDY

      DO J=1,NYT
      DO I=1,NXT
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
         SLP1(I1,J1)=SLP1(I1,J1)+                      &
                     HBWGT1(I,J,1)*SLP21(I,J,1)        &
                   + HBWGT1(I,J,2)*SLP21(I,J,2)        &
                   + HBWGT1(I,J,3)*SLP21(I,J,3)        &
                   + HBWGT1(I,J,4)*SLP21(I,J,4)
         DO K=1,KMX
            WRK1(K) =                                  &
                HBWGT1(I,J,1)*T21(I,J,K,1)             &
              + HBWGT1(I,J,2)*T21(I,J,K,2)             &
              + HBWGT1(I,J,3)*T21(I,J,K,3)             &
              + HBWGT1(I,J,4)*T21(I,J,K,4)
            WRK2(K) =                                  &
                HBWGT1(I,J,1)*Q21(I,J,K,1)             &
              + HBWGT1(I,J,2)*Q21(I,J,K,2)             &
              + HBWGT1(I,J,3)*Q21(I,J,K,3)             &
              + HBWGT1(I,J,4)*Q21(I,J,K,4)
         END DO
         DO N=1,NZ
            TENV1 = T1(I1,J1,N)
            QENV1 = Q1(I1,J1,N)
            IF (PMID1(I1,J1,N).GE.PCST(1)) THEN        ! Below PCST(1)
               T1(I1,J1,N)=TENV1+WRK1(1)
               Q1(I1,J1,N)=QENV1+WRK2(1)
            ELSE IF(PMID1(I1,J1,N).LE.PCST(KMX))THEN
               T1(I1,J1,N)=TENV1+WRK1(KMX)
               Q1(I1,J1,N)=QENV1+WRK2(KMX)
            ELSE
               DO K=1,KMX-1
               if(PMID1(I1,J1,N).LE.PCST(K).and.PMID1(I1,J1,N).GT.PCST(K+1))then
                  W1=ALOG(1.*PCST(K+1))-ALOG(1.*PCST(K))
                  W=(ALOG(1.*PMID1(I1,J1,N))-ALOG(1.*PCST(K)))/W1
                  T1(I1,J1,N)=TENV1+WRK1(K)*(1.-W)+WRK1(K+1)*W
                  Q1(I1,J1,N)=QENV1+WRK2(K)*(1.-W)+WRK2(K+1)*W
                  GO TO 877
               end if
               END DO
            END IF
 877        CONTINUE
         END DO
 46      CONTINUE
      ENDDO
      ENDDO
           
! based on Ts, Zs, SLP1 ==> PS1  ==> P1
 
      DO J=1,NY
      DO I=1,NX
         ZSFC = ZS1(I,J)
         TSFC = TS1(I,J)*(1.+D608*QS1(I,J))
         A = (GAMMA * ZSFC) / TSFC
         P1(I,J,1) = SLP1(I,J)/(1+A)**COEF2
!         PD1(I,J)=P1(I,J,1)-PDTOP-PT
         PD1(I,J)=P1(I,J,1)-PT
      ENDDO
      ENDDO

! PD(I,J)=P1(I,J,1)-PDTOP-PT=PSFC(I,J)-PDTOP-PT
      DO J=1,NY
      DO I=1,NX
        P1(I,J,NZ1)=PT
        DO K=1,NZ
         P1(I,J,NZ1-K)=P1(I,J,NZ1-K+1)+PDTOP*ETA1(K)+PD1(I,J)*ETA2(K)     ! PD(I,J) changed
        ENDDO
      ENDDO
      ENDDO
 
      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
         PMID1(I,J,K)=EXP((ALOG(P1(I,J,K))+ALOG(P1(I,J,K+1)))*0.5)
      ENDDO
      ENDDO
      ENDDO
 
      DO j = 1,ny
      DO i = 1,nx
         Z1(I,J,1)=ZS1(I,J)
         DO L=2,nz+1
            Z1(I,J,L)=Z1(I,J,L-1)+T1(I,J,L-1)*             &
                     (Q1(I,J,L-1)*0.608+1.0)*287.04*       &
                     (ALOG(P1(I,J,L-1))-ALOG(P1(I,J,L)))/G
         ENDDO
      ENDDO
      ENDDO

      PMV1=PMID1
 
      DO J=2,NY-1
      DO K=1,NZ
      DO I=2,NX-1
         PMV1(I,J,K)=0.25*(PMID1(I,J,K)+PMID1(I+1,J,K)+  &
                           PMID1(I,J+1,K)+PMID1(I+1,J+1,K))
      END DO
      END DO
      ENDDO

      DO J=1,NYT
      DO I=1,NXT
         DO N1=1,4
            IV1(N1)=IIV1(I,J,N1)
            JV1(N1)=JJV1(I,J,N1)
         END DO
         DO N1=1,4
         DO K=1,KMX
            U21(I,J,K,N1)=U_1(IV1(N1),JV1(N1),K)
            V21(I,J,K,N1)=V_1(IV1(N1),JV1(N1),K)
         END DO
         END DO
      ENDDO
      ENDDO

      DO J=1,NYT
      DO I=1,NXT
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
         USCM(I1,J1) = USCM(I1,J1)                     &
                     + VBWGT1(I,J,1)*U21(I,J,1,1)      &
                     + VBWGT1(I,J,2)*U21(I,J,1,2)      &
                     + VBWGT1(I,J,3)*U21(I,J,1,3)      &
                     + VBWGT1(I,J,4)*U21(I,J,1,4)
         VSCM(I1,J1) = VSCM(I1,J1)                     &
                     + VBWGT1(I,J,1)*V21(I,J,1,1)      &
                     + VBWGT1(I,J,2)*V21(I,J,1,2)      &
                     + VBWGT1(I,J,3)*V21(I,J,1,3)      &
                     + VBWGT1(I,J,4)*V21(I,J,1,4)
 56     CONTINUE
      END DO
      END DO

      DO J=1,NYT
      DO I=1,NXT
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
         DO K=1,KMX
            WRK1(K) =                          &
               VBWGT1(I,J,1)*U21(I,J,K,1)      &
             + VBWGT1(I,J,2)*U21(I,J,K,2)      &
             + VBWGT1(I,J,3)*U21(I,J,K,3)      &
             + VBWGT1(I,J,4)*U21(I,J,K,4)
            WRK2(K) =                          &
               VBWGT1(I,J,1)*V21(I,J,K,1)      &
             + VBWGT1(I,J,2)*V21(I,J,K,2)      &
             + VBWGT1(I,J,3)*V21(I,J,K,3)      &
             + VBWGT1(I,J,4)*V21(I,J,K,4)
         END DO
         DO N=1,NZ
            IF (PMV1(I1,J1,N).GE.PCST(1)) THEN            ! Below PCST(1)
               U1(I1,J1,N)=U1(I1,J1,N)+WRK1(1)
               V1(I1,J1,N)=V1(I1,J1,N)+WRK2(1)
            ELSE IF(PMV1(I1,J1,N).LE.PCST(KMX)) THEN
               U1(I1,J1,N)=U1(I1,J1,N)+WRK1(KMX)
               V1(I1,J1,N)=V1(I1,J1,N)+WRK2(KMX)
            ELSE
               DO K=1,KMX-1
               if(PMV1(I1,J1,N).LE.PCST(K).and.PMV1(I1,J1,N).GT.PCST(K+1))then
                  W1=ALOG(1.*PCST(K+1))-ALOG(1.*PCST(K))
                  W=(ALOG(1.*PMV1(I1,J1,N))-ALOG(1.*PCST(K)))/W1
                  U1(I1,J1,N)=U1(I1,J1,N)+WRK1(K)*(1.-W)+WRK1(K+1)*W
                  V1(I1,J1,N)=V1(I1,J1,N)+WRK2(K)*(1.-W)+WRK2(K+1)*W
                  GO TO 881
               end if
               END DO
            END IF
 881        CONTINUE
         END DO
 57      CONTINUE
      ENDDO
      ENDDO
 
      DO J=1,NY
      DO I=1,NX
         U1(I,J,1)=USC2(I,J)
         V1(I,J,1)=VSC2(I,J)
      END DO
      END DO

! WRITE 4x area env. data HWRF

      IUNIT=30+ITIM
 
      WRITE(IUNIT) NX,NY,NZ
 
      WRITE(IUNIT) DLMD,DPHD,CENTRAL_LON,CENTRAL_LAT
      WRITE(IUNIT) PT,PDTOP,WBD,SBD
      WRITE(IUNIT) T1
      WRITE(IUNIT) Q1
      WRITE(IUNIT) U1
      WRITE(IUNIT) V1
      WRITE(IUNIT) Z1
!     WRITE(IUNIT) GLON,GLAT
      WRITE(IUNIT) HLON,HLAT,VLON,VLAT
      WRITE(IUNIT) P1
      WRITE(IUNIT) PD1
      WRITE(IUNIT) ETA1
      WRITE(IUNIT) ETA2
      WRITE(IUNIT) USCM    ! no longer use
      WRITE(IUNIT) VSCM
 
      CLOSE(IUNIT)
     
      print*,'NX,NY,NZ=',NX,NY,NZ

      print*,'PDTOP,PT=',PDTOP,PT
 
!     WRITE(61)((SLP1(I,J),I=1,NX),J=1,NY)
!     DO K=1,NZ
!        WRITE(61)((Z1(I,J,K),I=1,NX),J=1,NY)
!     END DO
!     DO K=1,NZ
!        WRITE(61)((USC2(I,J),I=1,NX),J=1,NY)
!     END DO
!     DO K=1,NZ
!        WRITE(61)((VSC2(I,J),I=1,NX),J=1,NY)
!     END DO
!     DO K=1,NZ
!        WRITE(61)((Q1(I,J,K),I=1,NX),J=1,NY)
!     END DO
!     DO K=1,NZ
!        WRITE(61)((U1(I,J,K),I=1,NX),J=1,NY)
!     END DO
!     DO K=1,NZ
!        WRITE(61)((V1(I,J,K),I=1,NX),J=1,NY)
!     END DO
   
      STOP

   END IF !* <---------------------------------------------------------------

      PRINT*, 'Recycled vortex will NOT be enhanced by bogus vortex.' 

   IF(IFLAG_NEST.NE.0)THEN

      imax1=i_max
      jmax1=j_max

      print*,'I,J,vmax=',imax1,jmax1,sqrt(vmax1)
 
      iter=0
      beta=1.0

      vobs=vobs_o/(C101(i_max,j_max)+1.E-10)

      print*,'vobs_o,vobs,C101=',vobs_o,vobs,C101(i_max,j_max)

 876  CONTINUE

      VMAX=0.
!      DO J=1,NYT
!      DO I=1,NXT
!         J1=J+JWMIN1-1+MDY
!         I1=I+IWMIN1-1+MDX
         I1=i_max
         J1=j_max
         UUT=beta*USC1(I1,J1)+USCM(I1,J1)
         VVT=beta*VSC1(I1,J1)+VSCM(I1,J1)
         FF=UUT*UUT+VVT*VVT
!         R05=SQRT((VLAT(I1,J1)-CLAT_NHC)**2+((VLON(I1,J1)-CLON_NHC)*cost)**2)
         R05=dist(i,j)
!         IF(VMAX.LT.FF.and.R05.LT.RMX_d)THEN
            VMAX=FF
            IMV=I1
            JMV=J1
!         END IF
!      END DO
!      END DO

      PRINT*,'I1,J1,I,J,VMAX=',i_max,j_max,IMV,JMV,SQRT(VMAX)

      IMV=I1
      JMV=J1

      UU11=beta*USC1(IMV,JMV)
      VV11=beta*VSC1(IMV,JMV)
      UUM1=USCM(IMV,JMV)
      VVM1=VSCM(IMV,JMV)
      QQ=sqrt((uu11**2+vv11**2)*vobs**2-(vv11*uum1-uu11*vvm1)**2)
        
      beta1=(-(uum1*uu11+vvm1*vv11)+QQ)/(uu11**2+vv11**2+1.E-20)

      print*,'UU11,VV11,UUM1,VVM1,QQ,beta1=',UU11,VV11,UUM1,VVM1,QQ,beta1

      if(beta1.gt.1.0.or.beta1.lt.0.)beta1=1.0

      print*,'UU11,VV11,UUM1,VVM1,QQ,beta1=',UU11,VV11,UUM1,VVM1,QQ,beta1

      beta=beta*beta1
      iter=iter+1
      
      print*,'iter,beta=',iter,beta
 
!     IF(iter.lt.2)go to 876

!     if(beta.gt.1.25) beta=1.25
!     if(beta.lt.0.75) beta=0.75

!     beta=1.0         ! test

!      IF(beta.lt.0.75.and.vobs.lt.24.)THEN
!      IF(beta.lt.0.5.and.vobs.lt.24.)THEN
      IF(beta.lt.0.7.and.vobs.lt.v64kt)THEN
        IFLAG=2
        PRINT*,'IFLAG2=',IFLAG
        OPEN(69,file='flag_file2',form='formatted')
          WRITE(69,*)IFLAG
          write(69,*)K850
!          write(69,*)TWMAX,RWMAX1,fact,psfc_obs 
        CLOSE(69)
        STOP
      END IF

      uv21=sqrt(uu11**2+vv11**2)
      
      uv22=beta*TWMAX

      print*,'max hurricane pert, beta*TWMAX=',uv21,uv22
      print*,'TWMAX=',TWMAX

!     IFLAG=0

!     IF (uv22.LT.(0.5*vobs)) IFLAG=1
!     IF (TWMAX.LT.10.)THEN
!        IF(beta.gt.1.2.or.beta.lt.0.8) IFLAG=1
!     ELSE
!        cut1=min(2.5,1.2+(TWMAX-10.)/15.)
!        cut2=max(0.5,0.8-(TWMAX-10.)/50.)
!        IF(beta.gt.cut1.or.beta.lt.cut2) IFLAG=1
!     END IF 
   
!     beta=1.

      beta=max(beta,10./vmax_s)                ! beta*vmax_s > 10 m/s
      if(beta.gt.1.)beta=1.0

      print*,'new beta,vmax_s=',beta,vmax_s

! now modify the horricane component (by beta)

      T_4=T_1
      Q_4=Q_1    
 
      DO J=1,NY,10
         write(*,33)(HLAT2(I,J),I=1,NX,10)
      END DO
      DO J=1,NY,10
         write(*,33)(VLAT2(I,J),I=1,NX,10)
      END DO
 33   format(15F8.1)

      print*,'CLON_NEW,CLAT_NEW=',CLON_NEW,CLAT_NEW

!     PW=1.0          ! no vertical weighting for storm pert

      ddr = RWM(2)-RWM(1)

      TH1=0.
      do i=1,100
         TH1(i)=TWM(i)
      end do

      do i=1,200
         RP(i)=i*ddr*DST1
      end do

! correct wind
      END IF   ! IFLAG_NEST

      beta=1.


      DO K=1,KMX
      DO J=1,NY
      DO I=1,NX
         U_1(I,J,K)=U_1(I,J,K)*beta
         V_1(I,J,K)=V_1(I,J,K)*beta
      END DO
      END DO
      END DO
!
      DO J=1,NY
      DO I=1,NX
         USC2(I,J)=beta*USC1(I,J)+USCM(I,J)
         VSC2(I,J)=beta*VSC1(I,J)+VSCM(I,J)
      END DO
      END DO

      IF(IFLAG_NEST.NE.0)THEN
!
! correct temp, water vapor and sea-level press

      ics=2
!     fact=1.0

      CALL CORT_MAT_1(IR1,NX,NY,NZ,KMX,                   &
		      T_X,Q_X,SLP_X,HLON,HLAT,VLON,VLAT,  &
		      CLON_NEW,CLAT_NEW,                  &
                      beta,fact,aaa,bbb,iparam,ics) 

      DO J=1,NY
      DO I=1,NX
         SLP_1(I,J)=SLP_1(I,J)*beta+SLP_X(I,J)
      END DO
      END DO

! correct surface press

!     sum1=0.
!     pct_m=0.
!     do j=jmn1,jmx1
!     do i=imn1,imx1
!	 sum1=sum1+1.
!        pct_m=pct_m+SLP_1(I,J)
!     end do
!     end do
!     pct_m=pct_m/sum1
        
      pct_m=0.
      DO J=1,NY
      DO I=1,NX
	 IF (pct_m.GT.SLP_1(I,J)) THEN
	    pct_m=SLP_1(I,J)
	 END IF
      END DO
      END DO
      pct_m=pct_m+50.

!     ps_rat=0.5*(pct_m+psfc_obs1)/pct_m
!     ps_rat=min(pct_m,psfc_obs1)/pct_m
      ps_rat=psfc_obs1/(pct_m-1.E-20)

      print*,'ps_rat=',ps_rat,psfc_obs1,pct_m
      
      if(ps_rat.gt.10.)ps_rat=10.0
      if(ps_rat.lt.(-10.))ps_rat=-10.0
 
      DO J=1,NY
      DO I=1,NX
         SLP_1(I,J)=SLP_1(I,J)*ps_rat
      END DO
      END DO

      DO J=1,NY
      DO I=1,NX
      DO K=1,KMX
	 TEK1=TENV(I,J,K)+T_1(I,J,K)
         T_1(I,J,K)=(T_1(I,J,K)*beta+T_X(I,J,K))*ps_rat
!        Q_1(I,J,K)=Q_1(I,J,K)*beta+Q_X(I,J,K)
	 TEK2=TENV(I,J,K)+T_1(I,J,K)
	 ESRR=exp(4302.645*(TEK2-TEK1)/((TEK2-29.66)*(TEK1-29.66)))
	 Q_1(I,J,K)=ESRR*Q_1(I,J,K)
      END DO
      END DO
      END DO

!     T_1=T_4
!     Q_1=Q_4

!     if(zmax.gt.250.)then
!        T_4=0.
!        Q_4=0.
!     end if

      PRINT*,'complete CORT'  ! --------------------- 

      END IF   ! IFLAG_NEST

      beta=1.

      SLP21=0.
      U21=0.                  ! working array for T_4
      V21=0.                  ! working array for Q_4
      T21=0.
      Q21=0.

      DO J=1,NYT
      DO I=1,NXT
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
         DO N1=1,4
           IH1(N1)=IIH1(I,J,N1)
           JH1(N1)=JJH1(I,J,N1)
         END DO
         DO N1=1,4
            SLP21(I,J,N1)=SLP_1(IH1(N1),JH1(N1))
         END DO
         K=1             ! surface
         DO N1=1,4
            IF (P1(I1,J1,K).GE.PCST(1)) THEN
               U21(I,J,K,N1)=T_4(IH1(N1),JH1(N1),1)         ! old temp
               V21(I,J,K,N1)=Q_4(IH1(N1),JH1(N1),1)         ! old q
               T21(I,J,K,N1)=T_1(IH1(N1),JH1(N1),1)
!              Q21(I,J,K,N1)=Q_1(IH1(N1),JH1(N1),1)
            ELSE IF(P1(I1,J1,K).LE.PCST(KMX))THEN
               U21(I,J,K,N1)=T_4(IH1(N1),JH1(N1),KMX)
               V21(I,J,K,N1)=Q_4(IH1(N1),JH1(N1),KMX)
               T21(I,J,K,N1)=T_1(IH1(N1),JH1(N1),KMX)
!              Q21(I,J,K,N1)=Q_1(IH1(N1),JH1(N1),KMX)
            ELSE
               DO K1=2,KMX
                  PDIF1=P1(I1,J1,K)-PCST(K1)
 !CWH                   if (ZDIF1.GE.0.) then
                  if (PDIF1.GE.0.) then
                     FACT1=PDIF1/(PCST(K1-1)-PCST(K1))
                     U21(I,J,K,N1)=T_4(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                                  +T_4(IH1(N1),JH1(N1),K1-1)*FACT1
                     V21(I,J,K,N1)=Q_4(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                                  +Q_4(IH1(N1),JH1(N1),K1-1)*FACT1
                     T21(I,J,K,N1)=T_1(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
                                  +T_1(IH1(N1),JH1(N1),K1-1)*FACT1
!                    Q21(I,J,K,N1)=Q_1(IH1(N1),JH1(N1),K1)*(1.-FACT1)     &
!                                 +Q_1(IH1(N1),JH1(N1),K1-1)*FACT1
                     GO TO 58
                  end if
               END DO
            END IF
 58         CONTINUE
         END DO
 66      CONTINUE
      ENDDO
      ENDDO

      DO J=1,NYT
      DO I=1,NXT
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
         SLP1(I1,J1) = SLP1(I1,J1)                          &
                   + HBWGT1(I,J,1)*SLP21(I,J,1)             &
                   + HBWGT1(I,J,2)*SLP21(I,J,2)             &
                   + HBWGT1(I,J,3)*SLP21(I,J,3)             &
                   + HBWGT1(I,J,4)*SLP21(I,J,4)
         TENV1       = TS1(I1,J1)
         TS1(I1,J1)  = TENV1                                &
                   + HBWGT1(I,J,1)*T21(I,J,1,1)             &
                   + HBWGT1(I,J,2)*T21(I,J,1,2)             &
                   + HBWGT1(I,J,3)*T21(I,J,1,3)             &
                   + HBWGT1(I,J,4)*T21(I,J,1,4)
         T_OLD       = TENV1                                &
                   + HBWGT1(I,J,1)*U21(I,J,1,1)             &
                   + HBWGT1(I,J,2)*U21(I,J,1,2)             &
                   + HBWGT1(I,J,3)*U21(I,J,1,3)             &
                   + HBWGT1(I,J,4)*U21(I,J,1,4)
         Q_OLD       = QS1(I1,J1)                           &
                   + HBWGT1(I,J,1)*V21(I,J,1,1)             &
                   + HBWGT1(I,J,2)*V21(I,J,1,2)             &
                   + HBWGT1(I,J,3)*V21(I,J,1,3)             &
                   + HBWGT1(I,J,4)*V21(I,J,1,4)
         ESRR        = exp(4302.645*(TS1(I1,J1)-T_OLD)/     &
                     ((TS1(I1,J1)-29.66)*(T_OLD-29.66))) ! 4302.645=17.67*243.5
         QS1(I1,J1)  = Q_OLD + (ESRR-1.)*Q_OLD ! Assume RH=CONST before & after
 67      CONTINUE
      ENDDO
      ENDDO

!     WRITE(25)((SLP1(I,J),I=1,NX),J=1,NY,2)

! based on Ts, Zs, SLP1 ==> PS1  ==> P1

      DO J=1,NY
      DO I=1,NX
         ZSFC = ZS1(I,J)
         TSFC = TS1(I,J)*(1.+D608*QS1(I,J))
         A = (GAMMA * ZSFC) / TSFC
         P1(I,J,1) = SLP1(I,J)/(1+A)**COEF2
!         PD1(I,J)=P1(I,J,1)-PDTOP-PT
         PD1(I,J)=P1(I,J,1)-PT
      ENDDO
      ENDDO
 
! PD(I,J)=P1(I,J,1)-PDTOP-PT=PSFC(I,J)-PDTOP-PT
      DO J=1,NY
      DO I=1,NX
        P1(I,J,NZ1)=PT
      DO K=1,NZ
         P1(I,J,NZ1-K)=P1(I,J,NZ1-K+1)+PDTOP*ETA1(K)+PD1(I,J)*ETA2(K)
      ENDDO
      ENDDO
      ENDDO

      DO K=1,NZ
      DO J=1,NY
      DO I=1,NX
         PMID1(I,J,K)=EXP((ALOG(P1(I,J,K))+ALOG(P1(I,J,K+1)))*0.5)
      ENDDO
      ENDDO
      ENDDO

      do j = 1,ny
      do i = 1,nx
         Z1(I,J,1)=ZS1(I,J)
         DO L=2,nz+1
            Z1(I,J,L)=Z1(I,J,L-1)+T1(I,J,L-1)*              &
                      (Q1(I,J,L-1)*0.608+1.0)*287.04*       &
                      (ALOG(P1(I,J,L-1))-ALOG(P1(I,J,L)))/G
         ENDDO
      ENDDO
      END DO

! interpolate vertically to P level in new coordinate  (V Points)
 
      PMV1=PMID1
 
      DO J=2,NY-1
      DO K=1,NZ
      DO I=2,NX-1
         PMV1(I,J,K)=0.25*(PMID1(I,J,K)+PMID1(I+1,J,K)+    &
                           PMID1(I,J+1,K)+PMID1(I+1,J+1,K))
      END DO
      END DO
      END DO

!     WRITE(63)((SLP1(I,J),I=1,NX),J=1,NY,2)
!     DO K=1,NZ+1
!        WRITE(63)((Z1(I,J,K),I=1,NX),J=1,NY,2)
!     END DO
!     DO K=1,NZ+1
!        WRITE(63)((P1(I,J,K),I=1,NX),J=1,NY,2)
!     END DO
!     DO K=1,NZ
!        WRITE(63)((T1(I,J,K),I=1,NX),J=1,NY,2)
!     END DO
!     DO K=1,NZ
!        WRITE(63)((Q1(I,J,K),I=1,NX),J=1,NY,2)
!     END DO
!     DO K=1,NZ
!        WRITE(63)((U1(I,J,K),I=1,NX),J=1,NY,2)
!     END DO
!     DO K=1,NZ
!        WRITE(63)((V1(I,J,K),I=1,NX),J=1,NY,2)
!     END DO

      PRINT*,'test0'

! add hurricane components
        
      SLP21=0.
      U21=0.
      V21=0.
      T21=0.
      Q21=0.

      DO J=1,NYT
      DO I=1,NXT
         DO N1=1,4
            IH1(N1)=IIH1(I,J,N1)
            JH1(N1)=JJH1(I,J,N1)
         END DO
         DO N1=1,4
            SLP21(I,J,N1)=SLP_1(IH1(N1),JH1(N1))
            DO K=1,KMX
               T21(I,J,K,N1)=T_1(IH1(N1),JH1(N1),K)
               Q21(I,J,K,N1)=Q_4(IH1(N1),JH1(N1),K)
               U21(I,J,K,N1)=T_4(IH1(N1),JH1(N1),K)
!              Q21(I,J,K,N1)=Q_1(IH1(N1),JH1(N1),K)
            END DO
         END DO
 76      CONTINUE
      ENDDO
      ENDDO

      DO J=1,NYT
      DO I=1,NXT
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
         SLP1(I1,J1)=SLP1(I1,J1)       ! already add the storm pert before
!        SLP1(I1,J1)=SLP1(I1,J1)+                      &
!                    HBWGT1(I,J,1)*SLP21(I,J,1)        &
!                  + HBWGT1(I,J,2)*SLP21(I,J,2)        &
!                  + HBWGT1(I,J,3)*SLP21(I,J,3)        &
!                  + HBWGT1(I,J,4)*SLP21(I,J,4)
         DO K=1,KMX
            WRK1(K) =                                  &
                HBWGT1(I,J,1)*T21(I,J,K,1)             &
              + HBWGT1(I,J,2)*T21(I,J,K,2)             &
              + HBWGT1(I,J,3)*T21(I,J,K,3)             &
              + HBWGT1(I,J,4)*T21(I,J,K,4)
            WRK2(K) =                                  &
                HBWGT1(I,J,1)*Q21(I,J,K,1)             &
              + HBWGT1(I,J,2)*Q21(I,J,K,2)             &
              + HBWGT1(I,J,3)*Q21(I,J,K,3)             &
              + HBWGT1(I,J,4)*Q21(I,J,K,4)
            WRK3(K) =                                  &
                HBWGT1(I,J,1)*U21(I,J,K,1)             &
              + HBWGT1(I,J,2)*U21(I,J,K,2)             &
              + HBWGT1(I,J,3)*U21(I,J,K,3)             &
              + HBWGT1(I,J,4)*U21(I,J,K,4)
         END DO
         DO N=1,NZ
            TENV1 = T1(I1,J1,N)
            QENV1 = Q1(I1,J1,N)
            IF (PMID1(I1,J1,N).GE.PCST(1)) THEN            ! Below PCST(1)
               T1(I1,J1,N)=TENV1+WRK1(1)
               Q_OLD=QENV1+WRK2(1)
               T_OLD=TENV1+WRK3(1)
!              Q1(I1,J1,N)=QENV1+WRK2(1)
            ELSE IF(PMID1(I1,J1,N).LE.PCST(KMX)) THEN
               T1(I1,J1,N)=TENV1+WRK1(KMX)
!              Q1(I1,J1,N)=QENV1+WRK2(KMX)
               Q_OLD=QENV1+WRK2(KMX)
               T_OLD=TENV1+WRK3(KMX)
            ELSE
               DO K=1,KMX-1
               if(PMID1(I1,J1,N).LE.PCST(K).and.PMID1(I1,J1,N).GT.PCST(K+1))then
                  W1=ALOG(1.*PCST(K+1))-ALOG(1.*PCST(K))
                  W=(ALOG(1.*PMID1(I1,J1,N))-ALOG(1.*PCST(K)))/W1
                  T1(I1,J1,N)=TENV1+WRK1(K)*(1.-W)+WRK1(K+1)*W
!                 Q1(I1,J1,N)=QENV1+WRK2(K)*(1.-W)+WRK2(K+1)*W
                  Q_OLD=QENV1+WRK2(K)*(1.-W)+WRK2(K+1)*W
                  T_OLD=TENV1+WRK3(K)*(1.-W)+WRK3(K+1)*W
                  GO TO 887
               end if
               END DO
            END IF
 887        CONTINUE
            ESRR = exp(4302.645*(T1(I1,J1,N)-T_OLD)/     &
                  ((T1(I1,J1,N)-29.66)*(T_OLD-29.66)))   ! 4302.645=17.67*243.5
            Q1(I1,J1,N) = Q_OLD+(ESRR-1.)*Q_OLD ! Assume RH=CONST before & after
         END DO
 77      CONTINUE
      ENDDO
      ENDDO

      PRINT*,'test01'

! based on Ts, Zs, SLP1 ==> PS1  ==> P1
      PRINT*,'test02'
 
      DO J=1,NY
      DO I=1,NX
         ZSFC = ZS1(I,J)
         TSFC = TS1(I,J)*(1.+D608*QS1(I,J))
         A = (GAMMA * ZSFC) / TSFC
         P1(I,J,1) = SLP1(I,J)/(1+A)**COEF2
!         PD1(I,J)=P1(I,J,1)-PDTOP-PT
         PD1(I,J)=P1(I,J,1)-PT
      ENDDO
      ENDDO
 
! PD(I,J)=P1(I,J,1)-PDTOP-PT=PSFC(I,J)-PDTOP-PT
      DO J=1,NY
      DO I=1,NX
        P1(I,J,NZ1)=PT
      DO K=1,NZ
         P1(I,J,NZ1-K)=P1(I,J,NZ1-K+1)+PDTOP*ETA1(K)+PD1(I,J)*ETA2(K)
      ENDDO
      ENDDO
      ENDDO

      DO j = 1,ny
      DO i = 1,nx
         Z1(I,J,1)=ZS1(I,J)
         DO L=2,nz+1
            Z1(I,J,L)=Z1(I,J,L-1)+T1(I,J,L-1)*              &
                      (Q1(I,J,L-1)*0.608+1.0)*287.04*       &
                      (ALOG(P1(I,J,L-1))-ALOG(P1(I,J,L)))/G
         ENDDO
      ENDDO
      ENDDO
        
      U21=0.
      V21=0.
 
      DO J=1,NYT
      DO I=1,NXT
         DO N1=1,4
            IV1(N1)=IIV1(I,J,N1)
            JV1(N1)=JJV1(I,J,N1)
         END DO
         DO N1=1,4
         DO K=1,KMX
            U21(I,J,K,N1)=U_1(IV1(N1),JV1(N1),K)
            V21(I,J,K,N1)=V_1(IV1(N1),JV1(N1),K)
         END DO
         END DO
 86      CONTINUE
      ENDDO
      ENDDO

      PRINT*,'test03'

      DO J=1,NYT
      DO I=1,NXT
         J1=J+JWMIN1-1+MDY
         I1=I+IWMIN1-1+MDX
         DO K=1,KMX
            WRK1(K) =                          &
               VBWGT1(I,J,1)*U21(I,J,K,1)      &
             + VBWGT1(I,J,2)*U21(I,J,K,2)      &
             + VBWGT1(I,J,3)*U21(I,J,K,3)      &
             + VBWGT1(I,J,4)*U21(I,J,K,4)
            WRK2(K) =                          &
               VBWGT1(I,J,1)*V21(I,J,K,1)      &
             + VBWGT1(I,J,2)*V21(I,J,K,2)      &
             + VBWGT1(I,J,3)*V21(I,J,K,3)      &
             + VBWGT1(I,J,4)*V21(I,J,K,4)
         END DO

         DO N=1,NZ

            IF (PMV1(I1,J1,N).GE.PCST(1)) THEN            ! Below PCST(1)
               U1(I1,J1,N)=U1(I1,J1,N)+WRK1(1)
               V1(I1,J1,N)=V1(I1,J1,N)+WRK2(1)
            ELSE IF (PMV1(I1,J1,N).LE.PCST(KMX)) THEN
               U1(I1,J1,N)=U1(I1,J1,N)+WRK1(KMX)
               V1(I1,J1,N)=V1(I1,J1,N)+WRK2(KMX)
            ELSE
               DO K=1,KMX-1
               if(PMV1(I1,J1,N).LE.PCST(K).and.PMV1(I1,J1,N).GT.PCST(K+1))then
                  W1=ALOG(1.*PCST(K+1))-ALOG(1.*PCST(K))
                  W=(ALOG(1.*PMV1(I1,J1,N))-ALOG(1.*PCST(K)))/W1
                  U1(I1,J1,N)=U1(I1,J1,N)+WRK1(K)*(1.-W)+WRK1(K+1)*W
                  V1(I1,J1,N)=V1(I1,J1,N)+WRK2(K)*(1.-W)+WRK2(K+1)*W
                  GO TO 888
               end if
               END DO
            END IF
 888        CONTINUE
         END DO
 87      CONTINUE
      ENDDO
      ENDDO

      DEALLOCATE ( IIH1, JJH1, IIV1, JJV1, HBWGT1, VBWGT1 )

      DO J=1,NY
      DO I=1,NX
         U1(I,J,1)=USC2(I,J)
         V1(I,J,1)=VSC2(I,J)
      END DO
      END DO
     print*,'anl_4x_step2 fort.64 NX,NY= ', NX,NY
 
!     WRITE(64)((SLP1(I,J),I=1,NX),J=1,NY)
!     DO K=1,NZ+1
!        WRITE(64)((Z1(I,J,K),I=1,NX),J=1,NY)
!     END DO
!     DO K=1,NZ+1
!        WRITE(64)((P1(I,J,K),I=1,NX),J=1,NY)
!     END DO
!     DO K=1,NZ
!        WRITE(64)((T1(I,J,K),I=1,NX),J=1,NY)
!     END DO
!     DO K=1,NZ
!        WRITE(64)((Q1(I,J,K),I=1,NX),J=1,NY)
!     END DO
!     DO K=1,NZ
!        WRITE(64)((U1(I,J,K),I=1,NX),J=1,NY)
!     END DO
!     DO K=1,NZ
!        WRITE(64)((V1(I,J,K),I=1,NX),J=1,NY)
!     END DO
!     WRITE(64)((USCM(I,J),I=1,NX),J=1,NY)
!     WRITE(64)((VSCM(I,J),I=1,NX),J=1,NY)

      IUNIT=50+ITIM
 
      WRITE(IUNIT) NX,NY,NZ
      WRITE(IUNIT) DLMD,DPHD,CENTRAL_LON,CENTRAL_LAT
      WRITE(IUNIT) PT,PDTOP,WBD,SBD
      WRITE(IUNIT) T1
      WRITE(IUNIT) Q1
      WRITE(IUNIT) U1
      WRITE(IUNIT) V1
      WRITE(IUNIT) Z1
!     WRITE(IUNIT) GLON,GLAT
      WRITE(IUNIT) HLON,HLAT,VLON,VLAT
      WRITE(IUNIT) P1
      WRITE(IUNIT) PD1
      WRITE(IUNIT) ETA1
      WRITE(IUNIT) ETA2

      CLOSE(IUNIT)

      END 

!==============================================================================

      SUBROUTINE FIND_NEWCT1(IX,JX,UD,VD,GLON2,GLAT2,    &
                             CLON_NEW1,CLAT_NEW1)
                                                    
!     PARAMETER (IR=100,IT=24,IX=254,JX=254)
      PARAMETER (IR=60,IT=24)
      PARAMETER (ID=31,JD=31,DTX=0.1,DTY=0.1)    ! Search x-Domain (ID-1)*DTX
      REAL (4) UD(IX,JX),VD(IX,JX),GLON2(IX,JX),GLAT2(IX,JX)
!     DIMENSION RWM(IR+1),TWM(IR+1)
      DIMENSION TNMX(ID,JD),RX(ID,JD),WTM(IR)
      REAL (4) CLON_NEW1,CLAT_NEW1

      PI=ASIN(1.)*2.
      RAD=PI/180.

      ddr=0.1

      pi180=RAD
      cost=cos(clat_new*pi180)

      ix2=ix/2
      jx2=jx/2
      DDS=(((GLON2(ix2+1,jx2)-GLON2(ix2,jx2))*cost)**2+     &
          (GLAT2(ix2,jx2+1)-GLAT2(ix2,jx2))**2)*1.5
 
 
       print*,'ix,jx,ix2,jx2=',ix,jx,ix2,jx2
       print*,'CLON_NEW,CLAT_NEW=',CLON_NEW1,CLAT_NEW1
       print*,'GLON2,GLAT2=',GLON2(1,1),GLAT2(1,1)
 
 
      XLAT = CLAT_NEW1-(JD-1)*DTY/2.
      XLON = CLON_NEW1-(ID-1)*DTX/2.
 
!     print *,'STARTING LAT, LON AT FIND NEW CENTER ',XLAT,XLON

      DO J=1,JD
      DO I=1,ID
      TNMX(I,J) = 0.
      RX(i,j)=0.
      BLON = XLON + (I-1)*DTX
      BLAT = XLAT + (J-1)*DTY

!..   CALCULATE TANGENTIAL WIND EVERY 0.2 deg INTERVAL
!..   10*10 deg AROUND 1ST GUESS VORTEX CENTER

      DO 10 JL=1,IR
      WTS= 0.
      DO 20 IL=1,IT
      DR = JL*ddr
!     DR = JL
      DD = (IL-1)*15*RAD
      DLON = DR*COS(DD)
      DLAT = DR*SIN(DD)
      TLON = BLON + DLON
      TLAT = BLAT + DLAT
       
!..   INTERPOLATION U, V AT TLON,TLAT AND CLACULATE TANGENTIAL WIND
       
      u1=0.
      v1=0.
      sum1=0.
      DO j1=1,JX
      DO i1=1,IX
         dist=(((GLON2(i1,j1)-TLON)*cost)**2+(GLAT2(i1,j1)-TLAT)**2)
         IF (dist.LT.DDS) THEN
            dist1=1./dist
            sum1=sum1+dist1
            u1=u1+UD(i1,j1)*dist1
            v1=v1+VD(i1,j1)*dist1
         END IF
      END DO
      END DO

      UT=u1/sum1
      VT=v1/sum1
       
!..   TANGENTIAL WIND
      WT = -SIN(DD)*UT + COS(DD)*VT
      WTS = WTS+WT
20    CONTINUE
      WTM(JL) = WTS/24.
10    CONTINUE

!     Southern Hemisphere
      IF (CLAT_NEW.LT.0) THEN
         DO JL=1,IR
            WTM(JL)=-WTM(JL)
         END DO
      END IF
! EnD SH
       
!     print*,'test1'
       
      TX = -10000000.
      DO KL = 1,IR
      IF(WTM(KL).GE.TX) THEN
      TX = WTM(KL)
      RRX = KL*ddr
      ENDIF
      ENDDO
!        DO KL=1,IR
!           TWM(KL)=WTM(KL)
!           RWM(KL)=KL*ddr
!        END DO
!        TWM(IR+1)=TX
!        RWM(IR+1)=RRX

      TNMX(I,J) = TX
      RX(I,J)=RRX
      ENDDO
      ENDDO
!C..  FIND NEW CENTER
      TTX = -1000000.
      DO I=1,ID
      DO J=1,JD
      IF(TNMX(I,J).GE.TTX) THEN
      TTX = TNMX(I,J)
      NIC = I
      NJC = J
      ENDIF
      ENDDO
      ENDDO
       
! QLIU test
!     print*,XLAT+30*DTY,XLON+30*DTX,TNMX(30,30)
      print*,'max WTM=',TNMX(30,30),RX(30,30)

      CLAT_NEW1 = XLAT + (NJC-1)*DTY
      CLON_NEW1 = XLON + (NIC-1)*DTX

!     print *,'NEW CENTER,  I, J IS   ',NIC,NJC
      print *,'NEW CENTER, LAT,LON IS ',CLAT_NEW1,CLON_NEW1
!     print *,'MAX TAN. WIND AT NEW CENTER IS ',TTX

      RETURN
      END

!==============================================================================
