!-----------------------------------------------------------------------
!
      PROGRAM CONVERT_LATLON_TO_IJ
!
!-----------------------------------------------------------------------
!
!***  Given the geographic latitude and longitude on an Arakawa B grid,
!***  compute the nearest mass point's (I,J).
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
!-----------------------------------------------------------------------
      INTEGER, PARAMETER :: DOUBLE=SELECTED_REAL_KIND(P=13,R=200)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!***  In the following line insert the desired latitude/longitude 
!***  on an B grid mass point.
!
!-----------------------------------------------------------------------
!
      REAL :: GLATD                                   !<-- Geographic latitude, degrees, positive north
      REAL :: GLOND                                   !<-- Geographic longitude, degrees, positive east
!
!-----------------------------------------------------------------------
!
!***  In the following lines set:
!***    (1) The global I and J extent of the grid
!***    (2) The geographic latitude (TPH0D) (degrees, positive north)
!***        and geographic longitude (TLM0D) (degrees, positive east) 
!***        for the central point on the grid.
!***    (3) The angular distance (degrees) from the domain's central
!***        point to the western boundary (WBD) and the same for
!***        the southern boundary (SBD).  Both values will be <0.
!***        of the domain.
!
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: IM=954,JM=835                                   !<-- Full grid dimensions
!
      REAL(kind=DOUBLE),PARAMETER :: TPH0D=54.00,TLM0D=-106.00             !<-- Central lat/lon
!
      REAL(kind=DOUBLE),PARAMETER :: WBD=-60.039                        &  !<-- Rotated longitude of western boundary
                                    ,SBD=-45.036                           !<-- Rotated latitude of southern boundary
!
      REAL(kind=DOUBLE),PARAMETER :: DPHD=-2.*SBD/(JM-1)                &  !<-- Rotated latitude grid increment (degrees)
                                    ,DLMD=-2.*WBD/(IM-1)                   !<-- Rotated longitude grid increment (degrees)
!
!--------------------
!*** Local Variables
!--------------------
!
      INTEGER :: I_NEAR,J_NEAR,N,N_SAVE,NCOL,NROW,INSIDE
!
      REAL(KIND=DOUBLE) :: ALPHA,ARG,BETA,COL,CROSS                     &
                          ,D2R,DIST_SAVE,DLM,DLON,DPH,GLAT,GLON         &
     &                    ,ONE,PI,PI_H,R2D,ROW,SB,TLAT0,TLON0           &
     &                    ,TPH0,TLM0,WB,X,Y,Z
!
      REAL(KIND=DOUBLE),DIMENSION(4) :: DISTANCE,TLAT,TLON
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  
!***  Convert from geographic to transformed coordinates (degrees).
!***
!-----------------------------------------------------------------------
!
      open(5,file='itag')
      read(5,*) GLATD
      read(5,*) GLOND
!
      ONE=1.
      PI=DACOS(-ONE)
      PI_H=ACOS(0.)
      D2R=PI/180.
      R2D=1./D2R
!
      GLAT=GLATD*D2R
      GLON=GLOND*D2R
      DPH=DPHD*D2R
      DLM=DLMD*D2R
      TPH0=TPH0D*D2R
      TLM0=TLM0D*D2R
      SB=SBD*D2R
      WB=WBD*D2R
!
      X=COS(TPH0)*COS(GLAT)*COS(GLON-TLM0)+SIN(TPH0)*SIN(GLAT)
      Y=COS(GLAT)*SIN(GLON-TLM0)
      Z=COS(TPH0)*SIN(GLAT)-SIN(TPH0)*COS(GLAT)*COS(GLON-TLM0)
      TLAT0=R2D*ATAN(Z/SQRT(X*X+Y*Y))
      TLON0=R2D*ATAN(Y/X)
!
      WRITE(0,50)TLAT0,TLON0
   50 FORMAT(' Transformed Latitude is',F8.3                            &
            ,4X,'Longitude is',F8.3)
!
!-----------------------------------------------------------------------
!***  REAL row and column of the location.
!-----------------------------------------------------------------------
!
      ROW=TLAT0/DPHD+0.5*(JM+1)
      COL=TLON0/DLMD+0.5*(IM+1)
!
      TLAT0=TLAT0*D2R
      TLON0=TLON0*D2R
!
!-----------------------------------------------------------------------
!***
!***     H3     H4
!***
!***
!***        H0
!***     H1     H2
!***
!-----------------------------------------------------------------------
!***  Rotated lat/lon of the four surrounding points.
!-----------------------------------------------------------------------
!
      NROW=MAX(INT(ROW),1)
      NCOL=MAX(INT(COL),1)
!
      TLAT(1)=(NROW-1)*DPH+SB
      TLAT(2)=TLAT(1)
      TLAT(3)=TLAT(1)+DPH
      TLAT(4)=TLAT(3)
      TLON(1)=(NCOL-1)*DLM+WB
      TLON(2)=TLON(1)+DLM
      TLON(3)=TLON(1)
      TLON(4)=TLON(2)
!
!-----------------------------------------------------------------------
!***  Distance to the four surrounding points.
!-----------------------------------------------------------------------
!
      N_SAVE=0
      DIST_SAVE=1.E10
!
      DO N=1,4
!
        DLON=TLON(N)-TLON0
        CROSS=ACOS(COS(DLON)*COS(TLAT(N)))
        ARG=TAN(TLAT(N))/SIN(DLON)
        ALPHA=ATAN(ARG)
        IF(DLON<0.)ALPHA=-ALPHA
        BETA=PI_H-ALPHA
        DISTANCE(N)=ACOS(COS(TLAT0)*COS(TLAT(N))*COS(DLON)              &
                        +SIN(TLAT0)*SIN(CROSS)*COS(BETA))
        IF(DISTANCE(N)<DIST_SAVE)THEN
          DIST_SAVE=DISTANCE(N)
          N_SAVE=N
        ENDIF
      ENDDO
!
!-----------------------------------------------------------------------
!
      I_NEAR=NINT((TLON(N_SAVE)-WB)/DLM+1)
      J_NEAR=NINT((TLAT(N_SAVE)-SB)/DPH+1)
!
      WRITE(0,100)I_NEAR,J_NEAR
  100 FORMAT(' Nearest Mass Point At I=',I4,' J=',I4)
!
      INSIDE=0
      IF(I_NEAR.GT.10 .AND. I_NEAR.LE.IM-10) THEN
        INSIDE=1
      ELSE
        INSIDE=0
        GO TO 102
      ENDIF
      IF(J_NEAR.GT.10 .AND. J_NEAR.LE.JM-10) THEN
        INSIDE=1
      ELSE
        INSIDE=0
      ENDIF
  102 CONTINUE
      WRITE(0,101)INSIDE
  101 FORMAT(' INSIDE = ',I1)
      WRITE(51,103) INSIDE
  103 FORMAT(I1)
!
!-----------------------------------------------------------------------
!
      END PROGRAM CONVERT_LATLON_TO_IJ
!
!-----------------------------------------------------------------------
