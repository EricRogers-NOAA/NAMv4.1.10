MODULE INTERP_COEF_BGRID

CONTAINS

SUBROUTINE EARTH_LATLON_BGRID ( HLAT,HLON,VLAT,VLON,     & !Earth lat,lon at H and V points
                                DLMD,DPHD,WBD,SBD,       & !input res,west & south boundaries,
                                CENTRAL_LAT,CENTRAL_LON, & ! central lat,lon, all in degrees   
                                IM,JM)
!
!============================================================================
!
 IMPLICIT NONE
 INTEGER,    INTENT(IN)                 :: IM,JM
 REAL(4),    INTENT(IN)                 :: DLMD,DPHD,WBD,SBD
 REAL(4),    INTENT(IN)                 :: CENTRAL_LAT,CENTRAL_LON
 REAL(4), DIMENSION(IM,JM), INTENT(OUT) :: HLAT,HLON,VLAT,VLON

! local

 INTEGER,PARAMETER                 :: KNUM=SELECTED_REAL_KIND(13) 
 INTEGER                           :: I,J
 REAL(KIND=KNUM)                   :: WB,SB,DLM,DPH,TPH0,STPH0,CTPH0
 REAL(KIND=KNUM)                   :: TLMH,TLMV,TLMH0,TLMV0,TPHH,TPHV,DTR
 REAL(KIND=KNUM)                   :: STPH,CTPH,STPV,CTPV,PI_2
 REAL(KIND=KNUM)                   :: SPHH,CLMH,FACTH,SPHV,CLMV,FACTV
 REAL(KIND=KNUM), DIMENSION(IM,JM) :: GLATH,GLONH,GLATV,GLONV
!-------------------------------------------------------------------------

      PI_2 = ACOS(0.)
      DTR  = PI_2/90.
      WB   = WBD * DTR                 ! WB:   western boundary in radians
      SB   = SBD * DTR                 ! SB:   southern boundary in radians
      DLM  = DLMD * DTR                ! DLM:  dlamda in radians 
      DPH  = DPHD * DTR                ! DPH:  dphi   in radians

!     For earth lat lon only

      TPH0  = CENTRAL_LAT*DTR          ! TPH0: central lat in radians 
      STPH0 = SIN(TPH0)
      CTPH0 = COS(TPH0)

      DO J = 1,JM

         TLMH0 = WB - DLM
         TLMV0 = WB - DLM + 0.5*DLM
         TPHH = SB + (J-1)*DPH
         TPHV = TPHH + 0.5*DPH
         STPH = SIN(TPHH)
         CTPH = COS(TPHH)
         STPV = SIN(TPHV)
         CTPV = COS(TPHV)

         DO I = 1,IM
           TLMH = TLMH0 + I*DLM

           SPHH = CTPH0 * STPH + STPH0 * CTPH * COS(TLMH)
           GLATH(I,J)=ASIN(SPHH)
           CLMH = CTPH*COS(TLMH)/(COS(GLATH(I,J))*CTPH0) &
                - TAN(GLATH(I,J))*TAN(TPH0)
           IF(CLMH .GT. 1.) CLMH = 1.0
           IF(CLMH .LT. -1.) CLMH = -1.0
           FACTH = 1.
           IF(TLMH .GT. 0.) FACTH = -1.
           GLONH(I,J) = -CENTRAL_LON*DTR + FACTH*ACOS(CLMH)

         ENDDO                                    

         DO I = 1,IM
           TLMV = TLMV0 + I*DLM

           SPHV = CTPH0 * STPV + STPH0 * CTPV * COS(TLMV)
           GLATV(I,J) = ASIN(SPHV)
           CLMV = CTPV*COS(TLMV)/(COS(GLATV(I,J))*CTPH0) &
                - TAN(GLATV(I,J))*TAN(TPH0)
           IF(CLMV .GT. 1.) CLMV = 1.
           IF(CLMV .LT. -1.) CLMV = -1.
           FACTV = 1.
           IF(TLMV .GT. 0.) FACTV = -1.
           GLONV(I,J) = -CENTRAL_LON*DTR + FACTV*ACOS(CLMV)

         ENDDO

      ENDDO

!     Conversion to degrees (may not be required, eventually)

      DO J = 1, JM
       DO I = 1, IM
          HLAT(I,J) =  GLATH(I,J) / DTR
          HLON(I,J) = -GLONH(I,J) / DTR
          IF(HLON(I,J) .GT. 60.) HLON(I,J) = HLON(I,J)  - 360.
          IF(HLON(I,J) .LT. -300.) HLON(I,J) = HLON(I,J) + 360.
!          IF(HLON(I,J) .GT. 180.) HLON(I,J) = HLON(I,J)  - 360.
!          IF(HLON(I,J) .LT. -180.) HLON(I,J) = HLON(I,J) + 360.
!
          VLAT(I,J) =  GLATV(I,J) / DTR
          VLON(I,J) = -GLONV(I,J) / DTR
          IF(VLON(I,J) .GT. 60.) VLON(I,J) = VLON(I,J)  - 360.
          IF(VLON(I,J) .LT. -300.) VLON(I,J) = VLON(I,J) + 360.
!          IF(VLON(I,J) .GT. 180.) VLON(I,J) = VLON(I,J)  - 360.
!          IF(VLON(I,J) .LT. -180.) VLON(I,J) = VLON(I,J) + 360.

       ENDDO
      ENDDO

END SUBROUTINE EARTH_LATLON_BGRID

!-----------------------------------------------------------------------------  

SUBROUTINE G2T2H_BGRID( IIH,JJH,                     & ! output grid index
                        HBWGT,                       & ! output weights in terms of parent grid
                        HLAT,HLON,                   & ! target (nest) input lat lon in degrees
                        DLMD1,DPHD1,WBD1,SBD1,       & ! parent res, west and south boundaries
                        CENTRAL_LAT,CENTRAL_LON,     & ! parent central lat,lon, all in degrees
                        P_IM,P_JM,                   & ! parent imax and jmax
                        IM,JM)                         ! target (nest) dimensions
! 
!============================================================================
!
 IMPLICIT NONE
 INTEGER,    INTENT(IN)                       :: IM,JM
 INTEGER,    INTENT(IN)                       :: P_IM,P_JM
 REAL(4),    INTENT(IN)                       :: DLMD1,DPHD1,WBD1,SBD1
 REAL(4),    INTENT(IN)                       :: CENTRAL_LAT,CENTRAL_LON
 REAL(4),    DIMENSION(IM,JM),   INTENT(IN)   :: HLAT,HLON
 REAL(4),    DIMENSION(IM,JM,4), INTENT(OUT)  :: HBWGT
 INTEGER,    DIMENSION(IM,JM,4), INTENT(OUT)  :: IIH,JJH
! local

 INTEGER                           :: I,J
 INTEGER                           :: I1,I2,J1,J2
 REAL                              :: X,Y,XI,YI,XF,YF
!-------------------------------------------------------------------------

      DO J = 1,JM
       DO I = 1,IM

         CALL TLL(HLON(I,J),HLAT(I,J),X,Y,CENTRAL_LAT,CENTRAL_LON)
         
         XI = (X-WBD1)/DLMD1 + 1
         YI = (Y-SBD1)/DPHD1 + 1
         I1 = XI
         I2 = I1+1
         J1 = YI
         J2 = J1+1
         XF=XI-I1
         YF=YI-J1
         IIH(I,J,1)=I1
         IIH(I,J,2)=I2
         IIH(I,J,3)=I1
         IIH(I,J,4)=I2
         JJH(I,J,1)=J1
         JJH(I,J,2)=J1
         JJH(I,J,3)=J2
         JJH(I,J,4)=J2
         HBWGT(I,J,1)=(1-XF)*(1-YF)
         HBWGT(I,J,2)=XF*(1-YF)
         HBWGT(I,J,3)=(1-XF)*YF
         HBWGT(I,J,4)=XF*YF

       ENDDO
      ENDDO

END SUBROUTINE G2T2H_BGRID



SUBROUTINE G2T2V_BGRID( IIV,JJV,                     & ! output grid index and weights
                        VBWGT,                       & ! output weights in terms of parent grid
                        VLAT,VLON,                   & ! target (nest) input lat lon in degrees
                        DLMD1,DPHD1,WBD1,SBD1,       & ! parent res, west and south boundaries
                        CENTRAL_LAT,CENTRAL_LON,     & ! parent central lat,lon, all in degrees
                        P_IM,P_JM,                   & ! parent imax and jmax
                        IM,JM)                         ! target (nest) dimensions
!
!============================================================================
!
 IMPLICIT NONE
 INTEGER,    INTENT(IN)                       :: IM,JM
 INTEGER,    INTENT(IN)                       :: P_IM,P_JM
 REAL(4),    INTENT(IN)                       :: DLMD1,DPHD1,WBD1,SBD1
 REAL(4),    INTENT(IN)                       :: CENTRAL_LAT,CENTRAL_LON
 REAL(4),    DIMENSION(IM,JM),   INTENT(IN)   :: VLAT,VLON
 REAL(4),    DIMENSION(IM,JM,4), INTENT(OUT)  :: VBWGT
 INTEGER,    DIMENSION(IM,JM,4), INTENT(OUT)  :: IIV,JJV

! local

 INTEGER                           :: I,J
 INTEGER                           :: I1,I2,J1,J2
 REAL                              :: X,Y,XI,YI,XF,YF
!-------------------------------------------------------------------------

      DO J = 1,JM
       DO I = 1,IM
 
         CALL TLL(VLON(I,J),VLAT(I,J),X,Y,CENTRAL_LAT,CENTRAL_LON)
         
         XI = (X-WBD1-0.5*DLMD1)/DLMD1 + 1
         YI = (Y-SBD1-0.5*DPHD1)/DPHD1 + 1
         I1 = XI
         I2 = I1+1
         J1 = YI
         J2 = J1+1
         XF=XI-I1
         YF=YI-J1
         IIV(I,J,1)=I1
         IIV(I,J,2)=I2
         IIV(I,J,3)=I1
         IIV(I,J,4)=I2
         JJV(I,J,1)=J1
         JJV(I,J,2)=J1
         JJV(I,J,3)=J2
         JJV(I,J,4)=J2
         VBWGT(I,J,1)=(1-XF)*(1-YF)
         VBWGT(I,J,2)=XF*(1-YF)
         VBWGT(I,J,3)=(1-XF)*YF
         VBWGT(I,J,4)=XF*YF
         
       ENDDO
      ENDDO


 RETURN
 END SUBROUTINE G2T2V_BGRID

   subroutine tll(almd,aphd,tlmd,tphd,tph0d,tlm0d)
!-------------------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------------------
      real, intent(in) :: almd, aphd
      real, intent(out) :: tlmd, tphd
      real, intent(in) :: tph0d, tlm0d
!-------------------------------------------------------------------------------
      real, parameter :: pi=3.141592654
      real, parameter :: dtr=pi/180.0
!
      real :: tph0, ctph0, stph0, relm, srlm, crlm
      real :: aph, sph, cph, cc, anum, denom
!-------------------------------------------------------------------------------
!
      if (tlm0d==0.0.and.tph0d==0.0) then
      tlmd=almd
      tphd=aphd
      else

      tph0=tph0d*dtr
      ctph0=cos(tph0)
      stph0=sin(tph0)
!
      relm=(almd-tlm0d)*dtr
      srlm=sin(relm)
      crlm=cos(relm)
      aph=aphd*dtr
      sph=sin(aph)
      cph=cos(aph)
      cc=cph*crlm
      anum=cph*srlm
      denom=ctph0*cc+stph0*sph
!
      tlmd=atan2(anum,denom)/dtr
      tphd=asin(ctph0*sph-stph0*cc)/dtr

      end if
!
      return
!
   end subroutine tll

END MODULE INTERP_COEF_BGRID
