!??????????????????????????????????????????????????????????
      SUBROUTINE CREAT_41X41(ITIM,KST,KMX,MTV6,KS850,U850,V850,
     &                       SDAT,P2)
! SUBPROGRAM    
!   PRGRMMR 
!     
! ABSTRACT
!     
!     DECLARE VARIABLES
!
      INTEGER I,J,K,NX,NY,NZ
!
!      PARAMETER (NX=215,NY=431,NZ=42)
      PARAMETER (NST=10,IRX=41,JRX=41,MAXNGRID=10000)
      PARAMETER (GAMMA=6.5E-3,G=9.8,GI=1./G,D608=0.608)
!     
! variable for outer nest

      REAL(4), ALLOCATABLE :: T3(:,:,:),Q3(:,:,:)
      REAL(4), ALLOCATABLE :: U3(:,:,:),V3(:,:,:)
      REAL(4), ALLOCATABLE :: Z3(:,:,:),P3(:,:,:)
      REAL(4), ALLOCATABLE :: SLP3(:,:),HP3(:,:,:),RH3(:,:)

      REAL(4), ALLOCATABLE :: ZM3(:,:,:),PM3(:,:,:)
      REAL(4), ALLOCATABLE :: PMV3(:,:,:),TV3(:,:,:) 
                                                                                                                                
      REAL(4), ALLOCATABLE :: HLON3(:,:),HLAT3(:,:)
      REAL(4), ALLOCATABLE :: VLON3(:,:),VLAT3(:,:)

      REAL(4) DLMD3,DPHD3,CLON3,CLAT3,PT3,PDTOP3

      real(4),allocatable :: UT850(:,:),VT850(:,:)
      real(4),allocatable :: wrk1(:,:,:),wrk2(:,:,:)

      integer(4) IGD(IRX,JRX,MAXNGRID),JGD(IRX,JRX,MAXNGRID)
      integer(4) NSUM(IRX,JRX),NGRID,KST
      real  WGD(IRX,JRX,MAXNGRID),WSUM(IRX,JRX)

      REAL(4) BLON(IRX), BLAT(JRX) 

      REAL(4) U850(IRX,JRX),V850(IRX,JRX)

      REAL(4) P2(KMX),SDAT(IRX,JRX,MTV6)

! Variables on P surface 

      real(4) SLON_N,SLAT_N,CLON_N,CLAT_N
                                                                                                       
      COMMON /NHC/ KSTM,IC_N(NST),JC_N(NST)
      COMMON /NHC1/ SLON_N(NST),SLAT_N(NST),CLON_N(NST),CLAT_N(NST)
      CHARACTER ST_NAME(NST)*3,STMNAME(NST)*3,TCVT(NST)*95
      COMMON /STNAME/ST_NAME,STMNAME
      COMMON /TCVIT/TCVT
      COMMON /RSFC/STRPSF(NST),STVMAX(NST)
!    
      COEF3=287.05*GI*GAMMA
      COEF2=1./COEF3

      pi=4.*atan(1.)
      pi_deg=180./pi 
      pi180=pi/180.
! read data
    
!      READ(5,*)ITIM,IBGS

! newly added
! compute SDATA from outter nest domain

      IUNIT=40+ITIM

      READ(IUNIT)KX,KY,KZ
                                                                                                                 
      KZ1=KZ+1
                                                                                                                 
      print*,'KX,KY,KZ=',KX,KY,KZ

      max_kz_kmx=kz
      if(kmx .gt. kz) then
         max_kz_kmx=kmx
      end if

      ALLOCATE ( HLON3(KX,KY),HLAT3(KX,KY) )
      ALLOCATE ( VLON3(KX,KY),VLAT3(KX,KY) )
                                                                                                                 
      ALLOCATE ( T3(KX,KY,KZ),Q3(KX,KY,KZ) )
      ALLOCATE ( U3(KX,KY,KZ),V3(KX,KY,KZ) )
      ALLOCATE ( Z3(KX,KY,KZ1),P3(KX,KY,KZ1) )
                                                                                                                 
      ALLOCATE ( SLP3(KX,KY),ZM3(KX,KY,max_kz_kmx) )
      ALLOCATE ( HP3(KX,KY,max_kz_kmx),RH3(KX,KY) )
      ALLOCATE ( PM3(KX,KY,max_kz_kmx),PMV3(KX,KY,max_kz_kmx) )
      ALLOCATE ( TV3(KX,KY,max_kz_kmx) )

      ALLOCATE ( wrk1(KX,KY,KMX), wrk2(KX,KY,KMX) )

      ALLOCATE ( UT850(KX,KY),VT850(KX,KY) )
                                                                                                                 
      READ(IUNIT) DLMD3,DPHD3,CLON3,CLAT3
      READ(IUNIT) PT3,PDTOP3
      READ(IUNIT) T3
      READ(IUNIT) Q3
      READ(IUNIT) U3
      READ(IUNIT) V3
      READ(IUNIT) Z3
      READ(IUNIT) HLON3,HLAT3,VLON3,VLAT3
      READ(IUNIT) P3
      READ(IUNIT) 
      READ(IUNIT) 
      READ(IUNIT) 
                                                                                                                 
                                                                                                                 
      CLOSE(IUNIT)

      DO J=1,KY
      DO I=1,KX
	IF(HLON3(I,J).GT.60.)HLON3(I,J)=HLON3(I,J)-360.
	IF(VLON3(I,J).GT.60.)VLON3(I,J)=VLON3(I,J)-360.
      END DO
      END DO

      print*,'finish reading data'

! LON & LAT at U,V
                                                                                                                                            
      IF(HLON3(1,1).LT.0.)THEN
        DO J=1,KY
        DO I=1,KX
          HLON3(I,J)=HLON3(I,J)+360.
        END DO
        END DO
      END IF
      IF(VLON3(1,1).LT.0.)THEN
        DO J=1,KY
        DO I=1,KX
          VLON3(I,J)=VLON3(I,J)+360.
        END DO
        END DO
      END IF
!!

      DO K=1,KZ
      DO J=1,KY
      DO I=1,KX
        ZM3(I,J,K)=(Z3(I,J,K)+Z3(I,J,K+1))*0.5
        PM3(I,J,K)=EXP((ALOG(1.*P3(I,J,K))+ALOG(1.*P3(I,J,K+1)))*0.5)
        TV3(I,J,K)=T3(I,J,K)*(1.+D608*Q3(I,J,K))
      END DO
      END DO
      END DO

!C        COMPUTE SEA LEVEL PRESSURE.
!C
       DO J=1,KY
       DO I=1,KX
         ZSF1 = ZM3(I,J,1)
         PSF1 = PM3(I,J,1)
         A = (GAMMA * ZSF1) / TV3(I,J,1)
         SLP3(I,J) = PSF1*(1+A)**COEF2
       ENDDO
       ENDDO

       PMV3=PM3
       DO J=2,KY-1
        IF(MOD(J,2).NE.0.)THEN
           DO K=1,KZ
           DO I=1,KX-1
             PMV3(I,J,K)=0.25*(PM3(I,J,K)+PM3(I+1,J,K)+
     &                         PM3(I,J-1,K)+PM3(I,J+1,K))
           END DO
           END DO
        ELSE
           DO K=1,KZ
           DO I=2,KX
             PMV3(I,J,K)=0.25*(PM3(I-1,J,K)+PM3(I,J,K)+
     &                         PM3(I,J-1,K)+PM3(I,J+1,K))
           END DO
           END DO
        END IF
        END DO

! Height at P3 grids

        DO J=1,KY
        DO I=1,KX
        DO K=1,KMX
          IF(P2(K).GE.PM3(I,J,1))THEN
            HP3(I,J,K)=ZM3(I,J,1)+
     &        TV3(I,J,1)/GAMMA*(1.-(P2(K)/PM3(I,J,1))**COEF3)
          ELSE IF(P2(K).LE.PM3(I,J,KZ))THEN
            HP3(I,J,K)=ZM3(I,J,KZ)+
     &        TV3(I,J,KZ)/GAMMA*(1.-(P2(K)/PM3(I,J,KZ))**COEF3)
          ELSE
            DO N=1,KZ-1
              IF(P2(K).LE.PM3(I,J,N).and.P2(K).GT.PM3(I,J,N+1))THEN
                HP3(I,J,K)=ZM3(I,J,N)+
     &            TV3(I,J,N)/GAMMA*(1.-(P2(K)/PM3(I,J,N))**COEF3)
                GO TO 510
              END IF
            END DO
 510         CONTINUE
          END IF
        END DO
        END DO
        END DO


       do i=1,IRX
         BLON(i)=SLON_N(kst)+i-1.
       end do
       do j=1,JRX
         BLAT(j)=SLAT_N(kst)+j-1.
       end do


      DO J=1,KY
      DO I=1,KX
        UT850(I,J)=U3(I,J,KS850)
        VT850(I,J)=V3(I,J,KS850)
      END DO
      END DO


       print*,'compute weight'
                                                                                                                 
       IGD=0
       JGD=0
       WGD=0
!       RDST=0.75
       RDST=min(0.75,2.*abs(VLON3(2,1)-VLON3(1,1)))
       RDST1=RDST*RDST
       DO J=1,JRX
       DO I=1,IRX
         NGRID=0
         NSUM(I,J)=0.
         WSUM(I,J)=0.
         DO N=1,KY
         DO K=1,KX
           COST2=(COS((VLAT3(K,N)+BLAT(J))*0.5*pi180))**2
           DIST=COST2*(VLON3(K,N)-BLON(I))**2+
     &             (VLAT3(K,N)-BLAT(J))**2
           IF(DIST.LE.RDST1)THEN
              if(ngrid .ge. maxngrid) then
                 print*,'NGRID has passed MAXNGRID.  Increase &
     &MAXNGRID in create_rel_domain.f.'
                 print*,'NGRID = ',NGRID
                 print*,'MAXNGRID = ',MAXNGRID
              end if
             NGRID=NGRID+1
             IGD(I,J,NGRID)=K
             JGD(I,J,NGRID)=N
             WGD(I,J,NGRID)=EXP(-DIST/RDST1)
             WSUM(I,J)=WSUM(I,J)+WGD(I,J,NGRID)
             NSUM(I,J)=NGRID
           END IF
         END DO
         END DO
!         print*,'I,J=',I,J,NSUM(I,J),WSUM(I,J),BLON(I),BLAT(J)
       END DO
       END DO
                       
         DO J=1,JRX
         DO I=1,IRX
           SDAT1=0.
           DO N=1,NSUM(I,J)
             I1=IGD(I,J,N)
             J1=JGD(I,J,N)
             SDAT1=SDAT1+UT850(I1,J1)*WGD(I,J,N)
           END DO
           U850(I,J)=SDAT1/WSUM(I,J)
           SDAT1=0.
           DO N=1,NSUM(I,J)
             I1=IGD(I,J,N)
             J1=JGD(I,J,N)
             SDAT1=SDAT1+VT850(I1,J1)*WGD(I,J,N)
           END DO
           V850(I,J)=SDAT1/WSUM(I,J)
         ENDDO
         ENDDO

        DO J=1,KY
        DO I=1,KX
        DO K=1,KMX
          IF(P2(K).GE.PMV3(I,J,1))THEN            ! Below PMV1(I,J,1)
!            wrk1(I,J,K)=U3(I,J,1)*(1.-(P2(K)-PMV3(I,J,1))*1.4E-5)
!            wrk2(I,J,K)=V3(I,J,1)*(1.-(P2(K)-PMV3(I,J,1))*1.4E-5)
            wrk1(I,J,K)=U3(I,J,1)
            wrk2(I,J,K)=V3(I,J,1)
          ELSE IF(P2(K).LE.PMV3(I,J,KZ))THEN
            wrk1(I,J,K)=U3(I,J,KZ)
            wrk2(I,J,K)=V3(I,J,KZ)
          ELSE
            DO N=1,KZ-1
              IF(P2(K).LE.PMV3(I,J,N).and.P2(K).GT.PMV3(I,J,N+1))THEN
                W1=ALOG(1.*PMV3(I,J,N+1))-ALOG(1.*PMV3(I,J,N))
                W=(ALOG(1.*P2(K))-ALOG(1.*PMV3(I,J,N)))/W1
                wrk1(I,J,K)=U3(I,J,N)+(U3(I,J,N+1)-U3(I,J,N))*W
                wrk2(I,J,K)=V3(I,J,N)+(V3(I,J,N+1)-V3(I,J,N))*W
                GO TO 44
              END IF
            END DO
 44         CONTINUE
          END IF
        END DO
      END DO
      END DO


       DO K=1,KMX
         K2=4*(K-1)+3+KMX+2
         K3=4*(K-1)+4+KMX+2
         K5=4*(K-1)+1+KMX+2                                   ! div (here U,V, otherwise move to next loop)
         K6=4*(K-1)+2+KMX+2                                   ! vor
         DO J=1,JRX
         DO I=1,IRX
           SDAT1=0.
           DO N=1,NSUM(I,J)
             I1=IGD(I,J,N)
             J1=JGD(I,J,N)
             SDAT1=SDAT1+wrk1(I1,J1,K)*WGD(I,J,N)
           END DO
           SDAT(I,J,K2)=SDAT1/WSUM(I,J)
           SDAT1=0.
           DO N=1,NSUM(I,J)
             I1=IGD(I,J,N)
             J1=JGD(I,J,N)
             SDAT1=SDAT1+wrk2(I1,J1,K)*WGD(I,J,N)
           END DO
           SDAT(I,J,K3)=SDAT1/WSUM(I,J)
           SDAT(I,J,K5)=SDAT(I,J,K2)
           SDAT(I,J,K6)=SDAT(I,J,K3)
         ENDDO
         ENDDO
            
!        WRITE(83)((SDAT(I,J,K2),I=1,41),J=1,41)
            
       ENDDO
            
!       DO K=1,KMX
!         K3=4*(K-1)+4+KMX+2
!         WRITE(83)((SDAT(I,J,K3),I=1,41),J=1,41)
!       END DO

!*** FOR T PS points
            
            
       IGD=0
       JGD=0
       WGD=0
       RDST=0.75
       RDST1=RDST*RDST
       DO J=1,JRX
       DO I=1,IRX
         NGRID=0
         WSUM(I,J)=0.
         DO N=1,KY
         DO K=1,KX
           COST2=(COS((HLAT3(K,N)+BLAT(J))*0.5*pi180))**2
           DIST=COST2*(HLON3(K,N)-BLON(I))**2+
     &             (HLAT3(K,N)-BLAT(J))**2
           IF(DIST.LE.RDST1)THEN
             NGRID=NGRID+1
             IGD(I,J,NGRID)=K
             JGD(I,J,NGRID)=N
             WGD(I,J,NGRID)=EXP(-DIST/RDST1)
             WSUM(I,J)=WSUM(I,J)+WGD(I,J,NGRID)
             NSUM(I,J)=NGRID
           END IF
         END DO
         END DO
       END DO
       END DO
            
       DO J=1,JRX
       DO I=1,IRX
           SDAT1=0.
           DO N=1,NSUM(I,J)
             I1=IGD(I,J,N)
             J1=JGD(I,J,N)
             SDAT1=SDAT1+Z3(I1,J1,1)*WGD(I,J,N)
           END DO
           SDAT(I,J,1)=SDAT1/WSUM(I,J)
           SDAT1=0.
           DO N=1,NSUM(I,J)
             I1=IGD(I,J,N)
             J1=JGD(I,J,N)
             SDAT1=SDAT1+SLP3(I1,J1)*WGD(I,J,N)
           END DO
           SDAT(I,J,2)=SDAT1/WSUM(I,J)
       END DO
       END DO
            
!       WRITE(83)((SDAT(I,J,1),I=1,41),J=1,41)
!       WRITE(83)((SDAT(I,J,2),I=1,41),J=1,41)
        

        K=1
        DO J=1,KY
        DO I=1,KX
          DTEMP=T3(I,J,K)-273.15
          ES=611.2*EXP(17.67*DTEMP/(DTEMP+243.5))
          QS3=0.622*ES/(PM3(I,J,K)-0.378*ES)
          RH3(I,J)=MIN(MAX(Q3(I,J,K)/QS3,0.),1.0)
        END DO
        END DO
                                                                                                                                
! Iterpolation to constant P
                                                                                                                                
      DO J=1,KY
      DO I=1,KX
        DO K=1,KMX
          IF(P2(K).GE.PM3(I,J,1))THEN            ! Below PM1(I,J,1)
            wrk1(I,J,K)=T3(I,J,1)-GAMMA*(HP3(I,J,K)-ZM3(I,J,1))
            DTEMP=T3(I,J,K)-273.15
            ES=611.2*EXP(17.67*DTEMP/(DTEMP+243.5))
            QSK=0.622*ES/(P2(K)-0.378*ES)
            wrk2(I,J,K)=RH3(I,J)*QSK            ! constant RH below KZ=1
          ELSE IF(P2(K).LE.PM3(I,J,KZ))THEN
            wrk1(I,J,K)=T3(I,J,KZ)-GAMMA*(HP3(I,J,K)-ZM3(I,J,KZ))
            wrk2(I,J,K)=Q3(I,J,KZ)             ! very small
          ELSE
            DO N=1,KZ-1
              IF(P2(K).LE.PM3(I,J,N).and.P2(K).GT.PM3(I,J,N+1))THEN
                W1=ALOG(1.*PM3(I,J,N+1))-ALOG(1.*PM3(I,J,N))
                W=(ALOG(1.*P2(K))-ALOG(1.*PM3(I,J,N)))/W1
                wrk1(I,J,K)=T3(I,J,N)+(T3(I,J,N+1)-T3(I,J,N))*W
                wrk2(I,J,K)=Q3(I,J,N)+(Q3(I,J,N+1)-Q3(I,J,N))*W
                GO TO 42
              END IF
            END DO
 42         CONTINUE
          END IF
        END DO
      END DO
      END DO
    
       DO K=1,KMX
         K1=K+2
!         K2=4*(K-1)+3+KMX+2
!         K3=4*(K-1)+4+KMX+2
         K4=K+5*KMX+2
!         K5=4*(K-1)+1+KMX+2                                   ! div
!         K6=4*(K-1)+2+KMX+2                                   ! vor
         DO J=1,JRX
         DO I=1,IRX
           SDAT1=0.
           DO N=1,NSUM(I,J)
             I1=IGD(I,J,N)
             J1=JGD(I,J,N)
             SDAT1=SDAT1+wrk1(I1,J1,K)*WGD(I,J,N)
           END DO
           SDAT(I,J,K1)=SDAT1/WSUM(I,J)
           SDAT1=0.
           DO N=1,NSUM(I,J)
             I1=IGD(I,J,N)
             J1=JGD(I,J,N)
             SDAT1=SDAT1+wrk2(I1,J1,K)*WGD(I,J,N)
           END DO
           SDAT(I,J,K4)=SDAT1/WSUM(I,J)
         ENDDO
         ENDDO
            
!        WRITE(83)((SDAT(I,J,K1),I=1,41),J=1,41)
            
       ENDDO
            
!       DO K=1,KMX
!         K4=K+5*KMX+2
!         WRITE(83)((SDAT(I,J,K4),I=1,41),J=1,41)
!       END DO
                                                                                          
       DEALLOCATE (HLON3,HLAT3,VLON3,VLAT3)
       DEALLOCATE (T3,Q3,U3,V3,Z3,P3,SLP3)
       DEALLOCATE (ZM3,PM3,PMV3,TV3,HP3,RH3)
       DEALLOCATE (wrk1,wrk2,UT850,VT850)


! finish computing SDATA

       end
