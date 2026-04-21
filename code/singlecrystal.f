C     EASYCP
C     AN EASY CRYSTAL PLASTICITY CODE FOR DEVELOPMENT
C
C     ERALP DEMIR
C     14.11.2023
C     ERALP.DEMIR@ENG.OX.AC.UK
C
C
C     CRYSTAL PLASTICITY SOLVER ROUTINE
C     THE INNER LOOP STRESS BASED SOLVER
C     REF. https://doi.org/10.1016/j.ijplas.2006.10.013
C
      SUBROUTINE CRYSTALPLASTICITY(NSLIP,DT,L,
     + NPROPS,PROPS,NSTATV,STATEV,PNEWDT,
     + SIGMAVEC,DDSIGDDEPS)
      IMPLICIT NONE
C     INPUTS
      INTEGER, INTENT(IN) :: NSLIP
      REAL(8), INTENT(IN) :: DT
      REAL(8), DIMENSION(3,3), INTENT(IN) :: L
      INTEGER, INTENT(IN) :: NPROPS
      REAL(8), DIMENSION(NPROPS), INTENT(IN) :: PROPS
      INTEGER, INTENT(IN) :: NSTATV
C     BOTH INPUT AND OUTPUTS
      REAL(8), DIMENSION(NSTATV), INTENT(INOUT) :: STATEV
      REAL(8), INTENT(INOUT) :: PNEWDT
C     OUTPUTS
      REAL(8), DIMENSION(6), INTENT(OUT) :: SIGMAVEC
      REAL(8), DIMENSION(6,6), INTENT(OUT) :: DDSIGDDEPS
C     VARIBLES USED WITHIN THIS SUBROUTINE
      INTEGER :: IITERMAX, OITERMAX
      REAL(8) :: ITOL, OTOL, CRIT, CUTBACK
      REAL(8) :: C11, C12, C44, C13, C33
      REAL(8) :: GDOT0, N, TAUC0, H0, TAUCS, A, Q
      REAL(8) :: GINV_T(3,3), SUMGAMMA_T
      REAL(8) :: GINV_0(3,3), PHI1, PHI, PHI2
      REAL(8) :: TAUC_T(NSLIP), SIGMAVEC_T(6)
      REAL(8) :: DIRC(NSLIP,3), NORC(NSLIP,3)
      REAL(8) :: DIRS(NSLIP,3), NORS(NSLIP,3)
      REAL(8) :: NS, DS, SCHMIDVEC(NSLIP,6)
      REAL(8) :: DIRS_0(NSLIP,3), NORS_0(NSLIP,3)
      REAL(8) :: NS_0, DS_0, SCHMID_0(NSLIP,3,3)
      REAL(8) :: SCHMID(NSLIP,3,3)
      REAL(8) :: SCHMIDXSCHMID(NSLIP,6,6)
      REAL(8) :: NSIJ(3,3), SNIJ(3,3)
      REAL(8) :: NSI(6), SNI(6), PMAT(6,6)
      REAL(8) :: CMATC(6,6)
      REAL(8) :: CMATS(6,6), TMAT(6,6)
      REAL(8) :: DSIGDPSI(6,6)
      REAL(8) :: DUM6X6(6,6), DUM3X3(3,3), DUM
      REAL(8) :: D(3,3), W(3,3),DET
      REAL(8) :: I3X3(3,3), I6X6(6,6)
      REAL(8) :: SIGMATRVEC(6),SIGMATR(3,3)
      REAL(8) :: SIGMA(3,3), TAU(NSLIP)
      REAL(8) :: LP(3,3), DP(3,3), DSTRANVEC(6)
      REAL(8) :: DSTRANP(3,3), DSTRANPVEC(6)
      REAL(8) :: TAUC(NSLIP), SUMGAMMA
      REAL(8) :: GDOT(NSLIP), DGDOTDTAU(NSLIP)
      REAL(8) :: LE(3,3), WE(3,3), DR(3,3)
      REAL(8) :: GINV(3,3), TAUCRES(NSLIP)
      REAL(8) :: IRESIDUE, ORESIDUE, PSI(6)
      REAL(8) :: DSIGMAVEC(6), DTAUC(NSLIP)
      REAL(8) :: TRACE
      INTEGER :: I, J, K, M, IS, CPCONV
      INTEGER :: PHASEID, OITER, IITER
C
C
C     READ DEFAULT NUMERICAL PARAMETERS
      IITERMAX=0; OITERMAX=0; ITOL=0.; OTOL=0.
      CRIT=0.; CUTBACK=0.
      CALL DEFAULTPARAMETERS(IITERMAX,OITERMAX,ITOL,OTOL,
     + CRIT,CUTBACK)
C
C     CONVERGENCE FLAG
      CPCONV=1
C     
C     IDENTITY MATRICES
      I3X3=0.
      DO I=1,3
          I3X3(I,I)=1.
      END DO
C
      I6X6=0.
      DO I=1,6
          I6X6(I,I)=1.
      END DO
C
C     INITIAL ORIENTATION
      PHI1=PROPS(1)
      PHI=PROPS(2)
      PHI2=PROPS(3)
      CALL CRYSTAL2SAMPLE(PHI1,PHI,PHI2,GINV_0)
C
C     PHASE-ID
      PHASEID = INT(PROPS(5))
C
C     ELASTIC CONSTANTS
      C11=PROPS(7)
      C12=PROPS(8)
      C44=PROPS(9)
      C13=PROPS(10)
      C33=PROPS(11)
C
C     SLIP LAW CONSTANTS
      GDOT0 = PROPS(13)
      N = PROPS(14)
C
C     HARDENING LAW CONSTANTS
      H0 = PROPS(16)
      TAUCS = PROPS(17)
      A = PROPS(18)
      Q = PROPS(19)
C
C     CONSTRUCT UNDEFORMED ELASTICITY MATRIX
      CMATC=0.
C     BCC/FCC
      IF (PHASEID.LT.3) THEN
          CMATC(1,1:3) = ( / C11, C12, C12 / )
          CMATC(2,1:3) = ( / C12, C11, C12 / )
          CMATC(3,1:3) = ( / C12, C12, C11 / )
          CMATC(4,4) = C44
          CMATC(5,5) = C44
          CMATC(6,6) = C44
C     HCP
      ELSEIF (PHASEID.EQ.3) THEN
          CMATC(1,1:3) = ( / C11, C12, C13 / )
          CMATC(2,1:3) = ( / C12, C11, C13 / )
          CMATC(3,1:3) = ( / C13, C13, C33 / )
          CMATC(4,4) = C44
          CMATC(5,5) = C44
          CMATC(6,6) = (C11-C12)/2.
      END IF     
C     
C
C     GET THE FORMER VALUE OF SIGMA FROM STATEV
      K=0
      DO I=1,6
          K=K+1
          SIGMAVEC_T(I) = STATEV(K+1)
      END DO
C
C
C      
C     GET THE FORMER VALUE OF GINV FROM STATEV
      K=0
      DO I=1,3
          DO J=1,3
              K=K+1
              GINV_T(I,J) = STATEV(K+7)
          END DO
      END DO    
C      
C     GET THE FORMER VALUE OF SUMGAMMA FROM STATEV
      SUMGAMMA_T = STATEV(17)
C      
C     GET THE FORMER VALUE OF CRSS FROM STATEV
      DO IS=1,NSLIP
          TAUC_T(IS)=STATEV(17+IS)
      END DO      
C
C     CALCULATE UNDEFORMED SLIP VECTORS (NOT STORED)
      CALL SLIPVECTORS(PHASEID,NSLIP,PROPS,NPROPS,
     + DIRC,NORC)
C
C
C
C     TRANSFORM SLIP VECTORS THE DEFORMED SAMPLE REFERENCE
C     AND NORMALIZE TO MAKE SURE THEY ARE UNIT VECTORS
      DIRS=0.; NORS=0.
      DIRS_0=0.; NORS_0=0.
      DO IS=1,NSLIP
          NS=0.; DS=0.
          NS_0=0.; DS_0=0.
          DO J=1,3
              DO K=1,3
                  DIRS(IS,J) = DIRS(IS,J) + 
     + GINV_T(J,K)* DIRC(IS,K)
                  NORS(IS,J) = NORS(IS,J) + 
     + GINV_T(J,K)* NORC(IS,K)
C
                  DIRS_0(IS,J) = DIRS_0(IS,J) + 
     + GINV_0(J,K)* DIRC(IS,K)
                  NORS_0(IS,J) = NORS_0(IS,J) + 
     + GINV_0(J,K)* NORC(IS,K)
              END DO
              DS=DS+DIRS(IS,J)**2
              NS=NS+NORS(IS,J)**2
C
              DS_0=DS_0+DIRS_0(IS,J)**2
              NS_0=NS_0+NORS_0(IS,J)**2
          END DO
C         NORMALIZE
          DO J=1,3
              DIRS(IS,J)=DIRS(IS,J)/SQRT(DS)
              NORS(IS,J)=NORS(IS,J)/SQRT(NS)
C
              DIRS_0(IS,J)=DIRS_0(IS,J)/SQRT(DS_0)
              NORS_0(IS,J)=NORS_0(IS,J)/SQRT(NS_0)
          END DO
      END DO
C
C
C     CALCULATE SCHMID TENSOR AND ITS DYADIC PRODUCT, PMAT
C     AT THE DEFORMED CONFIGURATION
      SCHMID=0.;SCHMID_0=0.
      SCHMIDXSCHMID=0.
      DO IS=1,NSLIP
C
          DO I=1,3
              DO J=1,3
                  SNIJ(I,J) = DIRS(IS,I)*NORS(IS,J)
                  NSIJ(J,I) = NORS(IS,J)*DIRS(IS,I)
                  SCHMID(IS,I,J) = SNIJ(I,J)
                  SCHMID_0(IS,I,J) = DIRS_0(IS,I)*NORS_0(IS,J)
              ENDDO
          ENDDO
C
C
C
C
C         VECTORIZE SCHMID DYADIC
C         SHEAR TERMS ARE DOUBLED
          SNI=0.
          DO I=1,3
              SNI(I)=SNIJ(I,I)
          END DO
          SNI(4)=SNIJ(1,2)+SNIJ(2,1)
          SNI(5)=SNIJ(1,3)+SNIJ(3,1)
          SNI(6)=SNIJ(2,3)+SNIJ(3,2)
C
C         VECTORIZE THE OTHER DYADIC
C         SHEAR TERMS ARE DOUBLED
          NSI=0.
          DO I=1,3
              NSI(I)=NSIJ(I,I)
          END DO
          NSI(4)=NSIJ(1,2)+NSIJ(2,1)
          NSI(5)=NSIJ(1,3)+NSIJ(3,1)
          NSI(6)=NSIJ(2,3)+NSIJ(3,2)          
C
C         VECTORIZED SCHMID TENSOR
          SCHMIDVEC(IS,1:6) = SNI
C
C         SCHMID DYADIC
          DO I=1,6
              DO J=1,6
                  SCHMIDXSCHMID(IS,I,J)=SNI(I)*NSI(J)
              ENDDO
          ENDDO
C
      ENDDO      
C
C
C     TRANSFORMATION MATRIX FOR ELASTICITY
      CALL TRANSFORM6X6(GINV_T,TMAT)
C
C     TRANSFORM ELASTICITY TO THE DEFORMED SAMPLE FRAME
C     CMATS = TMAT * CMATC * TMAT^T
      CMATS=0.
      DO I=1,6
          DO J=1,6
              DO K=1,6
                  DO M=1,6
                      CMATS(I,J)=CMATS(I,J) +
     + CMATC(K,M) * TMAT(I,K) * TMAT(J,M)
                  END DO
              END DO
          END DO
      END DO
C
C
C
C
C     TOTAL STRAIN RATE
      DO I=1,3
          DO J=1,3
              D(I,J) = (L(I,J) + L(J,I)) / 2.
          END DO
      END DO
C
C     VECTORIZED TOTAL STRAIN INCREMENT
      DO I=1,3
          DSTRANVEC(I) = D(I,I)
      END DO
      DSTRANVEC(4) = D(1,2) + D(2,1)
      DSTRANVEC(5) = D(1,3) + D(3,1)
      DSTRANVEC(6) = D(2,3) + D(3,2)
C
      DSTRANVEC=DSTRANVEC*DT
C
C     TOTAL SPIN
      DO I=1,3
          DO J=1,3
              W(I,J) = (L(I,J) - L(J,I)) / 2.
          END DO
      END DO
C
C     TRIAL STRESS
C     FIRST ADD TOTAL STRAIN INCREMENT
      SIGMATRVEC=0.
      DO I=1,6
          DO J=1,6
              SIGMATRVEC(I) = SIGMATRVEC(I) +
     + CMATS(I,J) * DSTRANVEC(J)
          END DO
      END DO
C     SECOND ADD STRESS AT THE FORMER TIME
      SIGMATRVEC = SIGMATRVEC + SIGMAVEC_T
C
C     TRACE TERM (HILL AND RICE)
      TRACE=0.
      DO I=1,3
          TRACE=TRACE+DSTRANVEC(I)
      END DO
C
C     CONVERT TO 3X3 MATRIX
      DO I=1,3
          SIGMATR(I,I)=SIGMATRVEC(I)
      END DO
      SIGMATR(1,2)=SIGMATRVEC(4); SIGMATR(2,1)=SIGMATRVEC(4)
      SIGMATR(1,3)=SIGMATRVEC(5); SIGMATR(3,1)=SIGMATRVEC(5)
      SIGMATR(3,2)=SIGMATRVEC(6); SIGMATR(2,3)=SIGMATRVEC(6)
C
C     STRESS GUESS (FULLY-PLASTIC)
      SIGMAVEC = SIGMAVEC_T
C
C     CONVERT STRESS TO A MATRIX
      DO I=1,3
          SIGMA(I,I) = SIGMAVEC(I)
      END DO
      SIGMA(1,2)=SIGMAVEC(4); SIGMA(2,1)=SIGMAVEC(4)
      SIGMA(1,3)=SIGMAVEC(5); SIGMA(3,1)=SIGMAVEC(5)
      SIGMA(3,2)=SIGMAVEC(6); SIGMA(2,3)=SIGMAVEC(6)      
!C      
!C     CO-ROTATIONAL TRIAL STRESS
!C     SIGMA = SIGMA + W * SIGMA - SIGMA * W
!      DO I=1,3
!          DO J=1,3
!              DO K=1,3
!                  SIGMA(I,J)=SIGMA(I,J)+
!     + W(I,K) * SIGMA(K,J) * DT
!              END DO
!          END DO
!      END DO
!      DO I=1,3
!          DO J=1,3
!              DO K=1,3
!                  SIGMA(I,J)=SIGMA(I,J)-
!     + SIGMA(I,K) * W(K,J) * DT
!              END DO
!          END DO
!      END DO
!C
C     VECTORIZE TRIAL STRESS      
      DO I=1,3
          SIGMATRVEC(I)=SIGMATR(I,I)
      END DO
      SIGMATRVEC(4) = SIGMATR(1,2)
      SIGMATRVEC(5) = SIGMATR(3,1)
      SIGMATRVEC(6) = SIGMATR(2,3)
C
C     VECTORIZESTRESS      
      DO I=1,3
          SIGMAVEC(I)=SIGMA(I,I)
      END DO
      SIGMAVEC(4) = SIGMA(1,2)
      SIGMAVEC(5) = SIGMA(3,1)
      SIGMAVEC(6) = SIGMA(2,3)  
C
C     ASSIGN THE INITIAL STATES
      TAUC = TAUC_T
C
C     ENTER SEMI-IMPLICIT SOLVER IF THERE IS NO ERROR UP TO THIS POINT
      IF (CPCONV.EQ.1) THEN
C
C
C         RESET THE VALUE OF THE RESIDUE
          ORESIDUE=1.D10
C
C         RESET OUTER LOOP ITERATION
          OITER=0
C         OUTER LOOP
          OUTER: DO OITER=1,OITERMAX
C
C             AVERAGE VALUE OF TAUC
              TAUC0=SUM(TAUC)/NSLIP
C
C             RESET THE VALUE OF THE RESIDUE
              IRESIDUE=1.D10              
C
C             RESET INNER LOOP ITERATION
              IITER=0
C             INNER LOOP
              INNER: DO IITER=1,IITERMAX
C
C
C                 COMPUTE RSS
                  TAU=0.
                  DO IS=1,NSLIP
                      DO I=1,6
                          TAU(IS) = TAU(IS)+
     + SCHMIDVEC(IS,I) * SIGMAVEC(I)
                      END DO
                  END DO                  
C
C
C
C                 FIND SLIP RATES
                  CALL POWERLAW(NSLIP,TAU,TAUC,N,GDOT0,
     + GDOT,DGDOTDTAU)
C
C                 PLASTIC VELOCITY GRADIENT AT THE INTERMEDIATE CONFIGURATION!!!
C                 
C                 PLASTIC STRAIN RATE AT THE DEFORMED CONFIGURATION!!!
                  LP=0.; DP=0.
                  DO I=1,3
                      DO J=1,3
                          DO IS=1,NSLIP
                              LP(I,J) = LP(I,J) +
     + GDOT(IS) * SCHMID_0(IS,I,J)
                              DP(I,J) = DP(I,J) +
     + GDOT(IS) * SCHMID(IS,I,J)
                          END DO
                      END DO
                  END DO
C
C                 PLASTIC STRAIN INCREMENT
                  DO I=1,3
                      DO J=1,3
                          DSTRANP(I,J) = (DP(I,J) + DP(J,I))/2.*DT
                      END DO
                  END DO
C                 VECTORIZE (WITH DOUBLE SHEAR COMPONENTS)
                  DO I=1,3
                      DSTRANPVEC(I)=DSTRANP(I,I)
                  END DO
                  DSTRANPVEC(4)=DSTRANP(1,2)+DSTRANP(2,1)
                  DSTRANPVEC(5)=DSTRANP(1,3)+DSTRANP(3,1)
                  DSTRANPVEC(6)=DSTRANP(3,2)+DSTRANP(2,3)
C
C
C
C                 PMAT
                  PMAT=0.
                  DO I=1,6
                      DO J=1,6
                          DO IS=1,NSLIP
                             PMAT(I,J) = PMAT(I,J) +
     + DGDOTDTAU(IS) * DT * SCHMIDXSCHMID(IS,I,J)
                          END DO
                      END DO
                  END DO
C
C
C                 PMAT IS SYMMETRIC
                  DO I=1,6
                      DO J=1,6
                          PMAT(I,J) = (PMAT(I,J) + PMAT(J,I))/2.
                      END DO
                  END DO
C
C
C                 CALCULATE TANGENT FOR NR-ITERATION
                  DUM6X6=0.
                  DO I=1,6
                      DO J=1,6
                          DO K=1,6
                              DUM6X6(I,J)=
     + DUM6X6(I,J) + CMATS(I,K) * PMAT(K,J)
                          END DO
                      END DO
                  END DO
C
C
C                 ADD IDENTITY WITH THE TRACE TERM (HILL AND RICE)
                  DUM6X6 = I6X6*(1.+TRACE) + DUM6X6
C
C
C                 TAKE THE INVERSE          
                  CALL INVERSENXN(DUM6X6,DSIGDPSI,6)
C
C                 CHECK THE INVERSION
                  IF(ANY(DSIGDPSI.NE.DSIGDPSI)) CPCONV = 0
C
C
C
C
C                 RESIDUAL
                  PSI=0.
                  DO I=1,6
                      DO J=1,6
                          PSI(I) = PSI(I) + CMATS(I,J)*DSTRANPVEC(J)
                      END DO
                  END DO
                  PSI = SIGMATRVEC - SIGMAVEC*(1.+TRACE) - PSI
C
C
                  IRESIDUE = 0.
                  DO I=1,6
                      IRESIDUE=IRESIDUE+PSI(I)**2
                  END DO
                  IRESIDUE=SQRT(IRESIDUE)
C
C                 EXIT THE INNER LOOP IF CONVERGED
                  IF (IRESIDUE<ITOL) EXIT INNER
C
C
C                 FIND STRESS INCREMENT IF NOT CONVERGED
                  DSIGMAVEC = 0.
                  DO I=1,6
                      DO J=1,6
                          DSIGMAVEC(I) = 
     + DSIGMAVEC(I) + DSIGDPSI(I,J) * PSI(J)
                      END DO
                  END DO
C
C                 CHECK THE STRESS INCREMENTS
C                 IF GREATER THAN CRITICAL CORRECT WITH LOWER VALUES
C                 REF. https://doi.org/10.1016/0022-5096(92)80003-9
                  DO I=1,6
                      IF (ABS(DSIGMAVEC(I)).GT.CRIT*TAUC0) THEN
                          DSIGMAVEC(I)=CRIT*TAUC0*SIGN(1.,DSIGMAVEC(I))
                      END IF
                  END DO
C
C                 UPDATE STRESS
                  SIGMAVEC = SIGMAVEC + DSIGMAVEC
C
C
C
C
C
C
              END DO INNER
C
C             COMPUTE THE STRAIN HARDENING 
              CALL VOCE(NSLIP,TAUC_T,GDOT,H0,TAUCS,A,Q,DT,DTAUC)
C
C
C
C             COMPUTE RESIDUAL FOR THE OUTER LEVEL LOOP
              TAUCRES = TAUC - TAUC_T - DTAUC
C
C             NORM OF THE RESIDUAL
              ORESIDUE=0.
              DO IS=1,NSLIP
                  ORESIDUE=ORESIDUE+TAUCRES(IS)**2
              END DO
              ORESIDUE= SQRT(ORESIDUE)
C
C             CONVERGENCE CHECK
              IF (ORESIDUE<OTOL) EXIT OUTER
C
C
C             UPDATE TAUC OTHERWISE
              TAUC = TAUC_T + DTAUC
C
C
          END DO OUTER
C
          IF (IITER.EQ.IITERMAX) CPCONV=0
          IF (OITER.EQ.OITERMAX) CPCONV=0
C
      END IF
C
C
C     IF CONVERGED CALCULATE JACOBIAN
      IF (CPCONV.EQ.1) THEN
C
C         DDSIGDDEPS = INV(I*(1+TR)+C*P)*(C+(SIGMA OTIMES I))
C
C         SET IT TO ELASTICITY INITIALLY
          DUM6X6=CMATS
C
C         SUBTRACT THE TERM "SIGMA OTIMES I"
C         IN DIRECT NOTATION (MATRIX FORM)
C         [S11 S22 S33 S12 S13 S23]^T * [1 1 1 0 0 0]
          DO I=1,3
              DO J=1,3
                  DUM6X6(I,J) = DUM6X6(I,J) - SIGMAVEC(I)
                  DUM6X6(I+3,J) = DUM6X6(I+3,J) - SIGMAVEC(I+3)
              END DO
          END DO
C
C         PRE-MULTIPLY WITH THE DSIGDPSI
          DDSIGDDEPS=0.
          DO I=1,6
              DO J=1,6
                  DO K=1,6
                      DDSIGDDEPS(I,J) = DDSIGDDEPS(I,J) +
     + DSIGDPSI(I,K) * DUM6X6(K,J)
                  END DO
              END DO
          END DO
C
C         MAKE JACOBIAN SYMMETRIC FOR NUMERICAL ERRORS
          DO I=1,6
              DO J=1,6
                  DDSIGDDEPS(I,J)=
     + (DDSIGDDEPS(I,J)+DDSIGDDEPS(J,I))/2.
              END DO
          END DO
C
C
      END IF   
C
C
C
C     IF CONVERGED UPDATE ORIENTATIONS
      IF (CPCONV.EQ.1) THEN
C
C         ELASTIC PART OF THE VELOCITY GRADIENT
          LE = L - LP
C
C         ELASTIC SPIN
          DO I=1,3
              DO J=1,3
                  WE(I,J) = (LE(I,J) - LE(J,I)) / 2.
              END DO
          END DO
C
C         ROTATION UPDATE
          DR = I3X3 - WE * DT
C
C         UPDATE THE ORIENTATION MATRIX
          GINV=0.
          DO I=1,3
              DO J=1,3
                  DO K=1,3
                      GINV(I,J)=
     + GINV(I,J) + DR(K,I) * GINV_T(K,J)
                  END DO
              END DO
          END DO
C
      END IF      
C
C
C     IF CONVERGED CALCULATE CAUCHY STRESS RATE
      IF (CPCONV.EQ.1) THEN
C
C
!C         CONVERT STRESS TO A MATRIX
!          DO I=1,3
!              SIGMA(I,I) = SIGMAVEC(I)
!          END DO
!          SIGMA(1,2)=SIGMAVEC(4); SIGMA(2,1)=SIGMAVEC(4)
!          SIGMA(1,3)=SIGMAVEC(5); SIGMA(3,1)=SIGMAVEC(5)
!          SIGMA(3,2)=SIGMAVEC(6); SIGMA(2,3)=SIGMAVEC(6)             
!C
!          DUM3X3=0.
!          DO I=1,3
!              DO J=1,3
!                  DO K=1,3
!                      DUM3X3(I,J)=
!     + DUM3X3(I,J) + WE(I,K) * SIGMA(K,J)
!                  END DO
!              END DO
!          END DO
!C
!          DO I=1,3
!              DO J=1,3
!                  DO K=1,3
!                      DUM3X3(I,J)=
!     + DUM3X3(I,J) - SIGMA(I,K) * WE(K,J)
!                  END DO
!              END DO
!          END DO
!C
!          SIGMA = SIGMA + DUM3X3 * DT
!C
!C
!C         VECTORIZE STRESS      
!          DO I=1,3
!              SIGMAVEC(I)=SIGMA(I,I)
!          END DO
!          SIGMAVEC(4) = SIGMA(1,2)
!          SIGMAVEC(5) = SIGMA(3,1)
!          SIGMAVEC(6) = SIGMA(2,3)             
!C
C
      END IF      
C
C
C     IF CONVERGED UPDATE THE STATES
      IF (CPCONV.EQ.1) THEN
C
C         ASSIGN CURRENT VALUE OF STRESS
          K=0
          DO I=1,6
              K=K+1
              STATEV(K+1)=SIGMAVEC(I)
          END DO
C
C
C
C      
C         ASSIGN CURRENT VALUE OF ORIENTATION
          K=0
          DO I=1,3
              DO J=1,3
                  K=K+1
                  STATEV(K+7)=GINV(I,J)
              END DO
          END DO    
C      
          DUM=0.
          DO IS=1,NSLIP
              DUM = DUM + ABS(GDOT(IS))*DT
          END DO
          SUMGAMMA = SUMGAMMA_T + DUM
C
C         ASSIGN CURRENT VALUE OF CUMULATIVE SLIP
          STATEV(17)=SUMGAMMA
C      
C         ASSIGN CURRENT VALUE OF CRSS
          DO IS=1,NSLIP
              STATEV(17+IS)=TAUC(IS)
          END DO             
          
      END IF
C
C
C     IF THERE IS NO CONVERGENCE AT THE END OF ITERATIONS
      IF (CPCONV.EQ.0) THEN
          PNEWDT = CUTBACK
          SIGMAVEC = 0.
          DDSIGDDEPS = 0.
      END IF
C
C
C
      RETURN
      END SUBROUTINE CRYSTALPLASTICITY
C
C
C
C     THIS SUBROUTINE CALCULATES TOTAL VELOCITY GRADIENT
C     DIFFERENTLY FOR SMALL STRAIN AND LARGE STRAIN CASES
      SUBROUTINE VELOCITYGRADIENT(NLGEOM,NTENS,NDI,NSHR,
     + DSTRAN,DFGRD0,DFGRD1,DT,L)
      IMPLICIT NONE
C     INPUTS
      INTEGER, INTENT(IN) :: NLGEOM
      INTEGER, INTENT(IN) :: NTENS
      INTEGER, INTENT(IN) :: NDI
      INTEGER, INTENT(IN) :: NSHR
      REAL(8), DIMENSION(NTENS), INTENT(IN) :: DSTRAN
      REAL(8), DIMENSION(3,3), INTENT(IN) :: DFGRD0
      REAL(8), DIMENSION(3,3), INTENT(IN) :: DFGRD1
      REAL(8), INTENT(IN) :: DT
      REAL(8), DIMENSION(3,3), INTENT(OUT) :: L
C     VARIABLES USED WITHIN
      REAL(8) :: I3X3(3,3)
      REAL(8) :: FDOT(3,3), FINV(3,3), DET
      INTEGER :: I, J, K
C
C     3X3 IDENTITY MATRIX
      I3X3=0.
      DO I=1,3
          I3X3(I,I)=1.
      END DO
C
      L=0.
C
C     CHECK FOR DT=0
      IF (DT.NE.0.) THEN            
C
C
C         FOR LARGE DEFORMATION CASE
          IF (NLGEOM.EQ.1) THEN
C
C
C             DEFORMATION GRADIENT TIME DERIVATIVE, FDOT
              FDOT = (DFGRD1 - DFGRD0) / DT
C
C             INVERSE OF THE DEFORMATION GRADIENT
              CALL INVERSE3X3(DFGRD1,FINV,DET)
C
C
C             VELOCITY GRADIENT
              DO I=1,3
                  DO J=1,3
                      DO K=1,3
                          L(I,J) = L(I,J) +
     + FDOT(I,K)*FINV(K,J)
                      END DO
                  END DO
              END DO
C
C
C             APPLY THE CHECK
              IF (DET.EQ.0.) L=0.
C
C
C
C
C         SMALL STRAIN CASE (SPIN=0)
          ELSE
C
C
C
C             NORMAL STRAIN TERMS TO DIAGONALS
              DO I=1,NDI
                  L(I,I)=DSTRAN(I)
              END DO
C
C             SHEAR TERMS (ABAQUS USES ENGINEERING SHEAR STRAINS)
              IF (NSHR.GT.0) THEN
C
C                 1-2 SHEAR TERM
                  IF (NSHR.GE.1) THEN
                      L(1,2)=DSTRAN(NDI+1)/2.
                      L(2,1)=DSTRAN(NDI+1)/2.
                  END IF
C
C                 1-3 SHEAR TERM
                  IF (NSHR.GE.2) THEN
                      L(1,3)=DSTRAN(NDI+2)/2.
                      L(3,1)=DSTRAN(NDI+2)/2.
                  END IF
C
C                 2-3 SHEAR TERM
                  IF (NSHR.GE.3) THEN
                      L(2,3)=DSTRAN(NDI+3)/2.
                      L(3,2)=DSTRAN(NDI+3)/2.
                  END IF               
C
C
C             DIVIDE BY TIME INCREMENT TO GET STRAIN RATES
C             ROTATION PART IS NOT COMPUTED SINCE
C             THE ROTATION INCREMENTS "DROT" IS EQUAL TO IDENTITY MATRIX FOR SMALL STRAINS
C             L = D + W = D (SMALL STRAINS)
              L = L / DT
C
C
C
              END IF
C
C
C         END OF SMALL/LARGE STRAIN CASE
          END IF
C
C         END OF DT=0 CASE
      END IF
C
C
C
C
C
      END SUBROUTINE VELOCITYGRADIENT
C
C