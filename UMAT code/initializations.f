C     EASYCP
C     AN EASY CRYSTAL PLASTICITY CODE FOR DEVELOPMENT
C
C     ERALP DEMIR
C     14.11.2023
C     ERALP.DEMIR@ENG.OX.AC.UK
C
C
C     INITIALIZATION ROUTINE TO:
C     CHECK THE INPUTS
C     PRINT THE VALUES OF ALL THE INPUTS TO LOG AND DAT FILES
C     ASSIGN STATEV VECTOR INITIALLY
C     
      SUBROUTINE INITIALIZE_ONCE(NPROPS,PROPS,NSTATV,STATEV)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPROPS
      REAL(8), DIMENSION(NPROPS), INTENT(IN) :: PROPS
      INTEGER, INTENT(IN) :: NSTATV
      REAL(8), DIMENSION(NSTATV), INTENT(INOUT) :: STATEV
C
      REAL(8) :: PHI1, PHI, PHI2
      REAL(8) :: GINV(3,3)
      INTEGER :: NSLIP, I, J, K
C
C     NUMBER OF SLIP SYSTEMS
      NSLIP = INT(PROPS(6))      
C 
C     CHECK INPUT ENTRIES FOR EACH PHASE
      CALL CHECK(NPROPS,PROPS,NSTATV,STATEV)
C
C
C
C     CALCULATE CRYSTAL2SAMPLE ROTATIONS
      PHI1=PROPS(1)
      PHI=PROPS(2)
      PHI2=PROPS(3)
      CALL CRYSTAL2SAMPLE(PHI1,PHI,PHI2,GINV)
          
C     STATE INITIALIZATION
C     CAUCHY STRESS          
      STATEV(2:7)=0.
C
C     INITIAL ORIENTATION MATRIX - STATEV(8-16)
      K=0
      DO I=1,3
          DO J=1,3
              K=K+1
              STATEV(K+7) = GINV(I,J)
          END DO
      END DO
C
C     INITIAL CUMULATIVE SLIP
      STATEV(17)=0.
C
C     INITIAL VALUE OF CRSS (FROM PROPS(20-25))
C     BASED ON THE SELECTED NUMBER OF SLIP SYSTEM
C     CRSS IS ASSIGNED TO EACH SLIP SYSTEM OF A SET
C
C     FOR BCC
      IF (INT(PROPS(5)).EQ.1) THEN
C         SET-1
          IF (INT(PROPS(6)).EQ.12) THEN
              DO I=1,12
                  STATEV(17+I)=PROPS(21)
              END DO
          END IF
C         SET-2
          IF (INT(PROPS(6)).EQ.24) THEN
              DO I=1,12
                  STATEV(17+I)=PROPS(21)
              END DO
              DO I=1,12
                  STATEV(29+I)=PROPS(22)
              END DO
          END IF
C         SET-3
          IF (INT(PROPS(6)).EQ.48) THEN
              DO I=1,12
                  STATEV(17+I)=PROPS(21)
              END DO
              DO I=1,12
                  STATEV(29+I)=PROPS(22)
              END DO
              DO I=1,24
                  STATEV(41+I)=PROPS(23)
              END DO
          END IF
C
C     FOR FCC (INCLUDING CUBIC SLIP)
      ELSEIF (INT(PROPS(5)).EQ.2) THEN
C         SET-1
          IF (INT(PROPS(6)).EQ.12) THEN
              DO I=1,12
                  STATEV(17+I)=PROPS(21)
              END DO
          END IF
C         SET-2    
          IF (INT(PROPS(6)).EQ.18) THEN 
              DO I=1,12
                  STATEV(17+I)=PROPS(21)
              END DO
              DO I=1,6
                  STATEV(29+I)=PROPS(22)
              END DO
          END IF
C
C     FOR HCP
      ELSEIF (INT(PROPS(5)).EQ.3) THEN
C         SET-1
          IF (INT(PROPS(6)).EQ.3) THEN
              DO I=1,3
                  STATEV(17+I)=PROPS(21)
              END DO
          END IF
C         SET-2
          IF (INT(PROPS(6)).EQ.6) THEN
              DO I=1,3
                  STATEV(17+I)=PROPS(21)
              END DO
              DO I=1,3
                  STATEV(20+I)=PROPS(22)
              END DO              
          END IF
C         SET-3
          IF (INT(PROPS(6)).EQ.12) THEN 
              DO I=1,3
                  STATEV(17+I)=PROPS(21)
              END DO
              DO I=1,3
                  STATEV(20+I)=PROPS(22)
              END DO
              DO I=1,6
                  STATEV(23+I)=PROPS(23)
              END DO
          END IF
C         SET-4
          IF (INT(PROPS(6)).EQ.24) THEN 
              DO I=1,3
                  STATEV(17+I)=PROPS(21)
              END DO
              DO I=1,3
                  STATEV(20+I)=PROPS(22)
              END DO
              DO I=1,6
                  STATEV(23+I)=PROPS(23)
              END DO
              DO I=1,12
                  STATEV(29+I)=PROPS(24)
              END DO
          END IF
C         SET-5
          IF (INT(PROPS(6)).EQ.30) THEN 
              DO I=1,3
                  STATEV(17+I)=PROPS(21)
              END DO
              DO I=1,3
                  STATEV(20+I)=PROPS(22)
              END DO
              DO I=1,6
                  STATEV(23+I)=PROPS(23)
              END DO
              DO I=1,12
                  STATEV(29+I)=PROPS(24)
              END DO
              DO I=1,6
                  STATEV(41+I)=PROPS(25)
              END DO
          END IF
C     END OF PHASE CHECK
      END IF
C
C
C     SET THE INTIALIZATION FLAG TO UNITY
      STATEV(1)=1.          
C
C      
C
      RETURN
      END SUBROUTINE INITIALIZE_ONCE
C
C
