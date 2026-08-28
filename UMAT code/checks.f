C     EASYCP
C     AN EASY CRYSTAL PLASTICITY CODE FOR DEVELOPMENT
C
C     ERALP DEMIR
C     14.11.2023
C     ERALP.DEMIR@ENG.OX.AC.UK
C
C
C     CHECK THE CONSISTENCY OF INPUTS
      SUBROUTINE CHECK(NPROPS,PROPS,NSTATV,STATEV)
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPROPS
      REAL(8), DIMENSION(NPROPS), INTENT(IN) :: PROPS
      INTEGER, INTENT(IN) :: NSTATV
      REAL(8), DIMENSION(NSTATV), INTENT(IN) :: STATEV
      INTEGER :: FLAG, SET
C
C     NPROPS  :   NUMBER OF ENTRIES IN PROPS VECTOR
C     PROPS   :   INPUT PROPS VECTOR
C     NSTATV  :   SIZE OF STATEV VECTOR
C     STATEV  :   STATE VARIABLES
C
C
C     CHECK IF THE NUMBER OF DEPVAR ENTRY IS RIGHT
      IF (NPROPS.LT.25) CALL ERROR(1)
C
C     CHECK THE PHASE ID
      IF ((PROPS(5).LE.0.).OR.(PROPS(5).GT.3.)) CALL ERROR(2)
C
C     CHECK THE NUMBER OF SLIP SYSTEMS
      FLAG=0
C     FOR BCC
      IF (INT(PROPS(5)).EQ.1) THEN
          IF (INT(PROPS(6)).EQ.12) THEN
              FLAG=1
              SET=1
          END IF
          
          IF (INT(PROPS(6)).EQ.24) THEN
              FLAG=1
              SET=2
          END IF
          
          IF (INT(PROPS(6)).EQ.48) THEN
              FLAG=1
              SET=3
          END IF

C     FOR FCC (INCLUDING CUBIC SLIP)
      ELSEIF (INT(PROPS(5)).EQ.2) THEN
          IF (INT(PROPS(6)).EQ.12) THEN 
              FLAG=1
              SET=1
          END IF
          
          IF (INT(PROPS(6)).EQ.18) THEN 
              FLAG=1
              SET=2
          END IF
C     FOR HCP
      ELSEIF (INT(PROPS(5)).EQ.3) THEN
          IF (INT(PROPS(6)).EQ.3) THEN
              FLAG=1
              SET=1
          END IF
          
          IF (INT(PROPS(6)).EQ.6) THEN
              FLAG=1
              SET=2
          END IF
          
          IF (INT(PROPS(6)).EQ.12) THEN 
              FLAG=1
              SET=3
          END IF
          
          IF (INT(PROPS(6)).EQ.24) THEN 
              FLAG=1
              SET=4
          END IF
          
          IF (INT(PROPS(6)).EQ.30) THEN 
              FLAG=1
              SET=5
          END IF
          
      END IF
      IF (FLAG.EQ.0) CALL ERROR(3)
C
C     CHECK THE ELASTIC CONSTANTS (POSITIVITY)
      IF (PROPS(7).LE.0.) CALL ERROR(4)
      IF (PROPS(8).LE.0.) CALL ERROR(4)
      IF (PROPS(9).LE.0.) CALL ERROR(4)
C     FOR HCP CASE ONLY
      IF (INT(PROPS(5)).EQ.3) THEN
          IF (PROPS(10).LE.0.) CALL ERROR(4)
          IF (PROPS(11).LE.0.) CALL ERROR(4)
      END IF
C
C     CHECK MATERIAL CONSTANTS
C     C/A RATIO
C     FOR HCP CASE ONLY
      IF (INT(PROPS(5)).EQ.3) THEN
          IF (PROPS(12).LE.0.) CALL ERROR(5)
      END IF
C     REFERENCE SLIP RATE
      IF (PROPS(13).LE.0.) CALL ERROR(5)
C     SLIP RATE SENSITIVITY EXPONENT
      IF (PROPS(14).LE.0.) CALL ERROR(5)
C
C
C     HARDENING MULTIPLIER
C     IT CAN BE SET TO ZERO IF NO HARDENING REQUIRED
      IF (PROPS(16).LT.0.) CALL ERROR(5)
C     SATURATION VALUE OF CRSS 
C     IT CAN NOT BE LESS THAN THE INITIAL VALUE
      IF (PROPS(17).LE.0.) CALL ERROR(5)
      IF (PROPS(17).LT.PROPS(15)) CALL ERROR(5)
C     HARDENING EXPONENT
C     CAN TAKE ZERO MAGNITUDE
      IF (PROPS(18).LT.0.) CALL ERROR(5)
C     LATENT HARDENING COEFFICIENT
C     CAN TAKE ZERO MAGNITUDE
      IF (PROPS(19).LT.0.) CALL ERROR(5)
C
C     INITIAL VALUE(S) OF CRSS
      FLAG=1
C     BASED ON THE NUMBER OF THE SLIP SYSTEMS
C     FOR BCC
      IF (INT(PROPS(5)).EQ.1) THEN
          IF (SET.EQ.1) THEN
              IF (PROPS(21).LE.0.) FLAG=0
          ELSEIF (SET.EQ.2) THEN
              IF (PROPS(21).LE.0.) FLAG=0
              IF (PROPS(22).LE.0.) FLAG=0
          ELSEIF (SET.EQ.3) THEN
              IF (PROPS(21).LE.0.) FLAG=0
              IF (PROPS(22).LE.0.) FLAG=0
              IF (PROPS(23).LE.0.) FLAG=0
          END IF
C     FOR FCC (INCLUDING CUBIC SLIP)
      ELSEIF (INT(PROPS(5)).EQ.2) THEN
          IF (SET.EQ.1) THEN
              IF (PROPS(21).LE.0.) FLAG=0
          ELSEIF (SET.EQ.2) THEN
              IF (PROPS(21).LE.0.) FLAG=0
              IF (PROPS(22).LE.0.) FLAG=0
          END IF
C     FOR HCP
      ELSEIF (INT(PROPS(5)).EQ.3) THEN
          IF (SET.EQ.1) THEN
              IF (PROPS(21).LE.0.) FLAG=0
          ELSEIF (SET.EQ.2) THEN
              IF (PROPS(21).LE.0.) FLAG=0
              IF (PROPS(22).LE.0.) FLAG=0
          ELSEIF (SET.EQ.3) THEN
              IF (PROPS(21).LE.0.) FLAG=0
              IF (PROPS(22).LE.0.) FLAG=0
              IF (PROPS(23).LE.0.) FLAG=0
          ELSEIF (SET.EQ.4) THEN
              IF (PROPS(21).LE.0.) FLAG=0
              IF (PROPS(22).LE.0.) FLAG=0
              IF (PROPS(23).LE.0.) FLAG=0
              IF (PROPS(24).LE.0.) FLAG=0
          ELSEIF (SET.EQ.5) THEN
              IF (PROPS(21).LE.0.) FLAG=0
              IF (PROPS(22).LE.0.) FLAG=0
              IF (PROPS(23).LE.0.) FLAG=0
              IF (PROPS(24).LE.0.) FLAG=0
              IF (PROPS(25).LE.0.) FLAG=0
          END IF
      END IF
C
C     ERROR AND EXIT IF NOT CA
      IF (FLAG.EQ.0) CALL ERROR(6)
C
C
C     CHECK IF THE NUMBER OF STATEV ENTRY IS RIGHT
      FLAG=0
C     FOR BCC
      IF (INT(PROPS(5)).EQ.1) THEN
          IF ((INT(PROPS(6)).EQ.12).AND.(INT(NSTATV).GE.29)) FLAG=1
          IF ((INT(PROPS(6)).EQ.24).AND.(INT(NSTATV).GE.41)) FLAG=1
          IF ((INT(PROPS(6)).EQ.48).AND.(INT(NSTATV).GE.65)) FLAG=1
C     FOR FCC (INCLUDING CUBIC SLIP)
      ELSEIF (INT(PROPS(5)).GE.2) THEN
          IF ((INT(PROPS(6)).EQ.12).AND.(INT(NSTATV).GE.29)) FLAG=1
          IF ((INT(PROPS(6)).EQ.18).AND.(INT(NSTATV).GE.35)) FLAG=1
C     FOR HCP
      ELSEIF (INT(PROPS(5)).GE.3) THEN
          IF ((INT(PROPS(6)).EQ.3).AND.(INT(NSTATV).GE.20)) FLAG=1
          IF ((INT(PROPS(6)).EQ.6).AND.(INT(NSTATV).GE.23)) FLAG=1
          IF ((INT(PROPS(6)).EQ.12).AND.(INT(NSTATV).GE.29)) FLAG=1
          IF ((INT(PROPS(6)).GE.24).AND.(INT(NSTATV).GE.41)) FLAG=1
          IF ((INT(PROPS(6)).GE.30).AND.(INT(NSTATV).GE.47)) FLAG=1
      END IF
      IF (FLAG.EQ.0) CALL ERROR(7)
C
C     
      RETURN
      END SUBROUTINE CHECK
C
C
C
C
C
C
C     PRINT THE INPUTS TO THE LOG AND DAT FILES
      SUBROUTINE PRINTINPUTS(NPROPS,PROPS)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPROPS
      REAL(8), DIMENSION(NPROPS), INTENT(IN) :: PROPS
C
C     WRITE TO LOG FILE
      WRITE(*,*) 'BUNGE ANGLE / PHI-1 [DEG]', PROPS(1)
      WRITE(*,*) 'BUNGE ANGLE / PHI-2 [DEG]', PROPS(2)
      WRITE(*,*) 'BUNGE ANGLE / PHI-3 [DEG]', PROPS(3)
      WRITE(*,*) 'GRAIN-ID', PROPS(4)
      WRITE(*,*) 'PHASE-ID', PROPS(5)
      WRITE(*,*) 'NUMBER OF SLIP SYSTEMS', PROPS(6)
      WRITE(*,*) 'ELASTIC CONSTANTS - C11 [MPA]', PROPS(7)
      WRITE(*,*) 'ELASTIC CONSTANTS - C12 [MPA]', PROPS(8)
      WRITE(*,*) 'ELASTIC CONSTANTS - C44 [MPA]', PROPS(9)
      WRITE(*,*) 'ELASTIC CONSTANTS (HCP) - C13 [MPA]', PROPS(10)
      WRITE(*,*) 'ELASTIC CONSTANTS (HCP) - C33 [MPA]', PROPS(11)
      WRITE(*,*) 'C/A RATIO (HCP)', PROPS(12)
      WRITE(*,*) 'REFERENCE SLIP RATE - GDOT0 [1/S]', PROPS(13)
      WRITE(*,*) 'SLIP RATE SENSITIVITY EXPONENT - N', PROPS(14)
      WRITE(*,*) 'STRAIN HARDENING - H0 [MPA]', PROPS(16)
      WRITE(*,*) 'SATURATION VALUE OF SLIP - TAUCS [MPA]', PROPS(17)
      WRITE(*,*) 'HARDENING EXPONENT - A', PROPS(18)
      WRITE(*,*) 'LATENT HARDENING COEFFICIENT - Q', PROPS(19)
      WRITE(6,*) 'INITIAL VALUES OF CRSS - TAUC0 [MPA]', PROPS(21:25)
C
C     WRITE TO DAT FILE
      WRITE(6,*) 'BUNGE ANGLE / PHI-1 [DEG]', PROPS(1)
      WRITE(6,*) 'BUNGE ANGLE / PHI-2 [DEG]', PROPS(2)
      WRITE(6,*) 'BUNGE ANGLE / PHI-3 [DEG]', PROPS(3)
      WRITE(6,*) 'GRAIN-ID', PROPS(4)
      WRITE(6,*) 'PHASE-ID', PROPS(5)
      WRITE(6,*) 'NUMBER OF SLIP SYSTEMS', PROPS(6)
      WRITE(6,*) 'ELASTIC CONSTANTS - C11 [MPA]', PROPS(7)
      WRITE(6,*) 'ELASTIC CONSTANTS - C12 [MPA]', PROPS(8)
      WRITE(6,*) 'ELASTIC CONSTANTS - C44 [MPA]', PROPS(9)
      WRITE(6,*) 'ELASTIC CONSTANTS (HCP) - C13 [MPA]', PROPS(10)
      WRITE(6,*) 'ELASTIC CONSTANTS (HCP) - C33 [MPA]', PROPS(11)
      WRITE(6,*) 'C/A RATIO (HCP)', PROPS(12)
      WRITE(6,*) 'REFERENCE SLIP RATE - GAMMADOT0 [1/S]', PROPS(13)
      WRITE(6,*) 'SLIP RATE SENSITIVITY EXPONENT - N', PROPS(14)
      WRITE(6,*) 'STRAIN HARDENING - H0 [MPA]', PROPS(16)
      WRITE(6,*) 'SATURATION VALUE OF SLIP - TAUCS [MPA]', PROPS(17)
      WRITE(6,*) 'HARDENING EXPONENT - A', PROPS(18)
      WRITE(6,*) 'LATENT HARDENING COEFFICIENT - Q', PROPS(19)
      WRITE(6,*) 'INITIAL VALUES OF CRSS - TAUC0 [MPA]', PROPS(21:25)
C
      END SUBROUTINE PRINTINPUTS
C
C
C
C
      SUBROUTINE PRINTOUTPUTS(NSLIP)
C
	IMPLICIT NONE
      INTEGER, INTENT(IN) :: NSLIP
C
      INTEGER COUNT, I
      CHARACTER*2 IJ
C
C     WRITE THE LEGEND OF THE OUTPUT VARIABLES (NSTATV)
C     OUTPUTS TO EXTRACT: 
C     1: CAUCHY STRESS (X6)
C     2: CRYSTAL ORIENTATION MATRIX (X9)
C     3: CUMULATIVE SLIP (X1)
C     4: CRSS (X NSLIP)
C
C
C
C
C
C
      COUNT = 1
      OPEN(100,FILE='../STATEV_LEGEND.TXT',ACTION='WRITE',
     + STATUS='REPLACE')
C
C
C
C
C
C     CAUCHY STRESS
C
      DO I = 1, 6
C
          COUNT = COUNT + 1
C
          SELECT CASE(I)
C
          CASE(1)
              IJ = '11'
          CASE(2)
              IJ = '22'
          CASE(3)
              IJ = '33'
          CASE(4)
              IJ = '12'
          CASE(5)
              IJ = '13'
          CASE(6)
              IJ = '23'                  
C
          END SELECT
C
C             
C
          WRITE(100,'(A7,I3,A30,A2,A6)')
     + 'STATEV-', COUNT, 
     + ':   CAUCHY STRESS COMPONENT - ',  IJ,  ' [MPA]'
C
      END DO
C
C
C
C
C
C
C
C
C     ORIENTATION MATRIX
C
      DO I = 1, 9
C
          COUNT = COUNT + 1
C
          SELECT CASE(I)
C
          CASE(1)
              IJ = '11'
          CASE(2)
              IJ = '12'
          CASE(3)
              IJ = '13'
          CASE(4)
              IJ = '21'
          CASE(5)
              IJ = '22'
          CASE(6)
              IJ = '23'
          CASE(7)
              IJ = '31'
          CASE(8)
              IJ = '32'
          CASE(9)
              IJ = '33'  
C
          END SELECT
C
C
C
          WRITE(100,'(A7,I3,A33,A2,A4)')
     + 'STATEV-', COUNT, 
     + ':   CRYSTAL ORIENTATION MATRIX - ',
     + IJ,  ' [-]'
C
      END DO
C
C
C
C
C
C
C 
C
C     CUMULATIVE SLIP
C
      COUNT = COUNT + 1
C
C
C
C
      WRITE(100,'(A7,I3,A19,A4)')
     + 'STATEV-', COUNT, 
     + ':   CUMULATIVE SLIP', ' [-]'
C
C
C
C
C     CRITICALLY-RESOLVED SHEAR STRESS (CRSS)
C
      DO I = 1, NSLIP
C
          COUNT = COUNT + 1
C
C
C
          WRITE(100,'(A7,I3,A25,I2,A6)')
     + 'STATEV-', COUNT, 
     + ':   CRSS ON SLIP SYSTEM -',  I,  ' [MPA]'
C
      END DO
C
C
C
C
C
C
      CLOSE(100)
C
C
C
C
C
C
	RETURN
      END SUBROUTINE PRINTOUTPUTS
C
C
C
