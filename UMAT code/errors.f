C     EASYCP
C     AN EASY CRYSTAL PLASTICITY CODE FOR DEVELOPMENT
C
C     ERALP DEMIR
C     14.11.2023
C     ERALP.DEMIR@ENG.OX.AC.UK
C
C
C     WRITE ERROR MESSAGE TO THE LOG FILE (I.E. JOB-1.TXT)
      SUBROUTINE ERROR(ERRORID)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ERRORID
C
C     ERRORID  :   NUMBER OF ENTRIES IN PROPS VECTOR
      SELECT CASE(ERRORID)
C
C     NPROPS LESS THAN 20
      CASE(1)
          WRITE(*,*) '********************ERROR********************'
          WRITE(*,*) 'ERROR-001'
          WRITE(*,*) 'NUMBER OF PROPERIES IN PROPS IS LESS THAN 25!'
          WRITE(*,*) 'CHECK THE ENTRIES IN PROPS (USER-MATERIAL)!'
          WRITE(*,*) 'EXITING!'
          WRITE(*,*) '*********************************************'
          CALL XIT
C
C     PHASE-ID NOT WITHIN POSSIBLE RANGE
      CASE(2)
          WRITE(*,*) '********************ERROR********************'
          WRITE(*,*) 'ERROR-002'
          WRITE(*,*) 'THE PHASE-ID IS NOT WITHIN THE PROPER RANGE!'
          WRITE(*,*) 'CHECK THE VALUE OF PROPS(5)!'
          WRITE(*,*) 'IT MUST BE EITHER 1 OR 2 OR 3!'
          WRITE(*,*) 'EXITING!'
          WRITE(*,*) '*********************************************'
          CALL XIT
C
C     NUMBER OF SLIP SYSTEM INCONSISTENT WITH PHASE-ID
      CASE(3)
          WRITE(*,*) '********************ERROR********************'
          WRITE(*,*) 'ERROR-003'
          WRITE(*,*) 'THE NUMBER OF SLIP SYSTEMS IS INCONSISTENT!'
          WRITE(*,*) 'CHECK THE VALUE OF PROPS(6)!'
          WRITE(*,*) 'IT MUST BE CONSISTENT WITH THE PHASE-ID!'
          WRITE(*,*) 'EXITING!'
          WRITE(*,*) '*********************************************'
          CALL XIT
C
C     ELASTIC CONSTANTS (POSITIVITY)
      CASE(4)
          WRITE(*,*) '********************ERROR********************'
          WRITE(*,*) 'ERROR-004'
          WRITE(*,*) 'ELASTIC CONSTANTS CANNOT BE NEGATIVE OR ZERO!'
          WRITE(*,*) 'CHECK THE VALUES OF PROPS(7-11)!'
          WRITE(*,*) 'PROPS(10-11) COULD BE ZERO FOR HCP!'
          WRITE(*,*) 'EXITING!'
          WRITE(*,*) '*********************************************'
          CALL XIT
C
C     MATERIAL CONSTANTS (POSITIVITY)
      CASE(5)
          WRITE(*,*) '********************ERROR********************'
          WRITE(*,*) 'ERROR-005'
          WRITE(*,*) 'NEGATIVE OR ZERO MATERIAL CONSTANT(S)!'
          WRITE(*,*) 'CHECK THE VALUES OF PROPS(12-14 & 16-19)!'
          WRITE(*,*) 'PROPS(12) CAN BE ZERO FOR BCC AND FCC!'
          WRITE(*,*) 'PROPS(16) CAN BE ZERO IN CASE OF NO HARDENING!'
          WRITE(*,*) 'EXITING!'
          WRITE(*,*) '*********************************************'
          CALL XIT
C
C     CRSS - A SPECIAL CASE FOR MATERIAL CONSTANTS (POSITIVITY)
C     CRSS CAN HAVE DIFFERENT VALUE FOR DIFFERENT SLIP SYSTEM SETS
C     NUMBER OF SLIPS SYSTEM SETS FOR DIFFERENT PHASES CAN BE:
C     BCC: 1-3
C     FCC: 1-2
C     HCP: 1-5
      CASE(6)
          WRITE(*,*) '********************ERROR********************'
          WRITE(*,*) 'ERROR-006'
          WRITE(*,*) 'NEGATIVE OR ZERO CRSS FOR SLIP STSTEM SET(S)!'
          WRITE(*,*) 'CHECK THE VALUES OF PROPS(21-25)!'
          WRITE(*,*) 'SOME OF PROPS(21-25) COULD TAKE ZERO VALUE!'
          WRITE(*,*) 'BASED ON NUMSLIP THAT IDENTIFIES SLIP SYTEM SETS!'
          WRITE(*,*) 'EXITING!'
          WRITE(*,*) '*********************************************'
          CALL XIT
C
C     SIZE OF THE STATEV VECTOR IS INCONSISTENT
      CASE(7)
          WRITE(*,*) '********************ERROR********************'
          WRITE(*,*) 'ERROR-007'
          WRITE(*,*) 'THE SIZE OF STATEV IS INCONSISTENT!'
          WRITE(*,*) 'CHECK THE SIZE OF DEPVAR!'
          WRITE(*,*) 'THE DEPVAR SIZE SHALL BE >= 26+NSLIP'
          WRITE(*,*) 'IT MUST BE CONSISTENT WITH PROPS(6)!'
          WRITE(*,*) 'EXITING!'
          WRITE(*,*) '*********************************************'
          CALL XIT         
C
          
      END SELECT
C
C
C
C
C
C
C
C
C
C     
      RETURN
      END SUBROUTINE ERROR