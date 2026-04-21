C     EASYCP
C     AN EASY CRYSTAL PLASTICITY CODE FOR DEVELOPMENT
C
C     ERALP DEMIR
C     14.11.2023
C     ERALP.DEMIR@ENG.OX.AC.UK
C
C
C     SLIP STRAIN HARDENING LAWS
C     VOCE TYPE HARDENING
      SUBROUTINE VOCE(NSLIP,TAUC,GDOT,H0,TAUCS,A,Q,DT,DTAUC)
      IMPLICIT NONE
C     INPUTS
      INTEGER, INTENT(IN) :: NSLIP
      REAL(8), DIMENSION(NSLIP), INTENT(IN) :: TAUC
      REAL(8), DIMENSION(NSLIP), INTENT(IN) :: GDOT
      REAL(8), INTENT(IN) :: H0
      REAL(8), INTENT(IN) :: TAUCS
      REAL(8), INTENT(IN) :: A
      REAL(8), INTENT(IN) :: Q
      REAL(8), INTENT(IN) :: DT
C     OUTPUT
      REAL(8), DIMENSION(NSLIP), INTENT(OUT) :: DTAUC
      REAL(8) :: SELF(NSLIP), HLTNT(NSLIP,NSLIP)
      
      INTEGER :: IS, JS, I, J, K
      
C     SELF-HARDENING 
      DO IS=1,NSLIP
          
          SELF(IS) = H0 * (1. - TAUC(IS)/TAUCS)**A * ABS(GDOT(IS)) * DT
          
      END DO
C
C     LATENT HARDENING MATRIX
      HLTNT = Q
      DO K = 1, INT(NSLIP/3.)       
	    DO I = 1, 3
              DO J = 1, 3
                  HLTNT(3*(K-1)+I, 3*(K-1)+J)=1.
              ENDDO
          ENDDO
      ENDDO
C
C
C     INCREMENT IN TAUC
      DTAUC=0.
      DO IS=1,NSLIP
          DO JS=1,NSLIP
              DTAUC(IS) = DTAUC(IS) + HLTNT(IS,JS)*SELF(JS)
          END DO
      END DO
C
      RETURN
      END SUBROUTINE VOCE
C
C
