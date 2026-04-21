C     EASYCP
C     AN EASY CRYSTAL PLASTICITY CODE FOR DEVELOPMENT
C
C     ERALP DEMIR
C     14.11.2023
C     ERALP.DEMIR@ENG.OX.AC.UK
C
C
C     INVERSION OF A 3X3 MATRIX
      SUBROUTINE INVERSE3X3(A,AINV,DET)
      IMPLICIT NONE
      REAL(8), DIMENSION(3,3), INTENT(IN) :: A
      REAL(8), DIMENSION(3,3), INTENT(OUT) :: AINV
      REAL(8), INTENT(OUT) :: DET
      
C	FIRST CALCULATE THE DETERMINANT
	CALL DETERMINANT3X3(A,DET)
C
C	IF THE DETERMINANT IS GREATER THAN CERTAIN VALUE
	IF (ABS(DET) < 1.D-20) THEN
		AINV = 0.
          DET = 0.
	ELSE
		AINV(1,1)=((A(2,2)*A(3,3))-(A(2,3)*A(3,2)))/DET
		AINV(2,1)=-((A(2,1)*A(3,3))-(A(2,3)*A(3,1)))/DET
		AINV(3,1)=((A(2,1)*A(3,2))-(A(2,2)*A(3,1)))/DET
		AINV(1,2)=-((A(1,2)*A(3,3))-(A(1,3)*A(3,2)))/DET
		AINV(2,2)=((A(1,1)*A(3,3))-(A(1,3)*A(3,1)))/DET
		AINV(3,2)=-((A(1,1)*A(3,2))-(A(1,2)*A(3,1)))/DET
		AINV(1,3)=((A(1,2)*A(2,3))-(A(1,3)*A(2,2)))/DET
		AINV(2,3)=-((A(1,1)*A(2,3))-(A(2,1)*A(1,3)))/DET
		AINV(3,3)=((A(1,1)*A(2,2))-(A(1,2)*A(2,1)))/DET
	ENDIF      

C     
      RETURN
      END SUBROUTINE INVERSE3X3
C
C
C
C
      SUBROUTINE DETERMINANT3X3(A,DET)
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: A(3,3)
      REAL(8), INTENT(OUT) :: DET
C
      DET = 0.
      DET = A(1,1)*A(2,2)*A(3,3) + 
     + A(1,2)*A(2,3)*A(3,1) + 
     + A(2,1)*A(3,2)*A(1,3) -
     + A(1,3)*A(2,2)*A(3,1) -
     + A(1,1)*A(2,3)*A(3,2) -
     + A(1,2)*A(2,1)*A(3,3)
C
      RETURN
      END SUBROUTINE DETERMINANT3X3
      
      
      
      SUBROUTINE INVERSENXN(A,AINV,N)
C     ============================================================
C     INVERSE MATRIX
C     METHOD: BASED ON DOOLITTLE LU FACTORIZATION FOR AX=B
C     ALEX G. DECEMBER 2009
C     -----------------------------------------------------------
C     INPUT ...
C     A(N,N) - ARRAY OF COEFFICIENTS FOR MATRIX A
C     N      - DIMENSION
C     OUTPUT ...
C     AINV(N,N) - INVERSE MATRIX OF A
C
C     ===========================================================
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: N
      REAL(8), INTENT(OUT) :: AINV(N,N)
      REAL(8), INTENT(IN) :: A(N,N)
      REAL(8) :: L(N,N), U(N,N), B(N), D(N), X(N)
      REAL(8) :: COEFF, C(N,N)
      INTEGER :: I, J, K
C
C     STEP 0: INITIALIZATION FOR MATRICES L AND U AND B
C     FORTRAN 90/95 ALLOWS SUCH OPERATIONS ON MATRICES
      L=0.
      U=0.
      B=0.
      C=A
C
C     STEP 1: FORWARD ELIMINATION
      DO K=1,N-1
          DO I=K+1,N
              COEFF=C(I,K)/C(K,K)
              L(I,K) = COEFF
              DO J=K+1,N
                  C(I,J) = C(I,J)-COEFF*C(K,J)
              END DO
          END DO
      END DO
C
C     STEP 2: PREPARE L AND U MATRICES 
C     L MATRIX IS A MATRIX OF THE ELIMINATION COEFFICIENT
C     + THE DIAGONAL ELEMENTS ARE 1.0
      DO I=1,N
          L(I,I) = 1.0
      END DO  
C     U MATRIX IS THE UPPER TRIANGULAR PART OF A
      DO J=1,N
          DO I=1,J
              U(I,J) = C(I,J)
          END DO
      END DO
C
C     STEP 3: COMPUTE COLUMNS OF THE INVERSE MATRIX AINV
      DO K=1,N
          B(K)=1.0
          D(1) = B(1)
C     STEP 3A: SOLVE LD=B USING THE FORWARD SUBSTITUTION
          DO I=2,N
              D(I)=B(I)
              DO J=1,I-1
                  D(I) = D(I) - L(I,J)*D(J)
              END DO
          END DO
C         STEP 3B: SOLVE UX=D USING THE BACK SUBSTITUTION
          X(N)=D(N)/U(N,N)
          DO I = N-1,1,-1
              X(I) = D(I)
              DO J=N,I+1,-1
                  X(I)=X(I)-U(I,J)*X(J)
              END DO
              X(I) = X(I)/U(I,I)
          END DO
C         STEP 3C: FILL THE SOLUTIONS X(N) INTO COLUMN K OF AINV
          DO I=1,N
              AINV(I,K) = X(I)
          END DO
          B(K)=0.0
      END DO
C
      END SUBROUTINE INVERSENXN