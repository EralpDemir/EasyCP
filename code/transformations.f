C     EASYCP
C     AN EASY CRYSTAL PLASTICITY CODE FOR DEVELOPMENT
C
C     ERALP DEMIR
C     14.11.2023
C     ERALP.DEMIR@ENG.OX.AC.UK
C
C     CALCULATE CRYSTAL2SAMPLE TRANSFORMATION
C     ONLY ENTERED ONCE      
      SUBROUTINE CRYSTAL2SAMPLE(PH1,PH,PH2,R)
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: PH1
      REAL(8), INTENT(IN) :: PH
      REAL(8), INTENT(IN) :: PH2
      REAL(8), INTENT(OUT) :: R(3,3)
C
	REAL(8) :: PHI1, PHI2, PHI
      REAL(8) :: PI
C
      PI=4.*ATAN(1.)
C     CONVERT TO RADIANS
	PHI1=PH1*PI/180.
      PHI=PH*PI/180.
	PHI2=PH2*PI/180.
C
      R=0.
      R(1,1)=(COS(PHI1)*COS(PHI2))-(SIN(PHI1)*SIN(PHI2)*COS(PHI))
      R(2,1)=-(COS(PHI1)*SIN(PHI2))-(SIN(PHI1)*COS(PHI2)*COS(PHI))
      R(3,1)=SIN(PHI1)*SIN(PHI)
      R(1,2)=(SIN(PHI1)*COS(PHI2))+(COS(PHI1)*SIN(PHI2)*COS(PHI))
      R(2,2)=-(SIN(PHI1)*SIN(PHI2))+(COS(PHI1)*COS(PHI2)*COS(PHI))
      R(3,2)=-COS(PHI1)*SIN(PHI)
      R(1,3)=SIN(PHI2)*SIN(PHI)
      R(2,3)=COS(PHI2)*SIN(PHI)
      R(3,3)=COS(PHI)
      R = TRANSPOSE(R)
C
      RETURN
      END SUBROUTINE CRYSTAL2SAMPLE
C
C
C     CALCULATE SAMPLE2CRYSTAL TRANSFORMATION
C     ONLY ENTERED ONCE      
      SUBROUTINE SAMPLE2CRYSTAL(PH1,PH,PH2,R)
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: PH1
      REAL(8), INTENT(IN) :: PH
      REAL(8), INTENT(IN) :: PH2
      REAL(8), INTENT(OUT) :: R(3,3)
C
	REAL(8) :: PHI1, PHI2, PHI
      REAL(8) :: PI
C
      PI=4.*ATAN(1.)
C     CONVERT TO RADIANS
	PHI1=PH1*PI/180.
      PHI=PH*PI/180.
	PHI2=PH2*PI/180.
C	
C
      R=0.
      R(1,1)=(COS(PHI1)*COS(PHI2))-(SIN(PHI1)*SIN(PHI2)*COS(PHI))
      R(2,1)=-(COS(PHI1)*SIN(PHI2))-(SIN(PHI1)*COS(PHI2)*COS(PHI))
      R(3,1)=SIN(PHI1)*SIN(PHI)
      R(1,2)=(SIN(PHI1)*COS(PHI2))+(COS(PHI1)*SIN(PHI2)*COS(PHI))
      R(2,2)=-(SIN(PHI1)*SIN(PHI2))+(COS(PHI1)*COS(PHI2)*COS(PHI))
      R(3,2)=-COS(PHI1)*SIN(PHI)
      R(1,3)=SIN(PHI2)*SIN(PHI)
      R(2,3)=COS(PHI2)*SIN(PHI)
      R(3,3)=COS(PHI)      
      
      RETURN
      END SUBROUTINE SAMPLE2CRYSTAL
      
      
C
C     6X6 ROTATION MATRIX
C     CONVERTS A ROTATION MATRIX TO TRANSFORM 
C     4TH ORDER TRANSFORMATION MATRIX SPECIAL FOR STRESS AND STRAIN
C     VALID FOR ROTATING SYMMETRIC TENSORS
      SUBROUTINE TRANSFORM6X6(R,T) 
      IMPLICIT NONE
C     INPUT
C     ROTATION MATRIX
      REAL(8), INTENT(IN) :: R(3,3)
C     SPECIAL TRANSFORMATION MATRIX
      REAL(8), INTENT(OUT) :: T(6,6)
C
      T(1,1) = R(1,1)*R(1,1)
      T(1,2) = R(1,2)*R(1,2)
      T(1,3) = R(1,3)*R(1,3)
      T(1,4) = 2.*R(1,1)*R(1,2)
      T(1,6) = 2.*R(1,2)*R(1,3)
      T(1,5) = 2.*R(1,3)*R(1,1)
C
      T(2,1) = R(2,1)*R(2,1)
      T(2,2) = R(2,2)*R(2,2)
      T(2,3) = R(2,3)*R(2,3)
      T(2,4) = 2.*R(2,1)*R(2,2)
      T(2,6) = 2.*R(2,2)*R(2,3)
      T(2,5) = 2.*R(2,3)*R(2,1)
C
      T(3,1) = R(3,1)*R(3,1)
      T(3,2) = R(3,2)*R(3,2)
      T(3,3) = R(3,3)*R(3,3)
      T(3,4) = 2.*R(3,1)*R(3,2)
      T(3,6) = 2.*R(3,2)*R(3,3)
      T(3,5) = 2.*R(3,3)*R(3,1)
C
      T(4,1) = R(1,1)*R(2,1)
      T(4,2) = R(1,2)*R(2,2)
      T(4,3) = R(1,3)*R(2,3)
      T(4,4) = R(1,1)*R(2,2) + R(1,2)*R(2,1)
      T(4,6) = R(1,2)*R(2,3) + R(2,2)*R(1,3) 
      T(4,5) = R(1,3)*R(2,1) + R(2,3)*R(1,1)
C
      T(6,1) = R(2,1)*R(3,1)
      T(6,2) = R(2,2)*R(3,2)
      T(6,3) = R(2,3)*R(3,3)
      T(6,4) = R(2,1)*R(3,2) + R(3,1)*R(2,2)
      T(6,6) = R(2,2)*R(3,3) + R(2,3)*R(3,2) 
      T(6,5) = R(2,3)*R(3,1) + R(3,3)*R(2,1)
C
      T(5,1) = R(3,1)*R(1,1)
      T(5,2) = R(3,2)*R(1,2)
      T(5,3) = R(3,3)*R(1,3)
      T(5,4) = R(3,1)*R(1,2) + R(1,1)*R(3,2)
      T(5,6) = R(3,2)*R(1,3) + R(1,2)*R(3,3) 
      T(5,5) = R(3,3)*R(1,1) + R(3,1)*R(1,3)
C
      RETURN
      END SUBROUTINE TRANSFORM6X6
C
C