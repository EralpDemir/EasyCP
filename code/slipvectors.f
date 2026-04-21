C     EASYCP
C     AN EASY CRYSTAL PLASTICITY CODE FOR DEVELOPMENT
C
C     ERALP DEMIR
C     14.11.2023
C     ERALP.DEMIR@ENG.OX.AC.UK
C
C
C      
C      
C      
C     CALCULATE UNIT SLIP VECTORS
      SUBROUTINE SLIPVECTORS(PHASEID,NSLIP,PROPS,NPROPS,
     + DIRC0,NORC0)
C     INPUTS
      INTEGER, INTENT(IN) :: PHASEID
      INTEGER, INTENT(IN) :: NSLIP
C     INPUTS
      INTEGER, INTENT(IN) :: NPROPS
      REAL(8), DIMENSION(NPROPS), INTENT(IN) :: PROPS
C     OUTPUTS
      REAL(8), DIMENSION(NSLIP,3), INTENT(OUT) :: DIRC0
      REAL(8), DIMENSION(NSLIP,3), INTENT(OUT) :: NORC0
C     VARIABLES USED WITHIN
C     BCC
      REAL(8) ::  DIR1(48,3)
      REAL(8) ::  NOR1(48,3)
C     FCC
      REAL(8) ::  DIR2(18,3)
      REAL(8) ::  NOR2(18,3)
C     HCP
      REAL(8) ::  DIR3H(30,4)
      REAL(8) ::  NOR3H(30,4)
c
      REAL(8), DIMENSION(48,3) :: DIR, NOR
      REAL(8) :: CARATIO
      REAL(8) :: VAL
      INTEGER :: I, J, IS
C
C     RESET ARRAYS
      DIRC0=0.; NORC0=0.
      DIR=0.; NOR=0.
C
C
C     BCC SLIP DIRECTIONS     
C     <111> {110} SLIP FAMILY
      DIR1(1,:)  = ( /-1.,  1.,  1. / )
	DIR1(2,:)  = ( / 1., -1.,  1. / )
	DIR1(3,:)  = ( / 1.,  1.,  1. / )
	DIR1(4,:)  = ( / 1., -1.,  1. / )
	DIR1(5,:)  = ( / 1.,  1.,  1. / )
	DIR1(6,:)  = ( / 1.,  1., -1. / )
	DIR1(7,:)  = ( / 1.,  1.,  1. / )
	DIR1(8,:)  = ( /-1.,  1.,  1. / )
	DIR1(9,:)  = ( / 1., -1.,  1. / )
	DIR1(10,:)  = ( /1.,  1., -1. / )
	DIR1(11,:) = ( /-1.,  1.,  1. / )
	DIR1(12,:) = ( / 1.,  1., -1. / )
C     <111> {112} SLIP FAMILY 
	DIR1(13,:) = ( / 1.,  1., -1. / )
	DIR1(14,:) = ( / 1., -1.,  1. / )
	DIR1(15,:) = ( /-1.,  1.,  1. / )
	DIR1(16,:) = ( / 1.,  1.,  1. / )
	DIR1(17,:) = ( / 1., -1.,  1. / )
	DIR1(18,:) = ( / 1.,  1., -1. / )
	DIR1(19,:) = ( / 1.,  1.,  1. / )
	DIR1(20,:) = ( /-1.,  1.,  1. / )
	DIR1(21,:) = ( /-1.,  1.,  1. / )
	DIR1(22,:) = ( / 1.,  1.,  1. / )
	DIR1(23,:) = ( / 1.,  1., -1. / )
	DIR1(24,:) = ( / 1., -1.,  1. / )
C
C     <111> {123} SLIP FAMILY
      DIR1(25,:) = ( / 1.,  1., -1. / )
      DIR1(26,:) = ( / 1., -1.,  1. / )
      DIR1(27,:) = ( /-1.,  1.,  1. / )
      DIR1(28,:) = ( / 1.,  1.,  1. / )
      DIR1(29,:) = ( / 1., -1.,  1. / )
      DIR1(30,:) = ( / 1.,  1., -1. / )
      DIR1(31,:) = ( / 1.,  1.,  1. / )
      DIR1(32,:) = ( /-1.,  1.,  1. / )
      DIR1(33,:) = ( / 1.,  1., -1. / )
      DIR1(34,:) = ( / 1., -1.,  1. / )
      DIR1(35,:) = ( /-1.,  1.,  1. / )
      DIR1(36,:) = ( / 1.,  1.,  1. / )
      DIR1(37,:) = ( / 1., -1.,  1. / )
      DIR1(38,:) = ( / 1.,  1., -1. / )
      DIR1(39,:) = ( / 1.,  1.,  1. / )
      DIR1(40,:) = ( /-1.,  1.,  1. / )
      DIR1(41,:) = ( /-1.,  1.,  1. / )
      DIR1(42,:) = ( / 1.,  1.,  1. / )
      DIR1(43,:) = ( / 1.,  1., -1. / )
      DIR1(44,:) = ( / 1., -1.,  1. / )
      DIR1(45,:) = ( /-1.,  1.,  1. / )
      DIR1(46,:) = ( / 1.,  1.,  1. / )
      DIR1(47,:) = ( / 1.,  1., -1. / )
      DIR1(48,:) = ( / 1., -1.,  1. / )
C
C
C     BCC SLIP PLANE NORMALS
C     <111> {110} SLIP FAMILY
      NOR1(1,:)  = ( / 1.,  1.,  0. / )
      NOR1(2,:)  = ( / 1.,  1.,  0. / )
      NOR1(3,:)  = ( /-1.,  0.,  1. / )
      NOR1(4,:)  = ( /-1.,  0.,  1. / )
      NOR1(5,:)  = ( /-1.,  1.,  0. / )
      NOR1(6,:)  = ( /-1.,  1.,  0. / )
      NOR1(7,:)  = ( / 0.,  1., -1. / )
      NOR1(8,:)  = ( / 0.,  1., -1. / )
      NOR1(9,:)  = ( / 0.,  1.,  1. / )
      NOR1(10,:) = ( / 0.,  1.,  1. / )
      NOR1(11,:) = ( / 1.,  0.,  1. / )
      NOR1(12,:) = ( / 1.,  0.,  1. / )
C
C     <111> {112} SLIP FAMILY
      NOR1(13,:) = ( / 1.,  1.,  2. / )
      NOR1(14,:) = ( /-1.,  1.,  2. / )
      NOR1(15,:) = ( / 1., -1.,  2. / )
      NOR1(16,:) = ( / 1.,  1., -2. / )
      NOR1(17,:) = ( / 1.,  2.,  1. / )
      NOR1(18,:) = ( /-1.,  2.,  1. / )
      NOR1(19,:) = ( / 1., -2.,  1. / )
      NOR1(20,:) = ( / 1.,  2., -1. / )
      NOR1(21,:) = ( / 2.,  1.,  1. / )
      NOR1(22,:) = ( /-2.,  1.,  1. / )
      NOR1(23,:) = ( / 2., -1.,  1. / )
      NOR1(24,:) = ( / 2.,  1., -1. / )
C
C     <111> {123} SLIP FAMILY
      NOR1(25,:) = ( /  1.,  2.,  3. / )
      NOR1(26,:) = ( / -1.,  2.,  3. / )
      NOR1(27,:) = ( /  1., -2.,  3. / )
      NOR1(28,:) = ( /  1.,  2., -3. / )
      NOR1(29,:) = ( /  1.,  3.,  2. / )
      NOR1(30,:) = ( / -1.,  3.,  2. / )
      NOR1(31,:) = ( /  1., -3.,  2. / )
      NOR1(32,:) = ( /  1.,  3., -2. / )
      NOR1(33,:) = ( /  2.,  1.,  3. / )
      NOR1(34,:) = ( / -2.,  1.,  3. / )
      NOR1(35,:) = ( /  2., -1.,  3. / )
      NOR1(36,:) = ( /  2.,  1., -3. / )
      NOR1(37,:) = ( /  2.,  3.,  1. / )
      NOR1(38,:) = ( / -2.,  3.,  1. / )
      NOR1(39,:) = ( /  2., -3.,  1. / )
      NOR1(40,:) = ( /  2.,  3., -1. / )
      NOR1(41,:) = ( /  3.,  1.,  2. / )
      NOR1(42,:) = ( / -3.,  1.,  2. / )
      NOR1(43,:) = ( /  3., -1.,  2. / )
      NOR1(44,:) = ( /  3.,  1., -2. / )
      NOR1(45,:) = ( /  3.,  2.,  1. / )
      NOR1(46,:) = ( / -3.,  2.,  1. / )
      NOR1(47,:) = ( /  3., -2.,  1. / )
      NOR1(48,:) = ( /  3.,  2., -1. / )      
C
C
C
C     FCC SLIP DIRECTIONS
C     <110> {111} SLIP FAMILY
	DIR2(1,:) = ( / 1., -1.,  0. / )
	DIR2(2,:) = ( / 0.,  1., -1. / )
	DIR2(3,:) = ( / 1.,  0., -1. / )
	DIR2(4,:) = ( / 1.,  1.,  0. / )
	DIR2(5,:) = ( / 0.,  1., -1. / )
	DIR2(6,:) = ( / 1.,  0.,  1. / )
	DIR2(7,:) = ( / 1.,  1.,  0. / )
	DIR2(8,:) = ( / 0.,  1.,  1. / )
	DIR2(9,:) = ( / 1.,  0., -1. / )
	DIR2(10,:) = ( / 1., -1.,  0. / )
	DIR2(11,:) = ( / 0.,  1.,  1. / )
	DIR2(12,:) = ( / 1.,  0.,  1. / )
C
C     <110> {100} CUBIC SLIP FAMILY
      DIR2(13,:) = ( / 0.,  1.,  1. / )
      DIR2(14,:) = ( / 0.,  1., -1. / )
      DIR2(15,:) = ( / 1.,  0.,  1. / )
      DIR2(16,:) = ( / 1.,  0., -1. / )
      DIR2(17,:) = ( / 1.,  1.,  0. / )
      DIR2(18,:) = ( / 1., -1.,  0. / )
C
C
C     FCC SLIP PLANE NORMALS
C     <110> {111} SLIP FAMILY
	NOR2(1,:) = ( / 1.,  1.,  1. / )
	NOR2(2,:) = ( / 1.,  1.,  1. / )
	NOR2(3,:) = ( / 1.,  1.,  1. / )
	NOR2(4,:) = ( /-1.,  1.,  1. / )
	NOR2(5,:) = ( /-1.,  1.,  1. / )
	NOR2(6,:) = ( /-1.,  1.,  1. / )
	NOR2(7,:) = ( / 1., -1.,  1. / )
	NOR2(8,:) = ( / 1., -1.,  1. / )
	NOR2(9,:) = ( / 1., -1.,  1. / )
	NOR2(10,:) = ( / 1.,  1., -1. / )
	NOR2(11,:) = ( / 1.,  1., -1. / )
	NOR2(12,:) = ( / 1.,  1., -1. / )
C
C     <110> {100} CUBIC SLIP FAMILY
      NOR2(13,:) = ( / 1.,  0.,  0. / )
      NOR2(14,:) = ( / 1.,  0.,  0. / )
      NOR2(15,:) = ( / 0.,  1.,  0. / )
      NOR2(16,:) = ( / 0.,  1.,  0. / )
      NOR2(17,:) = ( / 0.,  0.,  1. / )
      NOR2(18,:) = ( / 0.,  0.,  1. / )
C
C
C     THE DIRECTIONS ARE CORRECTED BY ALVARO 23.02.2023
C     HCP SLIP DIRECTIONS
C     <-1-1.0>{00.1} / BASAL SYSTEMS (INDEPENDENT OF C/A-RATIO)
      DIR3H(1,:)  = ( / 2., -1., -1.,  0. / )
      DIR3H(2,:)  = ( /-1.,  2., -1.,  0. / )
      DIR3H(3,:)  = ( /-1., -1.,  2.,  0. / )
C
C     <-1-1.0>{1-1.0} / PRISMATIC SYSTEMS (INDEPENDENT OF C/A-RATIO)
      DIR3H(4,:)  = ( / 2., -1., -1.,  0. / )
      DIR3H(5,:)  = ( /-1.,  2., -1.,  0. / )
      DIR3H(6,:)  = ( /-1., -1.,  2.,  0. / )
C
C     <-1-1.0>{-11.1} / 1ST ORDER PYRAMIDAL <A> SYSTEMS (DIRECTION INDEPENDENT OF C/A-RATIO)
      DIR3H(7,:) = ( /-1.,  2., -1.,  0. / )
      DIR3H(8,:) = ( /-2.,  1.,  1.,  0. / )
      DIR3H(9,:) = ( /-1., -1.,  2.,  0. / )
      DIR3H(10,:) = ( / 1., -2.,  1.,  0. / )
      DIR3H(11,:) = ( / 2., -1., -1.,  0. / )
      DIR3H(12,:) = ( / 1.,  1., -2.,  0. / )
C
C     <11.3>{-10.1} / 1ST ORDER PYRAMIDAL <C+A> SYSTEMS (DIRECTION INDEPENDENT OF C/A-RATIO)
      DIR3H(13,:) = ( /-2.,  1.,  1.,  3. / )
      DIR3H(14,:) = ( /-1., -1.,  2.,  3. / )
      DIR3H(15,:) = ( /-1., -1.,  2.,  3. / )
      DIR3H(16,:) = ( / 1., -2.,  1.,  3. / )
      DIR3H(17,:) = ( / 1., -2.,  1.,  3. / )
      DIR3H(18,:) = ( / 2., -1., -1.,  3. / )
      DIR3H(19,:) = ( / 2., -1., -1.,  3. / )
      DIR3H(20,:) = ( / 1.,  1., -2.,  3. / )
      DIR3H(21,:) = ( / 1.,  1., -2.,  3. / )
      DIR3H(22,:) = ( /-1.,  2., -1.,  3. / )
      DIR3H(23,:) = ( /-1.,  2., -1.,  3. / )
      DIR3H(24,:) = ( /-2.,  1.,  1.,  3. / )
C
C     <11.3>{-1-1.2} / 2ND ORDER PYRAMIDAL <C+A> SYSTEMS
      DIR3H(25,:)  = ( /-1., -1.,  2.,  3. / )
      DIR3H(26,:)  = ( / 1., -2.,  1.,  3. / )
      DIR3H(27,:)  = ( / 2., -1., -1.,  3. / )
      DIR3H(28,:) = ( / 1.,  1., -2.,  3. / )
      DIR3H(29,:) = ( /-1.,  2., -1.,  3. / )
      DIR3H(30,:) = ( /-2.,  1.,  1.,  3. / )
C
C      
C     HCP SLIP PLANE NORMALS
C     <-1-1.0>{00.1} / BASAL SYSTEMS (INDEPENDENT OF C/A-RATIO)
      NOR3H(1,:) = ( /0.,  0.,  0.,  1. / )
      NOR3H(2,:) = ( /0.,  0.,  0.,  1. / )
      NOR3H(3,:) = ( /0.,  0.,  0.,  1. / )
C
C     <-1-1.0>{1-1.0} / PRISMATIC SYSTEMS (INDEPENDENT OF C/A-RATIO)
      NOR3H(4,:) = ( /0.,  1., -1.,  0. / )
      NOR3H(5,:) = ( /-1.,  0.,  1.,  0. / )
      NOR3H(6,:) = ( /1., -1.,  0.,  0. / )
C
C     <-1-1.0>{-11.1} / 1ST ORDER PYRAMIDAL <A> SYSTEMS (DIRECTION INDEPENDENT OF C/A-RATIO)
      NOR3H(7,:) = ( /1.,  0., -1.,  1. / )
      NOR3H(8,:) = ( /0.,  1., -1.,  1. / )
      NOR3H(9,:) = ( /-1.,  1.,  0., 1. / )
      NOR3H(10,:) = ( /-1.,  0.,  1., 1. / )
      NOR3H(11,:) = ( /0., -1.,  1.,  1. / )
      NOR3H(12,:) = ( /1., -1.,  0.,  1. / )
C
C     <11.3>{-10.1} / 1ST ORDER PYRAMIDAL <C+A> SYSTEMS (DIRECTION INDEPENDENT OF C/A-RATIO)
      NOR3H(13,:) = ( /1.,  0., -1.,  1. / )
      NOR3H(14,:) = ( /1.,  0., -1.,  1. / )
      NOR3H(15,:) = ( /0.,  1., -1.,  1. / )
      NOR3H(16,:) = ( /0.,  1., -1.,  1. / )
      NOR3H(17,:) = ( /-1.,  1.,  0.,  1. / )
      NOR3H(18,:) = ( /-1.,  1.,  0.,  1. / )
      NOR3H(19,:) = ( /-1.,  0.,  1.,  1. / )
      NOR3H(20,:) = ( /-1.,  0.,  1.,  1. / )
      NOR3H(21,:) = ( /0., -1.,  1.,  1. / )
      NOR3H(22,:) = ( /0., -1.,  1.,  1. / )
      NOR3H(23,:) = ( /1., -1.,  0.,  1. / )
      NOR3H(24,:) = ( /1., -1.,  0.,  1. / )
C
C     <11.3>{-1-1.2} / 2ND ORDER PYRAMIDAL <C+A> SYSTEMS
      NOR3H(25,:) = ( /1.,  1., -2.,  2. / )
      NOR3H(26,:) = ( /-1.,  2., -1.,  2. / )
      NOR3H(27,:) = ( /-2.,  1.,  1.,  2. / )
      NOR3H(28,:) = ( /-1., -1.,  2.,  2. / )
      NOR3H(29,:) = ( /1., -2.,  1.,  2. / )
      NOR3H(30,:) = ( /2., -1., -1.,  2. / )      
C
C     CARATIO (C/A) RATIO IS ONLY USED FOR HCP MATERIALS
C     SO, CARATIO CAN TAKE ANY VALUE FOR OTHER PHASES THAN HCP
C
C
C     BCC PHASE
      IF (PHASEID == 1) THEN
C
C
C
C         NORMALIZE SLIP DIRECTIONS AND NORMALS
          DO IS=1, 48
C
              VAL = SQRT(DIR1(IS,1)**2.+DIR1(IS,2)**2.+DIR1(IS,3)**2.)
              DO J=1,3
                  DIR(IS,J) = DIR1(IS,J)/VAL
              END DO
C
              VAL = SQRT(NOR1(IS,1)**2.+NOR1(IS,2)**2.+NOR1(IS,3)**2.)
              DO J=1,3
                  NOR(IS,J) = NOR1(IS,J)/VAL
              END DO
C
          END DO
C
C         ASSIGN THE SLIP SYSTEM
          DIRC0(1:NSLIP,1:3) = DIR(1:NSLIP,1:3)
          NORC0(1:NSLIP,1:3) = NOR(1:NSLIP,1:3)
C
C     FCC PHASE
      ELSEIF (PHASEID == 2) THEN
C
C
C
C
C         NORMALIZE SLIP DIRECTIONS AND NORMALS
          DO IS=1,18
              VAL = SQRT(DIR2(IS,1)**2.+DIR2(IS,2)**2.+DIR2(IS,3)**2.)
              DO J=1,3
                  DIR(IS,J) = DIR2(IS,J)/VAL
              END DO
C
              VAL = SQRT(NOR2(IS,1)**2.+NOR2(IS,2)**2.+NOR2(IS,3)**2.)
              DO J=1,3
                  NOR(IS,J) = NOR2(IS,J)/VAL
              END DO
C
          END DO
C
C
C         ASSIGN THE SLIP SYSTEM
          DIRC0(1:NSLIP,1:3) = DIR(1:NSLIP,1:3)
          NORC0(1:NSLIP,1:3) = NOR(1:NSLIP,1:3)
C
C
C
C
C
C     HCP PHASE
      ELSEIF (PHASEID == 3) THEN
C
C         ASSIGN C/A RATIO (ONLY USED FOR HCP MATERIALS)
          CARATIO = PROPS(12)
C
C         SLIP DIRECTION CONVERSION
C         [UVTW]->[3U/2 (U+2V)*SQRT(3)/2 W*(C/A)])
          DO IS=1,30
              DIR(IS,1) = 3.*DIR3H(IS,1)/2.
              DIR(IS,2) = (DIR3H(IS,1) + 2.*DIR3H(IS,2))*SQRT(3.)/2.
              DIR(IS,3) = DIR3H(IS,4)*CARATIO
          END DO
C
C
C         SLIP PLANE CONVERSION
C         (HKIL)->(H (H+2K)/SQRT(3) L/(C/A))
          DO IS=1,30
              NOR(IS,1) = NOR3H(IS,1)
              NOR(IS,2) = (NOR3H(IS,1) + 2.*NOR3H(IS,2))/SQRT(3.)
              NOR(IS,3) = NOR3H(IS,4)/CARATIO
          END DO
C
C         NORMALIZE SLIP DIRECTIONS AND NORMALS
          DO IS=1,30
              VAL = SQRT(DIR(IS,1)**2.+DIR(IS,2)**2.+DIR(IS,3)**2.)
              DO J=1,3
                  DIR(IS,J) = DIR(IS,J)/VAL
              END DO
C
              VAL = SQRT(NOR(IS,1)**2.+NOR(IS,2)**2.+NOR(IS,3)**2.)
              DO J=1,3
                  NOR(IS,J) = NOR(IS,J)/VAL
              END DO
C
          END DO
C
C         ASSIGN THE SLIP SYSTEM
          DIRC0(1:NSLIP,1:3) = DIR(1:NSLIP,1:3)
          NORC0(1:NSLIP,1:3) = NOR(1:NSLIP,1:3)
C
C
C
C
C
      END IF
C
C
C
      RETURN
      END
      
      