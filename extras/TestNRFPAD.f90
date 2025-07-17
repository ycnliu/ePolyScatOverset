PROGRAM TestNRFPAD
   IMPLICIT NONE

! ** Read in FLNNu and HNu functions and compute RFPAD

   INTEGER, PARAMETER :: XR = 8, InpLineL = 256
   REAL (KIND = XR), PARAMETER :: PI = 3.141592653589793E0_XR


   INTEGER, PARAMETER :: UnitFLNNu=9, UnitHNu=10

   CHARACTER (LEN = InpLineL) :: FileFLNNu, FileHNu
   CHARACTER (LEN  = InpLineL) :: LineIn
   REAL (KIND = XR) :: Chi, PhiIn

   COMPLEX (KIND = XR), DIMENSION(0:2,-1:1) :: PCoef   ! Coefficient to put in from of Legendre polynomial
                                ! PCoef(LL, Helicity)
   CHARACTER (LEN = InpLineL) :: HCoefsLabel
   INTEGER :: Helicity
   INTEGER :: NHNu, iHNu
   TYPE HNuTerm
      INTEGER :: Nu, pSin, pCos
      COMPLEX (KIND = XR) :: Coef
   END TYPE HNuTerm
   TYPE (HNuTerm), ALLOCATABLE, DIMENSION(:) :: HNu
   INTEGER :: Nu, pSin, pCos
   REAL (KIND = XR) :: FLR, FLI
   INTEGER :: ios
   INTEGER :: MaxHOrder
   INTEGER :: i
   INTEGER :: NTheIn
   REAL (KIND = XR), ALLOCATABLE, DIMENSION(:) :: TheIn, FLNNuIn
   INTEGER :: LL, N
   REAL (KIND = XR) :: Energy
   COMPLEX (KIND = XR), ALLOCATABLE, DIMENSION(:) :: TSum
   INTEGER :: iThe
   REAL (KIND = XR) :: The, Phi
   REAL (KIND = XR) :: ExpFac
   COMPLEX (KIND = XR) :: TTerm
   REAL (KIND = XR), EXTERNAL :: P
   INTEGER :: i360

   READ (UNIT = 5, FMT = *) FileFLNNu, FileHNu

   WRITE (UNIT = 6, FMT = "('File for FLNNu ', a)") TRIM(FileFLNNu)
   WRITE (UNIT = 6, fmt = "('File for HNu ', a)") TRIM(FileHNu)

   OPEN (UNIT = UnitFLNNu, FILE = TRIM(FileFLNNu))
   OPEN (UNIT = UnitHNu, FILE = TRIM(FileHNu))

   READ (UNIT = 5, FMT = *) Chi, PhiIn

   PCoef(0,:) = 1.0_XR
   PCoef(2,0) = 1.0_XR
   PCoef(2,-1) = -0.5_XR
   PCoef(2, 1) = -0.5_XR
   PCoef(1, 1) = (0.0_XR, -1.0_XR)
   PCoef(1, -1) = (0.0_XR, 1.0_XR)
   PCoef(1, 0) = 0.0_XR
   
   REWIND (UNIT = UnitHNu)
   READ (UNIT = UnitHNu, FMT = *) Helicity, HCoefsLabel
   NHNu = 0
   DO
      READ (UNIT = UnitHNu, FMT = *, IOSTAT = ios) Nu, pSin, pCos, FLR, FLI
      IF (ios /= 0 .OR. pSin < 0) THEN
         EXIT
      END IF
      NHNu = NHNu+1
   END DO

   ALLOCATE (HNu(NHNu))
   REWIND (UNIT = UnitHNu)
   READ (UNIT = UnitHNu, FMT = *) Helicity, HCoefsLabel
   WRITE (UNIT = 6, FMT = "('Helicity =', i5)") Helicity
   WRITE (UNIT = 6, FMT = "('Label - ', a)") TRIM(HCoefsLabel)
   iHNu = 0
   MaxHOrder = 0
   DO
      READ (UNIT = UnitHNu, FMT = *, IOSTAT = ios) Nu, pSin, pCos, FLR, FLI
      IF (ios /= 0 .OR. pSin < 0) THEN
         EXIT
      END IF
      iHNu = iHNu+1
      HNu(iHNu)%Nu = Nu
      HNu(iHNu)%pSin = pSin
      HNu(iHNu)%pCos = pCos
      HNu(iHNu)%Coef = CMPLX(FLR, FLI, KIND = XR)
      MaxHOrder = MAX(MaxHOrder, pSin+pCos)
   END DO

   REWIND (UNIT = UnitFLNNu)
   READ (UNIT = UnitFLNNu, FMT = "(a)") LineIn
   READ (UNIT = UnitFLNNu, FMT = *) NTheIn
   ALLOCATE (TheIn(NTheIn), FLNNuIn(NTheIn), TSum(NTheIn))
   DO i = 1, NTheIn
      READ (UNIT = UnitFLNNu, FMT = *) TheIn(i), FLNNuIn(i)
   END DO
   
   TSum = 0.0_XR
   REWIND (UNIT = UnitFLNNu)
   DO
      READ (UNIT = UnitFLNNu, FMT = "(a)", IOSTAT = ios) LineIn
      IF (LineIn(1:5) /= "FLNNu" .OR. ios /= 0) THEN
         EXIT
      END IF
      READ (UNIT = LineIn, FMT = "(10x, i2, 4x, i2, 5x, i4, 17x, f13.5)") LL, N, Nu, Energy
      READ (UNIT = UnitFLNNu, FMT = *) NTheIn
      DO i = 1, NTheIn
         READ (UNIT = UnitFLNNu, FMT = *) TheIn(i), FLNNuIn(i)
      END DO

      DO iThe = 1, NTheIn
         The = TheIn(iThe)
         Phi = PhiIn

         IF (TheIn(iThe) >= 0.0_XR) THEN
            i360 = INT(TheIn(iThe)/360._XR)
            The = TheIn(iThe)-REAL(i360*360, KIND = XR)
         ELSE
            i360 = INT(-TheIn(iThe)/360._XR)
            The = TheIn(iThe)+REAL((i360+1)*360, KIND = XR)
         END IF
         IF (The > 180._XR) THEN
            The = 360._XR-The
            Phi = Phi + 180._XR
         END IF


         DO  iHNu = 1, NHNu
! ** Correct for the fact that the H_nu were computed as cos(nu*phi) expansions so that one must divide
! ** by two for nu /= 0
            IF (HNu(iHNu)%Nu /= 0) THEN
               ExpFac = 0.5_XR
            ELSE
               ExpFac = 1.0_XR
            END IF

            IF (HNu(iHNu)%Nu == Nu) THEN
               TTerm = FLNNuIn(iThe)*PCoef(LL, Helicity)*P(LL, N, COS(Chi*PI/180.0_XR)) &
                    & *EXP((0.0_XR, 1.0_XR)*REAL(HNu(iHNu)%Nu+N, KIND = XR) &
                    & *Phi*PI/180.0_XR)
               IF (HNu(iHNu)%pSin /= 0 .AND. HNu(iHNu)%pCos /= 0) THEN
                  TTerm = TTerm*(SIN(Chi*PI/180.0_XR)**HNu(iHNu)%pSin) &
                       & *(COS(Chi*PI/180.0_XR)**HNu(iHNu)%pCos)
               ELSE IF (HNu(iHNu)%pSin /= 0) THEN
                  TTerm = TTerm*(SIN(Chi*PI/180.0_XR)**HNu(iHNu)%pSin)
               ELSE IF (HNu(iHNu)%pCos /= 0) THEN
                  TTerm = TTerm*(COS(Chi*PI/180.0_XR)**HNu(iHNu)%pCos)
               END IF
               TSum(iThe) = TSum(iThe)+TTerm*HNu(iHNu)%Coef*ExpFac
            END IF
            IF (-HNu(iHNu)%Nu == Nu .AND. HNu(iHNu)%Nu /= 0) THEN
               TTerm = FLNNuIn(iThe)*PCoef(LL, Helicity)*P(LL, N, COS(Chi*PI/180.0_XR)) &
                    & *EXP((0.0_XR, 1.0_XR)*REAL(-HNu(iHNu)%Nu+N, KIND = XR) &
                    & *Phi*PI/180.0_XR)
               IF (HNu(iHNu)%pSin /= 0 .AND. HNu(iHNu)%pCos /= 0) THEN
                  TTerm = TTerm*(SIN(Chi*PI/180.0_XR)**HNu(iHNu)%pSin) &
                       & *(COS(Chi*PI/180.0_XR)**HNu(iHNu)%pCos)
               ELSE IF (HNu(iHNu)%pSin /= 0) THEN
                  TTerm = TTerm*(SIN(Chi*PI/180.0_XR)**HNu(iHNu)%pSin)
               ELSE IF (HNu(iHNu)%pCos /= 0) THEN
                  TTerm = TTerm*(COS(Chi*PI/180.0_XR)**HNu(iHNu)%pCos)
               END IF
               TSum(iThe) = TSum(iThe)+TTerm*CONJG(HNu(iHNu)%Coef)*ExpFac
            END IF
         END DO
      END DO
   END DO

   WRITE (UNIT = 6, FMT = "('Maximum abs imag part', e17.8)") MAXVAL(ABS(AIMAG(TSum(:))))
   WRITE (UNIT = 6, FMT = "('RFPAD Chi =', f10.2, '  Phi =', f10.2)") Chi, PhiIn
   WRITE (UNIT = 6, FMT = "(i10)") NTheIn
   WRITE (UNIT = 6, FMT = "(f8.3, e17.8)") (TheIn(i), REAL(TSum(i), KIND = XR), i = 1, NTheIn)
END PROGRAM TestNRFPAD
FUNCTION P(L,M,X)
   IMPLICIT REAL (KIND = 8) (A-H, O-Z)
   INTEGER, PARAMETER :: XR = 8      
!   MADE DOUBLE PRECISION BY D LEVIN 8/17/76
!     FUNCTION TO CALCULATE ASSOCIATED LEGENDRE FUNCTIONS
!    FUNCTION P MODIFIED 12-73 TO CORRECT FORMAT
      IF(ABS(X).GT.1.0_XR) GO TO 250
      XX=X*X
      P=1.0_XR
      IF (M .NE. 0) GO TO 50
      IF (ABS(ABS(X)-1.0_XR) .LE. 0.5E-6_XR) GO TO 300
      IF(L .EQ. 0) RETURN
      IF(L .GT. 1) GO TO 10
      P=X
      RETURN
   10 CONTINUE
      P1=X
      P2=0.5_XR*(-1.0_XR+3.0_XR*X*X)
      IF (L .GT. 2) GO TO 20
      P=P2
      RETURN
   20 CONTINUE
      DO 30 I=3,L
      P3=((2*(I-1)+1)*X*P2-(I-1)*P1)/I
      P1=P2
      P2=P3
   30 CONTINUE
      P=P2
      RETURN
   50 IF (M .GT. 150 .AND. (L-M) .LT. 10) GO TO 200
      N=L-M
      IF (M .EQ. L .OR. N .EQ. 1 .OR. L .LE. 5) GO TO 100
      M1=M+1
      P1=PT(M,M,X)
      P2=PT(M1,M,X)
      M2=M+2
      DO 60 I=M2,L
!     P3=((2*(I-1)+1)*X*P2-(I-1+M)*P1)/(I-M)
      P3=(2.0_XR*(I-1.0_XR)+1.0_XR)*X*P2/(I-M)-(I-1.0_XR+M)*P1/(I-M)
      P1=P2
      P2=P3
   60 CONTINUE
      P=P2
      RETURN
  100 P=PT(L,M,X)
      RETURN
  200 CONTINUE
      WRITE(6,500) L,M
  500 FORMAT(1X,'ORDER AND DEGREE OF LEGENDRE FUNCTION ARE TOO LARGE', &
             I5, 1X, I5)
      STOP 'Problem in Legendre Polynomial'
  250 WRITE(UNIT = 6, FMT = "('The argument of the Function P  is ', e25.16, ' greater than 1.0')") X
      STOP 'Problem in Legendre 2'
 300  IF (X.GT.0.0_XR) RETURN
      IF(X.LE.0.0_XR.AND.MOD(L,2) .EQ. 0) RETURN
      IF(X.LE.0.0_XR.AND.MOD(L,2) .NE. 0) P=-1.0_XR
       RETURN
    END FUNCTION P
      FUNCTION PT(L,M,X)
         IMPLICIT REAL (KIND = 8) (A-H, O-Z)
         INTEGER, PARAMETER :: XR = 8      

      REAL (KIND = XR) :: L2

      IF (L .EQ. 0 .AND. M .EQ. 0) GO TO 90
      IF ( M.EQ. 0 .AND.(ABS(X).EQ.1.0_XR.OR.ABS(ABS(X)-1.0_XR).LE.1.0E-6_XR)) &
              GOTO 150
      IF (ABS(X).EQ.1.0_XR.OR.ABS(ABS(X)-1.0_XR).LE.1.0E-6_XR) GOTO 100
      IF (L.LT.M) GO TO 120
      NM=L-M
      XX=X*X
      P1=1.0_XR
      L2=L*2.0_XR
      IF (NM .LE. 1) GO TO 6
      NM1=NM-1
      TNM=NM
      DO 3 I=1,NM1
      P1=P1*L2/2.0_XR/TNM
      TNM=TNM-1.0_XR
    3 L2=L2-1.0_XR
      M1=M+1
      LL=M1
   6  CONTINUE
       IF (NM .LE. 1) LL=L
      DO 10 I=1,LL
      P1=P1*L2/2.0_XR
      L2=L2-1.0_XR
   10 CONTINUE
      IF (NM .LE. 1) GO TO 80
      TNM=NM
      TNM1=NM-1.0_XR
      L2=L*2-1.0_XR
      C=-TNM*TNM1/(2.0_XR*L2)
      P2=XX+C
      IF (NM .GT. 3) GO TO 30
      IF (NM .EQ. 3) P2=P2*X
      PT=P2*P1*(1.0_XR-XX)**(M/2.0_XR)
      RETURN
   30 NTERM=NM/2+1
!     NTERM = TOTAL TERM TO BE GENERATED
      LT=NTERM-2
      P2=P2*XX
      C1=C
      DO 40 I=1,LT
      TNM=TNM1-1.0_XR
      TNM1=TNM-1.0_XR
      L2=L2-2.0_XR
      C1=-C1*TNM*TNM1/(2.0_XR*(I+1.0_XR)*L2)
      P2=P2+C1
      IF (I.NE.LT) P2=P2*XX
      IF (I.EQ.LT.AND.MOD(NM,2) .NE. 0) P2=P2*X
   40 CONTINUE
      PT=P1*(1.0_XR-XX)**(M/2.0_XR)*P2
      RETURN
   80 PT=P1*(1.0_XR-XX)**(M/2.0_XR)
      IF (NM .EQ. 1) PT=PT*X
      RETURN
   90 PT=1.0_XR
      RETURN
  100 PT=0.0_XR
      RETURN
  120 WRITE(6,500) L,M
  500 FORMAT(1X, 'L IS LESS THAN M,ERROR RETURN FROM SUBROUTINE P', &
             ' L= ', I5,' M = ', I5)
      STOP 'Problem in PT'
  150  PT=1.0_XR
      IF(X .GT.0.0_XR) RETURN
      IF(X.LT.0.0_XR.AND.MOD(L,2) .EQ. 0) RETURN
      IF ( X .LT.0.0_XR .AND. MOD(L,2) .NE. 0) PT=-1.0_XR
      RETURN
 END FUNCTION PT

