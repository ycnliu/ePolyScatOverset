! ** 6789012  20 89012  30 89012  40 89012  50 89012  60 89012  70 89012  80 89012  90 89012 100 89012 110 89012 120 89012 130 890
MODULE RotationMat_SetUp
   IMPLICIT NONE
   INTEGER, PARAMETER :: XR = SELECTED_REAL_KIND(P = 12)
END MODULE RotationMat_SetUp
MODULE DS_FAC_MODULE
   USE RotationMat_SetUp
   IMPLICIT NONE
   
   INTEGER, PARAMETER :: NFAC1 = 150
   REAL (KIND = XR), DIMENSION(NFAC1+1) :: FAC1
   LOGICAL :: FAC1_Done
CONTAINS
   FUNCTION DS_FAC(N)
! **
! ** compute the log of N factorial
! ** 
! ** Initialize by setting FAC1_Done = .FALSE.
! **
      USE RotationMat_SetUp
      IMPLICIT NONE
      REAL (KIND = XR) :: DS_FAC
      INTEGER, INTENT(IN) :: N

      REAL (KIND = XR) :: FN
      INTEGER :: i
! **
      IF (N < 0) THEN
         DS_FAC = 0.0_XR
      ELSE
         IF (N > NFAC1) THEN
            FN=N+1
            DS_FAC=1._XR/(FN*FN)
            DS_FAC=.918938533195_XR-FN+(FN-.5_XR)*LOG(FN)+(.0833333333333_XR-DS_FAC*( &
                 .00277777777777_XR-DS_FAC*(.000793650793651_XR-.000595238095236_XR*DS_FAC)))/FN
         ELSE
            IF (.NOT. FAC1_Done) THEN
               FAC1(1)=0._XR
               DO I = 1, NFAC1
                  FAC1(I+1)=FAC1(I)+LOG(REAL(I, KIND = XR))
               END DO
            END IF
            DS_FAC = FAC1(N+1)
         END IF
      END IF
   END FUNCTION DS_FAC
END MODULE DS_FAC_MODULE
PROGRAM TestDS
   USE RotationMat_SetUp
   USE DS_FAC_MODULE
   IMPLICIT NONE

   REAL (KIND = XR), PARAMETER :: PI=3.14159265358979323846_XR
   INTEGER :: NI, J2, M12, M22, I
   REAL (KIND = XR), EXTERNAL :: DS

   FAC1_Done = .FALSE.
   
   NI = 5
   DO J2 = 1, 2
      DO M12 = -J2, J2, 2
         DO M22 = -J2, J2, 2
            WRITE (UNIT = 6, FMT = "('J =', f4.1, '  M1 =', F4.1, '  M2 =', F4.1)") J2/2.0_XR, M12/2.0_XR, M22/2.0_XR
            DO I = 1, NI
               WRITE (UNIT = 6, FMT = "('I =', i4, ' Ang deg =', F10.5, '  DS =', e16.8)") I, (I-1)&
                    & *360.0/REAL(NI-1, KIND = XR), DS(J2, M12, M22, (I-1)*PI/REAL(NI-1, KIND = XR))
            END DO
         END DO
      END DO
   END DO
END PROGRAM TestDS
FUNCTION DS(J2,MI2,MF2,BETA)
   USE RotationMat_SetUp      
   USE DS_FAC_MODULE
   IMPLICIT NONE
   REAL (KIND = XR) :: DS
   REAL (KIND = XR), INTENT(IN) :: BETA
   INTEGER, INTENT(IN) :: J2, MI2, MF2

! ** Compute the reduced rotation matrix d^J_(MI,MF)(BETA)
! ** On input J2 = 2*J, MI2 = 2*MI, MF2 = 2*MF, BETA is the angle in radians
! ** On output DS is the value of the reduced rotation matrix
! **
! ** Manuscript Title: The reduced rotation matrix.
! ** Authors: W.J. Braithwaite, J.G. Cramer
! ** Program title: DS
! ** Catalogue identifier: ABOR
! ** Journal reference: Comput. Phys. Commun. 3(1972)318
! ** Programming language: Fortran.
! ** Computer: IBM 360/91.
! ** Operating system: OS 360.
! ** RAM: 1K words
! ** Word size: 64
! ** Keywords: Atomic physics, General purpose, Nuclear physics, Rotation matrix, Rotation group, Correlation, 
! ** Euler angle, Symmetry, Helicity, Representation.
! ** Classification: 4.1.
! ** 
! ** Nature of problem:
! ** Subprogram DS is a Fortran IV double precision function which calculates the reduced matrix elements of finite 
! ** rotations in the angular momentum representation, using a standard phase convention. The four arguments of the
! ** FUNCTION are J2, twice the total angular momentum; MI2 and MF2, twice the z-projection of the total angular
! ** momentum in the inital and final coordinate systems, respectively and BETA, the Euler angle-of-rotation around y..
! ** 
! ** Solution method:
! ** A Wigner-closed-sum expression for djmmp(beta) is evaluated. Each term contains products of factorials. Using a 
! ** method similar to that of Wills (Comp. Phys. Commun. 2(1971)381), a common coefficient (containing factorials) is 
! ** evaluated by combining the logarithms of the factorials, followed by one expotentiation. The remaining expression 
! ** is written, without factorials, as a nested product. This method contracts well, in speed and accuracy, with methods
! ** that evaluate factorial products in the closed-sum coefficients term-by-term, before adding.
! ** 
! ** 

! **      SMALL D.  REDUCED ROTATION MATRIX.
! **

   INTEGER :: K1, K2, K3, K4, K5, KMAX, KMIN
   INTEGER :: K6, K7, K8, K9, I, J, K
   REAL (KIND = XR) :: ARG, S, C, P, F0
! **

   DS=0._XR

   K1=J2-MI2
   K2=J2+MI2
   K3=J2-MF2
   K4=J2+MF2
! **
   IF (MOD(K1, 2) == 0 .AND. MOD(K2, 2) == 0 .AND. MOD(K3, 2) == 0 .AND. MOD(K4, 2) == 0) THEN
      K1 = K1/2
      K2 = K2/2
      K3 = K3/2
      K4 = K4/2
      IF (MIN(K1, K2, K3, K4) >= 0) THEN
! **
         ARG = BETA/2.0_XR
         S = SIN(ARG)
         C = COS(ARG)
         K5=K2-K3
         KMAX=MIN(K1, K3)
         KMIN=MAX(0, -K5)
! **
! **      SPECIAL CASE ROTATIONS.
         IF (ABS(C) < SPACING(1.0_XR)) THEN
            KMIN = 0
            KMAX=KMIN
            C=1._XR
            IF (K5 /= 0) THEN
               RETURN
            END IF
         ELSE IF (ABS(S) < SPACING(1.0_XR)) THEN
            KMIN=K1
            KMAX=KMIN
            S=1._XR
            IF (K2-K4 /= 0) THEN
               RETURN
            END IF
         END IF
! **
         K5=K5+KMIN
         K6=K1-KMIN
         K7=K3-KMIN
         K8=K5+KMIN
         K9=K1+K2-K8
         DS=REAL((-1)**(K3-KMIN), KIND = XR)*(C**K8)*(S**K9)*EXP(.5_XR*(DS_FAC(K1)+DS_FAC(K2) &
              +DS_FAC(K3)+DS_FAC(K4))-(DS_FAC(K5)+DS_FAC(K6)+DS_FAC(K7)+DS_FAC(KMIN)))
         IF (KMAX > KMIN) THEN
            F0=(C/S)**2
            P=1._XR
            DO K=1,KMAX
               I=K-KMAX
               J=1-I
               P=1._XR-P*F0*REAL((K6+I)*(K7+I), KIND = XR)/REAL((K5+J)*(KMIN+J), KIND = XR)
            END DO
            DS=DS*P
         END IF
      END IF
   END IF
END FUNCTION DS
