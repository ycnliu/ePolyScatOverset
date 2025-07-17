PROGRAM Make2DDist
! ** 
! ** Read in a 1-D distribution for tunneling flux from a diatomic from a degenerate orbital and write out
! ** correct fluxes for the two orbitals using cos^2(gamma) and sin^2(gamma)
! **
   IMPLICIT NONE

   INTEGER, PARAMETER :: XR =8
   REAL (KIND = XR), PARAMETER :: PI = 3.141592653589793E0_XR

   INTEGER :: NThetaIn, NTheta, NPhi ! NTheta=181, NPhi=181
   INTEGER :: StrideTheta 

   REAL (KIND = XR), ALLOCATABLE, DIMENSION(:) :: Dist1DIn, ThetaValIn
   REAL (KIND = XR), ALLOCATABLE, DIMENSION(:) :: Dist1D, ThetaVal
   REAL (KIND = XR), ALLOCATABLE, DIMENSION(:) :: PhiVal
   INTEGER :: i, j

   READ (UNIT = 5, FMT = *) NThetaIn, NTheta, NPhi
   ALLOCATE (Dist1DIn(NThetaIn), ThetaValIn(NThetaIn))
   ALLOCATE (Dist1D(NTheta), ThetaVal(NTheta))
   ALLOCATE (PhiVal(NPhi))

   READ (UNIT = 5, FMT = *) (ThetaValIn(i), Dist1DIn(i), i = 1, NThetaIn)

   StrideTheta = (NThetaIn-1)/(NTheta-1)
   IF (StrideTheta*(NTheta-1)+1 /= NThetaIn) STOP 'Interpolation in Theta not programmed'

   ThetaVal = ThetaValIn(1:NThetaIn:StrideTheta)
   Dist1D = Dist1DIn(1:NThetaIn:StrideTheta)
   PhiVal = [ (REAL(i-1, KIND = XR)*360.0_XR/REAL(NPhi-1, KIND = XR), i = 1, NPhi) ]

   DO i = 1, NTheta
      DO j = 1, NPhi
         WRITE (UNIT = 6, FMT = "(3e17.8)") ThetaVal(i), PhiVal(j), Dist1D(i)*(COS(PhiVal(j)*PI/180.0_XR))**2
      END DO
      WRITE (UNIT = 6, FMT = "()")
   END DO

   DO i = 1, NTheta
      DO j = 1, NPhi
         WRITE (UNIT = 6, FMT = "(3e17.8)") ThetaVal(i), PhiVal(j), Dist1D(i)*(SIN(PhiVal(j)*PI/180.0_XR))**2
      END DO
      WRITE (UNIT = 6, FMT = "()")
   END DO
END PROGRAM Make2DDist

