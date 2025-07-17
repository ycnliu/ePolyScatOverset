MODULE SetUp
   IMPLICIT NONE
   INTEGER, PARAMETER :: XR = SELECTED_REAL_KIND(P = 12)
END MODULE SetUp
MODULE LMPI
   USE MPI
   IMPLICIT NONE
   INTEGER :: LMPI_ierr, LMPI_myid, LMPI_numprocs
   INTEGER :: LMPI_XR, LMPI_Deriv, LMPI_SUM_Deriv
   INTEGER, PARAMETER :: LMPI_master = 0
   DOUBLE PRECISION :: LMPI_TimeStart, LMPI_TimeLast
   INTEGER :: LMPI_comm_rev, LMPI_myid_rev
   INTEGER :: LMPI_tag
   INTEGER, DIMENSION(MPI_STATUS_SIZE) :: LMPI_status
   TYPE DistType
      INTEGER :: n0, m0, n1
   END TYPE DistType
   CONTAINS
      SUBROUTINE LMPI_StartTime
         USE SetUp
         IMPLICIT NONE
         LMPI_TimeStart = MPI_WTIME()
         LMPI_TimeLast = LMPI_TimeStart
      END SUBROUTINE LMPI_StartTime
      SUBROUTINE LMPI_TimeStamp(String)
         USE SetUp
         IMPLICIT NONE
         CHARACTER (LEN = *) :: String
         DOUBLE PRECISION :: Now
         Now = MPI_WTIME()
         IF (LMPI_myid == LMPI_master) THEN
            WRITE (UNIT = 6, FMT = "('Time Now =', f15.4, '  Delta time =', f15.4, ' ', a)") &
                 & Now-LMPI_TimeStart, Now-LMPI_TimeLast, TRIM(String)
         END IF
         LMPI_TimeLast = Now
      END SUBROUTINE LMPI_TimeStamp
      SUBROUTINE LMPI_myTimeStamp(String)
         USE SetUp
         IMPLICIT NONE
         CHARACTER (LEN = *) :: String
         DOUBLE PRECISION :: Now
         Now = MPI_WTIME()
         WRITE (UNIT = 6, FMT = "('myid is', i5, '  Time Now =', f15.4, '  Delta time =', f15.4, ' ', a)") &
              & LMPI_myid, Now-LMPI_TimeStart, Now-LMPI_TimeLast, TRIM(String)
         LMPI_TimeLast = Now
      END SUBROUTINE LMPI_myTimeStamp
      SUBROUTINE LMPI_FindDist(n0, n1, m0, MTot, N)
! ** distribute MTot tasks across N processors
! ** N = n0 + n1
! ** n0 processors should do m0 of the tasks
! ** n1 processors should do m0+1 of the tasks
         USE SetUp
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: MTot, N
         INTEGER, INTENT(OUT) :: n0, n1, m0
         m0 = MTot/N
         n1 = MTot-N*m0
         n0 = N-n1
      END SUBROUTINE LMPI_FindDist
      SUBROUTINE LMPI_FindDistT(DistTot, MTot, N)
! ** distribute MTot tasks across N processors
! ** N = n0 + n1
! ** n0 processors should do m0 of the tasks
! ** n1 processors should do m0+1 of the tasks
         USE SetUp
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: MTot, N
         TYPE (DistType), INTENT(OUT) :: DistTot
         DistTot%m0 = MTot/N
         DistTot%n1 = MTot-N*DistTot%m0
         DistTot%n0 = N-DistTot%n1
      END SUBROUTINE LMPI_FindDistT
      SUBROUTINE LMPI_FindRange(mFirst, mLast, n0, m0, myid)
         USE SetUp
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: n0, m0, myid
         INTEGER, INTENT(OUT) :: mFirst, mLast

         IF (myid <= n0-1) THEN
            mFirst = myid*m0+1
            mLast = (myid+1)*m0
         ELSE
            mFirst = m0*n0+(myid-n0)*(m0+1)+1
            mLast = MFirst+m0
         END IF
 
      END SUBROUTINE LMPI_FindRange
      INTEGER FUNCTION LMPI_myFirst(DistTot, myid)
         USE SetUp
         IMPLICIT NONE
         TYPE (DistType), INTENT(IN) :: DistTot
         INTEGER, INTENT(IN) :: myid
         INTEGER :: mFirst, mLast

         CALL LMPI_FindRangeT(mFirst, mLast, DistTot, myid)
         LMPI_myFirst = mFirst
      END FUNCTION LMPI_myFirst
      INTEGER FUNCTION LMPI_myLast(DistTot, myid)
         USE SetUp
         IMPLICIT NONE
         TYPE (DistType), INTENT(IN) :: DistTot
         INTEGER, INTENT(IN) :: myid
         INTEGER :: mFirst, mLast

         CALL LMPI_FindRangeT(mFirst, mLast, DistTot, myid)
         LMPI_myLast = mLast
      END FUNCTION LMPI_myLast
      INTEGER FUNCTION LMPI_myCnt(DistTot, myid)
         USE SetUp
         IMPLICIT NONE
         TYPE (DistType), INTENT(IN) :: DistTot
         INTEGER, INTENT(IN) :: myid

         LMPI_myCnt = LMPI_myLast(DistTot, myid)-LMPI_myFirst(DistTot, myid)+1
      END FUNCTION LMPI_myCnt
      SUBROUTINE LMPI_FindRangeT(mFirst, mLast, DistTot, myid)
         USE SetUp
         IMPLICIT NONE
         TYPE (DistType), INTENT(IN) :: DistTot
         INTEGER, INTENT(IN) :: myid
         INTEGER, INTENT(OUT) :: mFirst, mLast

         IF (myid <= DistTot%n0-1) THEN
            mFirst = myid*DistTot%m0+1
            mLast = (myid+1)*DistTot%m0
         ELSE
            mFirst = DistTot%m0*DistTot%n0+(myid-DistTot%n0)*(DistTot%m0+1)+1
            mLast = MFirst+DistTot%m0
         END IF
 
      END SUBROUTINE LMPI_FindRangeT
      SUBROUTINE LMPI_MakeGathervT(mCounts, mDisplacements, DistTot)
         IMPLICIT NONE
         INTEGER, INTENT(OUT), DIMENSION(:) :: mCounts, mDisplacements
         TYPE (DistType), INTENT(IN) :: DistTot

         CALL LMPI_MakeGatherv(mCounts, mDisplacements, DistTot%n0, DistTot%n1, DistTot%m0)

      END SUBROUTINE LMPI_MakeGathervT
      SUBROUTINE LMPI_MakeGatherv(mCounts, mDisplacements, n0, n1, m0)
         IMPLICIT NONE
         INTEGER, INTENT(OUT), DIMENSION(:) :: mCounts, mDisplacements
         INTEGER, INTENT(IN) :: n0, n1, m0

         INTEGER :: iDisplacements, i

         iDisplacements = 0
         DO i = 1, n0
            mDisplacements(i) = iDisplacements
            mCounts(i) = m0
            iDisplacements = iDisplacements+m0
         END DO
         DO i = 1, n1
            mDisplacements(i+n0) = iDisplacements
            mCounts(i+n0) = m0+1
            iDisplacements = iDisplacements+m0+1
         END DO
      END SUBROUTINE LMPI_MakeGatherv
      SUBROUTINE LMPI_Type_extend(OldType, Extent, NewType, ierr)

         IMPLICIT NONE
         INTEGER :: OldType, NewType, ierr
! **         INTEGER (KIND = MPI_ADDRESS_KIND) :: Extent
         INTEGER :: Extent
        
         INTEGER, DIMENSION(2) :: BlockLen(2), Disp(2), Types(2)

         BlockLen(1) = 1
         BlockLen(2) = 1

         Disp(1) = 0
         Disp(2) = Extent

         Types(1) = OldType
         Types(2) = MPI_UB

         CALL MPI_Type_struct(2, BlockLen, Disp, Types, NewType, ierr)
         
      END SUBROUTINE LMPI_Type_extend
END MODULE LMPI
PROGRAM Testmmp
! ** Testmmp - test parallel implementation with trivial parallelization and matrix multiply

   USE SetUp
   USE LMPI
   
   IMPLICIT NONE

   REAL (KIND = XR), ALLOCATABLE, DIMENSION(:,:) :: AMat, BMat, CMat
   REAL (KIND = XR), ALLOCATABLE, DIMENSION(:) :: Results, ResultsAll
   REAL (KIND = XR) :: Cpu_time0, Cpu_time1
   INTEGER :: i, i3, i2, i1, iTask, Count1, Count0, Count_Rate
   INTEGER :: ios
   INTEGER :: TaskSizeDo, NTask, NTaskIn
   TYPE (DistType) :: DistTot
   INTEGER, PARAMETER :: TypeCalc = 4
   INTEGER :: iTypeCalc
   CHARACTER (LEN = 8), PARAMETER, DIMENSION(TypeCalc) :: &
        & TypeCalcName = ["F90     ",  "F90_Tran", "DGEMM   ", "MATMUL  "]
   
   LMPI_XR = MPI_DOUBLE_PRECISION
   CALL MPI_INIT(LMPI_ierr)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, LMPI_myid, LMPI_ierr)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, LMPI_numprocs, LMPI_ierr)

   DO
      IF (LMPI_myid == LMPI_master) THEN
         READ (UNIT = 5, FMT = *, IOSTAT = ios) TaskSizeDo, NTask, NTaskIn
      END IF
      CALL MPI_BCAST(ios, 1, MPI_INTEGER, LMPI_master, MPI_COMM_WORLD, LMPI_ierr)
      IF (ios /= 0) THEN
         EXIT
      END IF
      
      CALL MPI_BCAST(TaskSizeDo, 1, MPI_INTEGER, LMPI_master, MPI_COMM_WORLD, LMPI_ierr)
      CALL MPI_BCAST(NTask, 1, MPI_INTEGER, LMPI_master, MPI_COMM_WORLD, LMPI_ierr)
      CALL MPI_BCAST(NTaskIn, 1, MPI_INTEGER, LMPI_master, MPI_COMM_WORLD, LMPI_ierr)
      
      ALLOCATE (AMat(TaskSizeDo, TaskSizeDo), BMat(TaskSizeDo, TaskSizeDo), CMat(TaskSizeDo, TaskSizeDo))
      ALLOCATE (Results(NTask), ResultsAll(NTask*LMPI_numprocs))

      IF (LMPI_myid == LMPI_master) THEN
         WRITE (UNIT = 6, FMT = "('numprocs ', i5, ' TaskSizeDo ', i8, '  NTask ', i8, ' NTaskIn ', i8)")&
              & LMPI_numprocs, TaskSizeDo, NTask, NTaskIn
      END IF
         
      BMat = 0.0_XR
      DO i1 = 1, TaskSizeDo
         DO i2 = 1, TaskSizeDo
            BMat(i2, i1) = REAL(i2+i1, KIND = XR)/REAL(TaskSizeDo*TaskSizeDo*TaskSizeDo*NTaskIn, KIND = XR)
         END DO
         BMat(i1, i1) = 1.0_XR
      END DO
         
      DO iTypeCalc = 1, TypeCalc
         Results = 0.0_XR
         ResultsAll = 0.0_XR
         CALL MPI_BARRIER(MPI_COMM_WORLD, LMPI_ierr)
! **          CALL CPU_TIME(Cpu_time0)
! **          CALL SYSTEM_CLOCK(Count0, Count_Rate)
         CPU_time0 = MPI_WTIME()
         CALL LMPI_FindDistT(DistTot, NTaskIn, LMPI_numprocs)
      
! ** WRITE (UNIT = 6, FMT = "('myid is', i5, ' LMPI_mycnt ', i8)") LMPI_myid,  LMPI_myCnt(DistTot, LMPI_myid)
         DO iTask = 1, NTask
            AMat = 1.0_XR
! **             IF (LMPI_myid == 0) THEN
! **                AMat = 1.0_XR
! **             ELSE
! **                AMat = 0.0_XR
! **                DO i1 = 1, TaskSizeDo
! **                   AMat(i1, i1) = 1.0_XR
! **                END DO
! **             END IF
            DO i = 1, LMPI_myCnt(DistTot, LMPI_myid)
               SELECT CASE (iTypeCalc)
               CASE (1) ! Standard Fortran loop
                  CMat = 0.0_XR
                  DO i3 = 1, TaskSizeDo
                     DO i2 = 1, TaskSizeDo
                        DO i1 = 1, TaskSizeDo
                           CMat(i2, i3) = CMat(i2, i3) + AMat(i2, i1)*BMat(i1, i3)
                        END DO
                     END DO
                  END DO
                  AMat = CMat
               CASE (2) ! Fortran loop with transpose matrix
                  CMat = 0.0_XR
                  DO i3 = 1, TaskSizeDo
                     DO i2 = 1, TaskSizeDo
                        DO i1 = 1, TaskSizeDo
                           CMat(i2, i3) = CMat(i2, i3) + AMat(i1, i2)*BMat(i1, i3)
                        END DO
                     END DO
                  END DO
                  AMat = CMat
               CASE (3) ! using optimaized BLAS routine
                  CALL DGEMM('N', 'N', TaskSizeDo, TaskSizeDo, TaskSizeDo, 1.0_XR, AMat, TaskSizeDo, &
                       & BMat, TaskSizeDo, 0.0_XR, CMAt, TaskSizeDo)
                  AMat = CMat
               CASE (4) ! Using MUTMUL fortran intrinsic
                  CMat = MATMUL(AMat, BMat)
                  AMat = CMat
               END SELECT
            END DO
         
            Results(iTask) = 0.0_XR
            DO i1 = 1, TaskSizeDo
               Results(iTask) = Results(iTask)+AMat(i1, i1)
            END DO
            Results(iTask) = Results(iTask)-REAL(TaskSizeDo, KIND = XR)
         END DO

! **          WRITE (UNIT = 6, FMT = "('myid is', i5, ' Results ', 10e17.8)") LMPI_myid, Results
         CALL MPI_GATHER(Results, NTask, LMPI_XR, ResultsAll, NTask, LMPI_XR, LMPI_master, &
              & MPI_COMM_WORLD, LMPI_ierr)

         CALL MPI_BARRIER(MPI_COMM_WORLD, LMPI_ierr)
         CPU_time1 = MPI_WTIME()
! **          CALL CPU_TIME(Cpu_time1)
! **          CALL SYSTEM_CLOCK(Count1, Count_Rate)
         IF (LMPI_myid == LMPI_master) THEN
! **             WRITE (UNIT = 6, FMT = "('On master ResultsAll ', 10e17.8)") ResultsAll
            WRITE (UNIT = 6, FMT = "(a8, ' time ', e16.8, ' secs Result', e17.8)") &
              & TypeCalcName(iTypeCalc), CPU_time1-CPU_time0, SUM(ResultsAll)/REAL(NTask*LMPI_numprocs, KIND = XR)
         END IF
      END DO

      DEALLOCATE (AMat, BMat, CMat)
      DEALLOCATE (Results, ResultsAll)
   END DO
   CALL MPI_FINALIZE(LMPI_ierr)
   
END PROGRAM Testmmp
