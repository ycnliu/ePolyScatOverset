PROGRAM CompareG0Debug
   IMPLICIT NONE
   INTEGER, PARAMETER :: XR = 8
   CHARACTER (LEN = 256) :: File1, File2
   CHARACTER (LEN = 256) :: Label1, Label2
   CHARACTER (LEN = 256) :: Line1, Line2
   CHARACTER (LEN = 256) :: Line1PMax, Line2PMax
   CHARACTER (LEN = 256) :: Line1MMax, Line2MMax

   INTEGER, PARAMETER :: Unit1 = 10, Unit2 = 12
   INTEGER :: ios
   REAL (KIND = XR) :: VR1, VI1, VR2, VI2, DiffP, DiffM
   REAL (KIND = XR) :: NewDiffP, NewDiffM
   INTEGER :: NComp
   
   READ (UNIT = 5, FMT = *) File1, File2
   WRITE (UNIT = 6, FMT = "('File1 ', a, /, 'File2 ', a)") TRIM(File1), TRIM(File2)
   OPEN (UNIT = Unit1, FILE = File1)
   OPEN (UNIT = Unit2, FILE = File2)

   DO
      REWIND (UNIT = Unit1)
      REWIND (UNIT = Unit2)
      
      READ (UNIT = 5, FMT = *, IOSTAT = ios ) Label1, Label2
      IF (ios /= 0) THEN
         EXIT
      END IF
      
      DiffP = -1.0_XR
      DiffM = -1.0_XR
      Line1PMax = ' '
      Line1MMax = ' '
      Line2PMax = ' '
      Line2MMax = ' '
      NComp = 0
      MainLoop: DO
         DO
            READ (UNIT = Unit1, FMT = "(a)", IOSTAT = ios) Line1
            IF (ios /= 0) THEN
               EXIT MainLoop
            END IF
            IF (Line1(1:LEN_TRIM(Label1)) == TRIM(Label1)) THEN
               EXIT
            END IF
         END DO
         
         DO
            READ (UNIT = Unit2, FMT = "(a)", IOSTAT = ios) Line2
            IF (ios /= 0) THEN
               EXIT MainLoop
            END IF
            IF (Line2(1:LEN_TRIM(Label2)) == TRIM(Label2)) THEN
               EXIT
            END IF
         END DO
         
         IF (Line1(1:10) == "G0Debug EG" .OR.Line1(1:10) == "G0Debug UG" ) THEN
! ** Compare partial wave output
            READ (UNIT = Line1, FMT = "(82x, e18.8, e17.8)") VR1, VI1
            READ (UNIT = Line2, FMT = "(82x, e18.8, e17.8)") VR2, VI2
         ELSE
! ** compare coordinate representation
            READ (UNIT = Line1, FMT = "(83x, e18.8, e17.8)") VR1, VI1
            READ (UNIT = Line2, FMT = "(83x, e18.8, e17.8)") VR2, VI2
         END IF
         NComp = NComp+1
         NewDiffP = ABS(SQRT((VR1+VR2)**2+(VI1+VI2)**2))
         IF (NewDiffP > DiffP) THEN
            DiffP = NewDiffP
            Line1PMax = Line1
            Line2PMax = Line2
         END IF
         NewDiffM = ABS(SQRT((VR1-VR2)**2+(VI1-VI2)**2))
         IF (NewDiffM > DiffM) THEN
            DiffM = NewDiffM
            Line1MMax = Line1
            Line2MMax = Line2
         END IF
         
      END DO MainLoop
      IF (NComp == 0) THEN
         WRITE (UNIT = 6, FMT = "('No comparisons ', a, 1x, a)") TRIM(Label1), TRIM(Label2)
      ELSE
         IF (DiffP < DiffM) THEN
            WRITE (UNIT = 6, FMT = "(a,/,a,/'P MaxDiff', e17.8)") TRIM(Line1PMax), TRIM(Line2PMax), DiffP
         ELSE
            WRITE (UNIT = 6, FMT = "(a,/,a,/'M MaxDiff', e17.8)") TRIM(Line1MMax), TRIM(Line2MMax), DiffM
         END IF
      END IF
   END DO
END PROGRAM CompareG0Debug
