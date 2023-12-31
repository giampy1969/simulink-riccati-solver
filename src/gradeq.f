      SUBROUTINE GRADEQ (N,MA,A,MB,B,LOW,IGH,CPERM,WK)
C
C     *****PARAMETERS:
      INTEGER IGH,LOW,MA,MB,N
      DOUBLE PRECISION A(MA,N),B(MB,N),CPERM(N),WK(N,2)
C
C     *****LOCAL VARIABLES:
      INTEGER I,IGHM1,IM,IP1,J,JM,JP1,K
      DOUBLE PRECISION CMAX,RMAX,SUMA,SUMB,TEMP
C
C     *****FORTRAN FUNCTIONS:
      DOUBLE PRECISION DABS
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ---------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE GRADES THE SUBMATRICES OF A AND B GIVEN BY
C     STARTING INDEX LOW AND ENDING INDEX IGH IN THE GENERALIZED
C     EIGENVALUE PROBLEM A*X = (LAMBDA)*B*X BY PERMUTING ROWS AND
C     COLUMNS SUCH THAT THE NORM OF THE I-TH ROW (COLUMN) OF THE
C     A SUBMATRIX DIVIDED BY THE NORM OF THE I-TH ROW (COLUMN) OF
C     THE B SUBMATRIX BECOMES SMALLER AS I INCREASES.
C     REF.:  WARD, R. C., BALANCING THE GENERALIZED EIGENVALUE
C     PROBLEM, SIAM J. SCI. STAT. COMPUT., VOL. 2, NO. 2, JUNE 1981,
C     141-152.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C       MA,MB   INTEGER
C               ROW DIMENSIONS OF THE ARRAYS CONTAINING MATRICES
C               A AND B RESPECTIVELY, AS DECLARED IN THE MAIN CALLING
C               PROGRAM DIMENSION STATEMENT;
C
C       N       INTEGER
C               ORDER OF THE MATRICES A AND B;
C
C       A       REAL(MA,N)
C               CONTAINS THE A MATRIX OF THE GENERALIZED EIGENPROBLEM
C               DEFINED ABOVE;
C
C       B       REAL(MB,N)
C               CONTAINS THE B MATRIX OF THE GENERALIZED EIGENPROBLEM
C               DEFINED ABOVE;
C
C       LOW     INTEGER
C               SPECIFIES THE BEGINNING INDEX FOR THE ROWS AND
C               COLUMNS OF A AND B TO BE GRADED;
C
C       IGH     INTEGER
C               SPECIFIES THE ENDING INDEX FOR THE ROWS AND COLUMNS
C               OF A AND B TO BE GRADED;
C
C       WK      REAL(N,2)
C               WORK ARRAY THAT MUST CONTAIN AT LEAST 2*N LOCATIONS.
C               ONLY LOCATIONS LOW THROUGH IGH AND N+LOW THROUGH
C               N+IGH ARE REFERENCED BY THIS SUBROUTINE.
C
C     ON OUTPUT:
C
C       A,B     CONTAIN THE PERMUTED AND GRADED A AND B MATRICES;
C
C       CPERM   REAL(N)
C               CONTAINS IN ITS LOW THROUGH IGH LOCATIONS THE
C               COLUMN PERMUTATIONS APPLIED IN GRADING THE
C               SUBMATRICES.  THE OTHER LOCATIONS ARE NOT REFERENCED
C               BY THIS SUBROUTINE;
C
C       WK      CONTAINS IN ITS LOW THROUGH IGH LOCATIONS THE ROW
C               PERMUTATIONS APPLIED IN GRADING THE SUBMATRICES.
C
C     *****ALGORITHM NOTES:
C     NONE.
C
C     *****HISTORY:
C     WRITTEN BY R. C. WARD.......
C
C     ---------------------------------------------------------------
C
      IF (LOW .EQ. IGH) GO TO 510
      IGHM1 = IGH-1
C
C     COMPUTE COLUMN NORMS OF A / THOSE OF B
C
      DO 420 J = LOW,IGH
         SUMA = 0.0D0
         SUMB = 0.0D0
         DO 410 I = LOW,IGH
            SUMA = SUMA + DABS(A(I,J))
            SUMB = SUMB + DABS(B(I,J))
  410    CONTINUE
         IF (SUMB .EQ. 0.0D0) GO TO 415
         WK(J,2) = SUMA / SUMB
         GO TO 420
  415    CONTINUE
         WK(J,2) = 1.0D38
  420 CONTINUE
C
C     PERMUTE COLUMNS TO ORDER THEM BY DECREASING QUOTIENTS
C
      DO 450 J = LOW,IGHM1
         CMAX = WK(J,2)
         JM = J
         JP1 = J+1
         DO 430 K = JP1,IGH
            IF (CMAX .GE. WK(K,2)) GO TO 430
            JM = K
            CMAX = WK(K,2)
  430    CONTINUE
         CPERM(J) = JM
         IF (JM .EQ. J) GO TO 450
         TEMP = WK(J,2)
         WK(J,2) = WK(JM,2)
         WK(JM,2) = TEMP
         DO 440 I = 1,IGH
            TEMP = B(I,J)
            B(I,J) = B(I,JM)
            B(I,JM) = TEMP
            TEMP = A(I,J)
            A(I,J) = A(I,JM)
            A(I,JM) = TEMP
  440    CONTINUE
  450 CONTINUE
      CPERM(IGH) = IGH
C
C     COMPUTE ROW NORMS OF A / THOSE OF B
C
      DO 470 I = LOW,IGH
         SUMA = 0.0D0
         SUMB = 0.0D0
         DO 460 J = LOW,IGH
            SUMA = SUMA + DABS(A(I,J))
            SUMB = SUMB + DABS(B(I,J))
  460    CONTINUE
         IF (SUMB .EQ. 0.0D0) GO TO 465
         WK(I,2) = SUMA / SUMB
         GO TO 470
  465    CONTINUE
         WK(I,2) = 1.0D38
C
C     PERMUTE ROWS TO ORDER THEM BY DECREASING QUOTIENTS
C
  470 CONTINUE
      DO 500 I = LOW,IGHM1
         RMAX = WK(I,2)
         IM = I
         IP1 = I+1
         DO 480 K = IP1,IGH
            IF (RMAX .GE. WK(K,2)) GO TO 480
            IM = K
            RMAX = WK(K,2)
  480    CONTINUE
         WK(I,1) = IM
         IF (IM .EQ. I) GO TO 500
         TEMP = WK(I,2)
         WK(I,2) = WK(IM,2)
         WK(IM,2) = TEMP
         DO 490 J = LOW,N
            TEMP = B(I,J)
            B(I,J) = B(IM,J)
            B(IM,J) = TEMP
            TEMP = A(I,J)
            A(I,J) = A(IM,J)
            A(IM,J) = TEMP
  490    CONTINUE
  500 CONTINUE
      WK(IGH,1) = IGH
  510 CONTINUE
      RETURN
C
C     LAST LINE OF GRADEQ
C
      END
