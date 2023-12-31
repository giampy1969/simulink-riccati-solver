      SUBROUTINE SCALBK (N,MZ,Z,M,LOW,IGH,CSCALE)
C
C     *****PARAMETERS:
      INTEGER IGH,LOW,M,MZ,N
      DOUBLE PRECISION CSCALE(N),Z(MZ,N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,IGHP1,II,IR,J,K,LOWM1
      DOUBLE PRECISION FI,TEMP
C
C     *****FORTRAN FUNCTIONS:
C     NONE
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ---------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE BACK TRANSFORMS THE EIGENVECTORS OF A
C     GENERALIZED EIGENVALUE PROBLEM A*X = (LAMBDA)*B*X BALANCED BY
C     SUBROUTINES REDUCE AND/OR SCALEG TO THOSE OF THE ORIGINAL
C     PROBLEM.
C     REF.:  WARD, R. C., BALANCING THE GENERALIZED EIGENVALUE
C     PROBLEM, SIAM J. SCI. STAT. COMPUT., VOL. 2, NO. 2, JUNE 1981,
C     141-152.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C       MZ      INTEGER
C               ROW DIMENSION OF THE ARRAY Z AS SPECIFIED IN THE MAIN
C               CALLING PROGRAM DIMENSION STATEMENT;
C
C       N       INTEGER
C               ORDER OF THE MATRICES A AND B IN THE EIGENPROBLEM;
C
C       M       INTEGER
C               SPECIFIES THE NUMBER OF EIGENVECTORS TO BE TRANS-
C               FORMED;
C
C       Z       REAL(MZ,N)
C               CONTAINS THE EIGENVECTORS TO BE TRANSFORMED;
C
C       LOW     INTEGER
C               SPECIFIES THE BEGINNING INDEX OF THE SUBMATRICES OF
C               A AND B WHICH WERE SCALED;
C
C       IGH     INTEGER
C               SPECIFIES THE ENDING INDEX OF THE SUBMATRICES OF
C               A AND B WHICH WERE SCALED;
C
C       CSCALE  REAL(N)
C               CONTAINS THE REDUCING COLUMN PERMUTATIONS AND SCALING
C               INFORMATION AS RETURNED FROM REDUCE AND SCALEG.  IF
C               REDUCE WAS NOT CALLED, SET LOW TO 1 AND IGH TO N.
C               IF SCALEG WAS NOT CALLED, SET THE SCALE FACTORS IN
C               CSCALE LOCATIONS LOW THROUGH IGH TO ZERO.
C
C     ON OUTPUT:
C
C       Z       CONTAINS THE TRANSFORMED EIGENVECTORS.
C
C     *****ALGORITHM NOTES:
C     NONE.
C
C     *****HISTORY:
C     WRITTEN BY R. C. WARD.......
C
C     ---------------------------------------------------------------
C
      IF (LOW .EQ. IGH) GO TO 545
C
C     APPLY SCALING TRANSFORMATION
C
      DO 540 I = LOW,IGH
         IR = CSCALE(I)
         FI = 2**IR
         DO 530 J = 1,M
            Z(I,J) = Z(I,J)*FI
  530    CONTINUE
  540 CONTINUE
  545 CONTINUE
C
C     APPLY REDUCING COLUMN PERMUTATIONS
C
      IF (LOW .EQ. 1) GO TO 570
      LOWM1 = LOW-1
      DO 560 II = 1,LOWM1
         I = LOW-II
         K = CSCALE(I)
         IF (K .EQ. I) GO TO 560
         DO 550 J = 1,M
            TEMP = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = TEMP
  550    CONTINUE
  560 CONTINUE
  570 CONTINUE
      IF (IGH .EQ. N) GO TO 600
      IGHP1 = IGH+1
      DO 590 I = IGHP1,N
         K = CSCALE(I)
         IF (K .EQ. I) GO TO 590
         DO 580 J = 1,M
            TEMP = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = TEMP
  580    CONTINUE
  590 CONTINUE
  600 CONTINUE
      RETURN
C
C     LAST LINE OF SCALBK
C
      END
