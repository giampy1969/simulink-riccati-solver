      SUBROUTINE REDUCE (N,MA,A,MB,B,LOW,IGH,CSCALE,WK)
C
C     *****PARAMETERS:
      INTEGER IGH,LOW,MA,MB,N
      DOUBLE PRECISION A(MA,N),B(MB,N),CSCALE(N),WK(N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,IFLOW,II,IP1,IS,J,JP1,K,L,LM1,M
      DOUBLE PRECISION F
C
C     *****FUNCTIONS:
C     NONE
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ---------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE REDUCES, IF POSSIBLE, THE ORDER OF THE
C     GENERALIZED EIGENVALUE PROBLEM A*X = (LAMBDA)*B*X BY PERMUTING
C     THE ROWS AND COLUMNS OF A AND B SO THAT THEY EACH HAVE THE
C     FORM
C                       U  X  Y
C                       0  C  Z
C                       0  0  R
C
C     WHERE U AND R ARE UPPER TRIANGULAR AND C, X, Y, AND Z ARE
C     ARBITRARY.  THUS, THE ISOLATED EIGENVALUES CORRESPONDING TO
C     THE TRIANGULAR MATRICES ARE OBTAINED BY A DIVISION, LEAVING
C     ONLY EIGENVALUES CORRESPONDING TO THE CENTER MATRICES TO BE
C     COMPUTED.
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
C               DEFINED ABOVE.
C
C     ON OUTPUT:
C
C       A,B     CONTAIN THE PERMUTED A AND B MATRICES;
C
C       LOW     INTEGER
C               BEGINNING INDEX OF THE SUBMATRICES OF A AND B
C               CONTAINING THE NON-ISOLATED EIGENVALUES;
C
C       IGH     INTEGER
C               ENDING INDEX OF THE SUBMATRICES OF A AND B
C               CONTAINING THE NON-ISOLATED EIGENVALUES.  IF
C               IGH = 1 (LOW = 1 ALSO), THE PERMUTED A AND B
C               MATRICES ARE UPPER TRIANGULAR;
C
C       CSCALE  REAL(N)
C               CONTAINS THE REQUIRED COLUMN PERMUTATIONS IN ITS
C               FIRST LOW-1 AND ITS IGH+1 THROUGH N LOCATIONS;
C
C       WK      REAL(N)
C               CONTAINS THE REQUIRED ROW PERMUTATIONS IN ITS FIRST
C               LOW-1 AND ITS IGH+1 THROUGH N LOCATIONS.
C
C     *****ALGORITHM NOTES:
C     NONE
C
C     *****HISTORY:
C     WRITTEN BY R. C. WARD.......
C
C     ---------------------------------------------------------------
C
      K = 1
      L = N
      GO TO 20
C
C     FIND ROW WITH ONE NONZERO IN COLUMNS 1 THROUGH L
C
   10 CONTINUE
      L = LM1
      IF (L .NE. 1) GO TO 20
      WK(1) = 1
      CSCALE(1) = 1
      GO TO 200
   20 CONTINUE
      LM1 = L-1
      DO 70 II = 1,L
         I = L+1-II
         DO 30 J = 1,LM1
            JP1 = J+1
            IF (A(I,J) .NE. 0.0D0 .OR. B(I,J) .NE. 0.0D0) GO TO 40
   30    CONTINUE
         J = L
         GO TO 60
   40    CONTINUE
         DO 50 J = JP1,L
            IF (A(I,J) .NE. 0.0D0 .OR. B(I,J) .NE. 0.0D0) GO TO 70
   50    CONTINUE
         J = JP1-1
   60    CONTINUE
         M = L
         IFLOW = 1
         GO TO 150
   70 CONTINUE
      GO TO 90
C
C     FIND COLUMN WITH ONE NONZERO IN ROWS K THROUGH N
C
   80 CONTINUE
      K = K+1
   90 CONTINUE
      DO 140 J = K,L
         DO 100 I = K,LM1
            IP1 = I+1
            IF (A(I,J) .NE. 0.0D0 .OR. B(I,J) .NE. 0.0D0) GO TO 110
  100    CONTINUE
         I = L
         GO TO 130
  110    CONTINUE
         DO 120 I = IP1,L
            IF (A(I,J) .NE. 0.0D0 .OR. B(I,J) .NE. 0.0D0) GO TO 140
  120    CONTINUE
         I = IP1-1
  130    CONTINUE
         M = K
         IFLOW = 2
         GO TO 150
  140 CONTINUE
      GO TO 200
C
C     PERMUTE ROWS M AND I
C
  150 CONTINUE
      WK(M) = I
      IF (I .EQ. M) GO TO 170
      DO 160 IS = K,N
         F = A(I,IS)
         A(I,IS) = A(M,IS)
         A(M,IS) = F
         F = B(I,IS)
         B(I,IS) = B(M,IS)
         B(M,IS) = F
  160 CONTINUE
C
C     PERMUTE COLUMNS M AND J
C
  170 CONTINUE
      CSCALE(M) = J
      IF (J .EQ. M) GO TO 190
      DO 180 IS = 1,L
         F = A(IS,J)
         A(IS,J) = A(IS,M)
         A(IS,M) = F
         F = B(IS,J)
         B(IS,J) = B(IS,M)
         B(IS,M) = F
  180 CONTINUE
  190 CONTINUE
      GO TO (10,80), IFLOW
  200 CONTINUE
      LOW = K
      IGH = L
      RETURN
C
C     LAST LINE OF REDUCE
C
      END
