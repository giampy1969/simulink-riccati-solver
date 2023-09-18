      SUBROUTINE EXCHQZ (A,B,Z,NMAX,N,L,LS1,LS2,EPS,IFAIL)
C
C     *****PARAMETERS:
      INTEGER IFAIL,L,LS1,LS2,N,NMAX
      DOUBLE PRECISION EPS,A(NMAX,N),B(NMAX,N),Z(NMAX,N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,IT1,IT2,J,LI,LJ,LL,L1,L2,L3
      DOUBLE PRECISION AMMBMM,AMNBNN,ANMBMM,ANNBNN,A11B11,A12B22,
     X                 A21B11,A22B22,BMNBNN,B12B22,D,E,F,G,SA,SB,
     X                 U(3,3)
      LOGICAL ALTB
C
C     *****FORTRAN FUNCTIONS:
      DOUBLE PRECISION DABS,DMAX1
C
C     *****SUBROUTINES CALLED:
C     GIV, ROTC, ROTR
C
C     ----------------------------------------------------------------
C
C     *****PURPOSE:
C     GIVEN THE UPPER TRIANGULAR MATRIX B AND THE UPPER HESSENBERG
C     MATRIX A WITH CONSECUTIVE LS1XLS1 AND LS2XLS2 DIAGONAL BLOCKS
C     (LS1,LS2 .LE. 2) STARTING AT ROW/COLUMN L, THIS SUBROUTINE
C     PRODUCES EQUIVALENCE TRANSFORMATIONS QT AND ZT THAT EXCHANGE
C     THE BLOCKS ALONG WITH THE GENERALIZED EIGENVALUES CORRESPONDING
C     TO THE REGULAR MATRIX PENCIL A - LAMBDA*B.
C
C     REF.:  VAN DOOREN, P., A GENERALIZED EIGENVALUE APPROACH FOR
C     SOLVING RICCATI EQUATIONS, SIAM J. SCI. STAT. COMPUT.,
C     VOL. 2, NO. 2, JUNE 1981, 121-135.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C       NMAX    INTEGER
C               ROW DIMENSION OF THE ARRAYS CONTAINING A,B,Z AS
C               DECLARED IN THE MAIN CALLING PROGRAM DIMENSION STATEMENT;
C
C       N       INTEGER
C               ORDER OF THE MATRICES A,B,Z;
C
C       A       REAL(NMAX,N)
C               UPPER HESSENBERG MATRIX WITH 1X1 OR 2X2 DIAGONAL
C               BLOCKS.  ELEMENTS OUTSIDE THE UPPER HESSENBERG
C               STRUCTURE ARE ARBITRARY;
C
C       B       REAL(NMAX,N)
C               UPPER TRIANGULAR MATRIX.  ELEMENTS OUTSIDE THE
C               UPPER TRIANGULAR STRUCTURE ARE ARBITRARY;
C
C       L       INTEGER
C               THE STARTING POSITION OF THE BLOCKS BEING INTERCHANGED;
C
C       LS1     INTEGER
C               THE SIZE OF THE FIRST BLOCK (.LE. 2);
C
C       LS2     INTEGER
C               THE SIZE OF THE SECOND BLOCK (.LE. 2);
C
C       EPS     REAL
C               REQUIRED ABSOLUTE ACCURACY OF THE RESULTS.  NORMALLY
C               EQUAL TO THE MACHINE PRECISION;
C
C     ON OUTPUT:
C
C       A,B     UPPER HESSENBERG MATRIX, UPPER TRIANGULAR MATRIX WITH
C               THE LS1 X LS1 AND LS2 X LS2 DIAGONAL BLOCKS STARTING
C               IN POSITION L INTERCHANGED;
C
C       Z       REAL(NMAX,N)
C               THIS ARRAY IS OVERWRITTEN BY THE PRODUCT OF THE
C               CONTENTS OF THE ARRAY Z(UPON ENTRY INTO THIS
C               SUBROUTINE), AND THE COLUMN TRANSFORMATIONS ZT
C               (CALCULATED BY THIS SUBROUTINE);
C
C       IFAIL   INTEGER
C               ERROR FLAG
C               = 1  INDICATES ATTEMPTED REORDERING FAILED
C               = 0  NORMAL RETURN.
C
C     *****ALGORITHM NOTES:
C     NONE
C
C     *****HISTORY:
C     WRITTEN BY P. VAN DOOREN("A GENERALIZED EIGENVALUE APPROACH
C     FOR SOLVING RICCATI EQUATIONS", INTERNAL REPORT NA-80-02,
C     DEPT. OF COMPUTER SCIENCE, STANFORD UNIVERSITY, 1980).
C
C     ----------------------------------------------------------------
C
      IFAIL = 0
      L1 = L+1
      LL = LS1+LS2
      IF (LL .GT. 2) GO TO 50
C
C     INTERCHANGE 1X1 AND 1X1 BLOCKS
C
      F = DMAX1(DABS(A(L1,L1)),DABS(B(L1,L1)))
      ALTB = .TRUE.
      IF (DABS(A(L1,L1)) .GE. F) ALTB = .FALSE.
      SA = A(L1,L1)/F
      SB = B(L1,L1)/F
      F = SA*B(L,L) - SB*A(L,L)
C
C     COMPUTE Z
C
      G = SA*B(L,L1) - SB*A(L,L1)
      CALL GIV (F,G,D,E)
      CALL ROTC (A,NMAX,N,L,L1,1,L1,D,E)
      CALL ROTC (B,NMAX,N,L,L1,1,L1,D,E)
      CALL ROTC (Z,NMAX,N,L,L1,1,N,D,E)
C
C     COMPUTE Q
C
      IF (ALTB) CALL GIV (B(L,L),B(L1,L),D,E)
      IF (.NOT. ALTB) CALL GIV (A(L,L),A(L1,L),D,E)
      CALL ROTR (A,NMAX,N,L,L1,L,N,D,E)
      CALL ROTR (B,NMAX,N,L,L1,L,N,D,E)
      A(L1,L) = 0.0D0
      B(L1,L) = 0.0D0
      RETURN
C
C     INTERCHANGE 1X1 AND 2X2 BLOCKS
C
   50 CONTINUE
      L2 = L+2
      IF (LS1 .EQ. 2) GO TO 100
      G = DMAX1(DABS(A(L,L)),DABS(B(L,L)))
      ALTB = .TRUE.
      IF (DABS(A(L,L)) .LT. G) GO TO 60
      ALTB = .FALSE.
      CALL GIV (A(L1,L1),A(L2,L1),D,E)
      CALL ROTR (A,NMAX,N,L1,L2,L1,N,D,E)
      CALL ROTR (B,NMAX,N,L1,L2,L1,N,D,E)
C
C     COMPUTE Q AND Z
C
   60 CONTINUE
      SA = A(L,L)/G
      SB = B(L,L)/G
      DO 80 J =1,2
         LJ = L+J
         DO 70 I = 1,3
            LI = L+I-1
            U(I,J) = SA*B(LI,LJ) - SB*A(LI,LJ)
   70    CONTINUE
   80 CONTINUE
      CALL GIV (U(3,1),U(3,2),D,E)
      CALL ROTC (U,3,3,1,2,1,3,D,E)
C
C     Q1
C
      CALL GIV (U(1,1),U(2,1),D,E)
      U(2,2) = -U(1,2)*E + U(2,2)*D
      CALL ROTR (A,NMAX,N,L,L1,L,N,D,E)
      CALL ROTR (B,NMAX,N,L,L1,L,N,D,E)
C
C     Z1
C
      IF (ALTB) CALL GIV (B(L1,L),B(L1,L1),D,E)
      IF (.NOT. ALTB) CALL GIV (A(L1,L),A(L1,L1),D,E)
      CALL ROTC (A,NMAX,N,L,L1,1,L2,D,E)
      CALL ROTC (B,NMAX,N,L,L1,1,L2,D,E)
      CALL ROTC (Z,NMAX,N,L,L1,1,N,D,E)
C
C     Q2
C
      CALL GIV (U(2,2),U(3,2),D,E)
      CALL ROTR (A,NMAX,N,L1,L2,L,N,D,E)
      CALL ROTR (B,NMAX,N,L1,L2,L,N,D,E)
C
C     Z2
C
      IF (ALTB) CALL GIV (B(L2,L1),B(L2,L2),D,E)
      IF (.NOT. ALTB) CALL GIV (A(L2,L1),A(L2,L2),D,E)
      CALL ROTC (A,NMAX,N,L1,L2,1,L2,D,E)
      CALL ROTC (B,NMAX,N,L1,L2,1,L2,D,E)
      CALL ROTC (Z,NMAX,N,L1,L2,1,N,D,E)
      IF (ALTB) GO TO 90
      CALL GIV (B(L,L),B(L1,L),D,E)
      CALL ROTR (A,NMAX,N,L,L1,L,N,D,E)
      CALL ROTR (B,NMAX,N,L,L1,L,N,D,E)
   90 CONTINUE
      A(L2,L) = 0.0D0
      A(L2,L1) = 0.0D0
      B(L1,L) = 0.0D0
      B(L2,L) = 0.0D0
      B(L2,L1) = 0.0D0
      RETURN
C
C     INTERCHANGE 2X2 AND 1X1 BLOCKS
C
  100 CONTINUE
      IF (LS2 .EQ. 2) GO TO 150
      G = DMAX1(DABS(A(L2,L2)),DABS(B(L2,L2)))
      ALTB = .TRUE.
      IF (DABS(A(L2,L2)) .LT. G) GO TO 120
      ALTB = .FALSE.
      CALL GIV (A(L,L),A(L1,L),D,E)
      CALL ROTR (A,NMAX,N,L,L1,L,N,D,E)
      CALL ROTR (B,NMAX,N,L,L1,L,N,D,E)
C
C     COMPUTE Q AND Z
C
  120 CONTINUE
      SA = A(L2,L2)/G
      SB = B(L2,L2)/G
      DO 130 I = 1,2
         LI = L+I-1
         DO 125 J = 1,3
            LJ = L+J-1
            U(I,J) = SA*B(LI,LJ) - SB*A(LI,LJ)
  125    CONTINUE
  130 CONTINUE
      CALL GIV (U(1,1),U(2,1),D,E)
      CALL ROTR (U,3,3,1,2,1,3,D,E)
C
C     Z1
C
      CALL GIV (U(2,2),U(2,3),D,E)
      U(1,2) = U(1,2)*E - U(1,3)*D
      CALL ROTC (A,NMAX,N,L1,L2,1,L2,D,E)
      CALL ROTC (B,NMAX,N,L1,L2,1,L2,D,E)
      CALL ROTC (Z,NMAX,N,L1,L2,1,N,D,E)
C
C     Q1
C
      IF (ALTB) CALL GIV (B(L1,L1),B(L2,L1),D,E)
      IF (.NOT. ALTB) CALL GIV (A(L1,L1),A(L2,L1),D,E)
      CALL ROTR (A,NMAX,N,L1,L2,L,N,D,E)
      CALL ROTR (B,NMAX,N,L1,L2,L,N,D,E)
C
C     Z2
C
      CALL GIV (U(1,1),U(1,2),D,E)
      CALL ROTC (A,NMAX,N,L,L1,1,L2,D,E)
      CALL ROTC (B,NMAX,N,L,L1,1,L2,D,E)
      CALL ROTC (Z,NMAX,N,L,L1,1,N,D,E)
C
C     Q2
C
      IF (ALTB) CALL GIV (B(L,L),B(L1,L),D,E)
      IF (.NOT. ALTB) CALL GIV (A(L,L),A(L1,L),D,E)
      CALL ROTR (A,NMAX,N,L,L1,L,N,D,E)
      CALL ROTR (B,NMAX,N,L,L1,L,N,D,E)
      IF (ALTB) GO TO 140
      CALL GIV (B(L1,L1),B(L2,L1),D,E)
      CALL ROTR (A,NMAX,N,L1,L2,L1,N,D,E)
      CALL ROTR (B,NMAX,N,L1,L2,L1,N,D,E)
  140 CONTINUE
      A(L1,L) = 0.0D0
      A(L2,L) = 0.0D0
      B(L1,L) = 0.0D0
      B(L2,L) = 0.0D0
      B(L2,L1) = 0.0D0
      RETURN
C
C     INTERCHANGE 2X2 AND 2X2 BLOCKS
C
  150 CONTINUE
      L3 = L+3
      AMMBMM = A(L,L)/B(L,L)
      ANMBMM = A(L1,L)/B(L,L)
      AMNBNN = A(L,L1)/B(L1,L1)
      ANNBNN = A(L1,L1)/B(L1,L1)
      BMNBNN = B(L,L1)/B(L1,L1)
      DO 180 IT1 = 1,3
         U(1,1) = 1.0D0
         U(2,1) = 1.0D0
         U(3,1) = 1.0D0
         DO 170 IT2 = 1,10
C
C     Q1, Q2
C
            CALL GIV (U(2,1),U(3,1),D,E)
            CALL ROTR (A,NMAX,N,L1,L2,L,N,D,E)
            CALL ROTR (B,NMAX,N,L1,L2,L1,N,D,E)
            U(2,1) = D*U(2,1) + E*U(3,1)
            CALL GIV (U(1,1),U(2,1),D,E)
            CALL ROTR (A,NMAX,N,L,L1,L,N,D,E)
            CALL ROTR (B,NMAX,N,L,L1,L,N,D,E)
C
C     Z1, Z2
C
            CALL GIV (B(L2,L1),B(L2,L2),D,E)
            CALL ROTC (A,NMAX,N,L1,L2,1,L3,D,E)
            CALL ROTC (B,NMAX,N,L1,L2,1,L2,D,E)
            CALL ROTC (Z,NMAX,N,L1,L2,1,N,D,E)
            CALL GIV (B(L1,L),B(L1,L1),D,E)
            CALL ROTC (A,NMAX,N,L,L1,1,L3,D,E)
            CALL ROTC (B,NMAX,N,L,L1,1,L1,D,E)
            CALL ROTC (Z,NMAX,N,L,L1,1,N,D,E)
C
C     Q3, Z3, Q4, Z4, Q5, Z5
C
            CALL GIV (A(L2,L),A(L3,L),D,E)
            CALL ROTR (A,NMAX,N,L2,L3,L,N,D,E)
            CALL ROTR (B,NMAX,N,L2,L3,L2,N,D,E)
            CALL GIV (B(L3,L2),B(L3,L3),D,E)
            CALL ROTC (A,NMAX,N,L2,L3,1,L3,D,E)
            CALL ROTC (B,NMAX,N,L2,L3,1,L3,D,E)
            CALL ROTC (Z,NMAX,N,L2,L3,1,N,D,E)
            CALL GIV (A(L1,L),A(L2,L),D,E)
            CALL ROTR (A,NMAX,N,L1,L2,L,N,D,E)
            CALL ROTR (B,NMAX,N,L1,L2,L1,N,D,E)
            CALL GIV (B(L2,L1),B(L2,L2),D,E)
            CALL ROTC (A,NMAX,N,L1,L2,1,L3,D,E)
            CALL ROTC (B,NMAX,N,L1,L2,1,L2,D,E)
            CALL ROTC (Z,NMAX,N,L1,L2,1,N,D,E)
            CALL GIV (A(L2,L1),A(L3,L1),D,E)
            CALL ROTR (A,NMAX,N,L2,L3,L1,N,D,E)
            CALL ROTR (B,NMAX,N,L2,L3,L2,N,D,E)
            CALL GIV (B(L3,L2),B(L3,L3),D,E)
            CALL ROTC (A,NMAX,N,L2,L3,1,L3,D,E)
            CALL ROTC (B,NMAX,N,L2,L3,1,L3,D,E)
            CALL ROTC (Z,NMAX,N,L2,L3,1,N,D,E)
C
C     TEST FOR CONVERGENCE
C
            IF (DABS(A(L2,L1)) .LE. EPS) GO TO 190
            A11B11 = A(L,L)/B(L,L)
            A12B22 = A(L,L1)/B(L1,L1)
            A21B11 = A(L1,L)/B(L,L)
            A22B22 = A(L1,L1)/B(L1,L1)
            B12B22 = B(L,L1)/B(L1,L1)
            U(1,1) = ((AMMBMM-A11B11)*(ANNBNN-A11B11) - AMNBNN*ANMBMM
     +               + ANMBMM*BMNBNN*A11B11)/A21B11+A12B22-
     +               A11B11*B12B22
            U(2,1) = (A22B22-A11B11)-A21B11*B12B22-(AMMBMM-A11B11)
     +               -(ANNBNN-A11B11) + ANMBMM*BMNBNN
            U(3,1) = A(L2,L1)/B(L1,L1)
  170    CONTINUE
  180 CONTINUE
      IFAIL = 1
      RETURN
  190 CONTINUE
      A(L2,L) = 0.0D0
      A(L2,L1) = 0.0D0
      A(L3,L) = 0.0D0
      A(L3,L1) = 0.0D0
      B(L1,L) = 0.0D0
      B(L2,L) = 0.0D0
      B(L2,L1) = 0.0D0
      B(L3,L) = 0.0D0
      B(L3,L1) = 0.0D0
      B(L3,L2) = 0.0D0
      RETURN
C
C     LAST LINE OF EXCHQZ
C
      END
