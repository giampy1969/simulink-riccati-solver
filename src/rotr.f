      SUBROUTINE ROTR (H,NMAX,N,L1,L2,M1,M2,S,C)
C
C     *****PARAMETERS:
      INTEGER L1,L2,M1,M2,N,NMAX
      DOUBLE PRECISION C,S,H(NMAX,N)
C
C     *****LOCAL VARIABLES:
      INTEGER J
      DOUBLE PRECISION T
C
C     *****FORTRAN FUNCTIONS:
C     NONE
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ________________________________________________________________
C
C     *****PURPOSE:
C     THIS ROUTINE PERFORMS THE GIVENS ROTATION:
C
C                       \ S  C \
C                       \-C  S \
C
C     ON ROWS L1 AND L2 OF H, THIS FROM COLUMNS M1 TO M2.
C
C     *****PARAMETER DESCRIPTION:
C     ON INPUT:
C       NMAX    INTEGER
C               ROW DIMENSION OF THE ARRAYS CONTAINING A,B,Z AS
C               DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT;
C
C       N       INTEGER
C               ORDER OF THE MATRICES A,B,Z;
C
C       H       REAL(NMAX,N)
C               THE ARRAY UPON WHICH THE GIVENS ROTATION IS TO BE
C               PERFORMED;
C
C       L1,L2   INTEGER
C               ROWS INVOLVED IN THE ROTATION;
C
C       M1,M2   INTEGER
C               COLUMNS INVOLVED IN THE ROTATION;
C
C       C,S     REAL
C               ROTATION PARAMETERS.
C
C     ON OUTPUT:
C
C       H       THE ARRAY WHOSE ROWS HAVE BEEN MODIFIED AS SPECIFIED
C               ABOVE.
C
C     *****HISTORY:
C     WRITTEN BY P. VAN DOOREN("A GENERALIZED EIGENVALUE APPROACH
C     FOR SOLVING RICCATI EQUATIONS", INTERNAL REPORT NA-80-02,
C     DEPT. OF COMPUTER SCIENCE, STANFORD UNIVERSITY, 1980).
C
C     ----------------------------------------------------------------
C
      DO 10 J = M1,M2
         T = S*H(L1,J) + C*H(L2,J)
         H(L2,J) = -C*H(L1,J) + S*H(L2,J)
         H(L1,J) = T
   10 CONTINUE
      RETURN
C
C     LAST LINE OF ROTR
C
      END
