      SUBROUTINE GIV (A,B,D,E)
C
C     *****PARAMETERS:
      DOUBLE PRECISION A,B,D,E
C
C     *****LOCAL VARIABLES:
      DOUBLE PRECISION C
C
C     *****FORTRAN FUNCTIONS:
      DOUBLE PRECISION DABS, DMAX1
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ----------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COMPUTES:
C
C               D = A / SQRT(A*A + B*B)
C
C               E = B / SQRT(A*A + B*B)
C
C     BY FIRST NORMALIZING BY THE LARGEST MAGNITUDE OF A OR B.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C       A       REAL;
C
C       B       REAL.
C
C     ON OUTPUT:
C
C       D       REAL
C               = A / SQRT(A*A + B*B) IF A AND B DON'T EQUAL ZERO;
C               = 1 IF A AND B EQUAL ZERO.
C
C       E       REAL
C               = B / SQRT(A*A + B*B) IF A AND B DON'T EQUAL ZERO;
C               = 0 IF A AND B EQUAL ZERO.
C
C     *****HISTORY:
C     WRITTEN BY P. VAN DOOREN("A GENERALIZED EIGENVALUE APPROACH
C     FOR SOLVING RICCATI EQUATIONS", INTERNAL REPORT NA-80-02,
C     DEPT. OF COMPUTER SCIENCE, STANFORD UNIVERSITY, 1980).
C
C     * MODIFIED 1/7/86 BY J. D. BIRDWELL, UNIVERSITY OF TENNESSEE,
C       TO DEFINE A ROTATION OF ZERO ANGLE IN CASE BOTH A AND B ARE ZERO.
C       THIS "CORRECTS" A PATHOLOGICAL CASE IN RICSOL.  NORMAL OPERATION
C       IS NOT MODIFIED.  THE MODIFICATION USES FORTRAN-77.
C
C     ----------------------------------------------------------------
C
      C = DMAX1(DABS(A),DABS(B))
      IF (C .NE. 0.D0) THEN
        D = A/C
        E = B/C
        C = DSQRT(D*D + E*E)
        D = D/C
        E = E/C
      ELSE
        D = 1.D0
        E = 0.D0
      END IF
      RETURN
C
C     LAST LINE OF GIV
C
      END
