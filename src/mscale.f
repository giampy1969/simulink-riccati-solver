      SUBROUTINE MSCALE (NA,M,N,ALPHA,A)
C
C     *****PARAMETERS:
      INTEGER NA,M,N
      DOUBLE PRECISION ALPHA,A(NA,N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE REPLACES THE M X N ARRAY A WITH (ALPHA*A)
C     WHERE ALPHA IS A (REAL) SCALAR.
C
C     *****PARAMETER DESCRIPTION:
C     ON INPUT:
C        NA               ROW DIMENSION OF THE ARRAY CONTAINING A AS
C                         DECLARED IN THE CALLING PROGRAM DIMENSION
C                         STATEMENT;
C
C        M                NUMBER OF ROWS OF THE MATRIX A;
C
C        N                NUMBER OF COLUMNS OF THE MATRIX A;
C
C        ALPHA            THE SCALAR MULTIPLIER;
C
C        A                AN M X N MATRIX.
C
C     ON OUTPUT:
C
C        A                THE M X N ARRAY CONTAINING ALPHA*A.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (ELEC. SYS. LAB., M.I.T., RM. 35-331,
C     CAMBRIDGE, MA 02139,  PH.: (617)-253-2125), SEPTEMBER 1977.
C     MOST RECENT VERSION: SEP. 21, 1977.
C
C     ------------------------------------------------------------------
C
      DO 20 J=1,N
         DO 10 I=1,M
            A(I,J)=ALPHA*A(I,J)
10       CONTINUE
20    CONTINUE
      RETURN
C
C     LAST LINE OF MSCALE
C
      END
