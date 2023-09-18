      SUBROUTINE MULA(NA,NB,N,M,L,A,B,WORK)
C
C     *****PARAMETERS:
      INTEGER NA,NB,N,M,L
      DOUBLE PRECISION A(NA,M),B(NB,L),WORK(L)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,K
C
C     *****FORTRAN FUNCTIONS:
C     NONE.
C
C     *****SUBROUTINES CALLED:
C     NONE.
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COMPUTES THE MATRIX PRODUCT A * B AND OVERWRITES
C     IT INTO THE ARRAY A.  WHERE A IS N BY M AND B IS M BY L AND M IS
C     GREATER THAN OR EQUAL TO L.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C       NA,NB   INTEGER
C               ROW DIMENSIONS OF THE ARRAYS CONTAINING A AND B,
C               RESPECTIVELY, AS DECLARED IN THE MAIN CALLING PROGRAM
C               DIMENSION STATEMENT;
C
C       N       INTEGER
C               ROW DIMENSION OF THE MATRIX A;
C
C       M       INTEGER
C               COLUMN DIMENSION OF THE MATRIX A AND ROW DIMENSION OF
C               THE MATRIX B;
C
C       L       INTEGER
C               COLUMN DIMENSION OF THE MATRIX B;
C
C       A       REAL(NA,M)
C               AN N BY M MATRIX;
C
C       B       REAL(NB,L)
C               AN M BY L MATRIX.
C
C     ON OUTPUT:
C
C       A       CONTAINS THE N BY L MATRIX PRODUCT A * B.
C
C     *****ALGORITHM NOTES:
C     NONE.
C
C     *****HISTORY:
C     THIS SUBROUTINE WAS WRITTEN BY W.F. ARNOLD, NAVAL WEAPONS CENTER,
C     CODE 35104, CHINA LAKE, CA  93555, AS PART OF THE SOFTWARE PACKAGE
C     RICPACK, SEPTEMBER 1983.
C
C     ------------------------------------------------------------------
C
      DO 40 I=1,N
          DO 20 J=1,L
              WORK(J) = 0.0D0
              DO 10 K=1,M
                  WORK(J) = WORK(J) + A(I,K)*B(K,J)
   10         CONTINUE
   20     CONTINUE
          DO 30 J=1,L
              A(I,J) = WORK(J)
   30     CONTINUE
   40 CONTINUE
      RETURN
C
C     LAST LINE OF MULA
C
      END
