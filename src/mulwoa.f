      SUBROUTINE MULWOA (NA,NB,N,A,B,WORK)
C
C     *****PARAMETERS:
      INTEGER NA,NB,N
      DOUBLE PRECISION A(NA,N),B(NB,N),WORK(N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,K
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE OVERWRITES THE ARRAY A WITH THE MATRIX PRODUCT
C        A*B.  BOTH A AND B ARE N X N ARRAYS AND MUST BE DISTINCT.
C
C     *****PARAMETER DESCRIPTION:
C     ON INPUT:
C        NA,NB            ROW DIMENSIONS OF THE ARRAYS CONTAINING A AND
C                         B, RESPECTIVELY, AS DECLARED IN THE CALLING
C                         PROGRAM DIMENSION STATEMENT;
C
C        N                ORDER OF THE MATRICES A AND B;
C
C        A                AN N X N MATRIX;
C
C        B                AN N X N MATRIX.
C
C     ON OUTPUT:
C
C        A                AN N X N ARRAY CONTAINING A*B;
C
C        WORK             A REAL SCRATCH VECTOR OF LENGTH N.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (ELEC. SYS. LAB., M.I.T., RM. 35-331,
C     CAMBRIDGE, MA 02139,  PH.: (617)-253-2125), SEPTEMBER 1977.
C     MOST RECENT VERSION: SEP. 21, 1977.
C
C     ------------------------------------------------------------------
C
      DO 40 I=1,N
         DO 20 J=1,N
            WORK(J)=0.0D0
            DO 10 K=1,N
               WORK(J)=WORK(J)+A(I,K)*B(K,J)
10          CONTINUE
20       CONTINUE
         DO 30 J=1,N
            A(I,J)=WORK(J)
30       CONTINUE
40    CONTINUE
      RETURN
C
C     LAST LINE OF MULWOA
C
      END
