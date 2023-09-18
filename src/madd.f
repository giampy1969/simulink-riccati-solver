      SUBROUTINE MADD (NA,NB,NC,M,N,A,B,C)
C
C     *****PARAMETERS:
      INTEGER NA,NB,NC,M,N
      DOUBLE PRECISION A(NA,N),B(NB,N),C(NC,N)
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
C     THIS SUBROUTINE COMPUTES THE MATRIX SUM A+B AND STORES THE
C     RESULT IN THE ARRAY C.  ALL MATRICES ARE M X N.  THE SUM
C     MAY BE OVERWRITTEN INTO A (B) BY DESIGNATING THE ARRAY C
C     TO BE A (B).
C
C     *****PARAMETER DESCRIPTION:
C     ON INPUT:
C        NA,NB,NC         ROW DIMENSIONS OF THE ARRAYS CONTAINING A,B,
C                         AND C,RESPECTIVELY, AS DECLARED IN THE CALLING
C                         PROGRAM DIMENSION STATEMENT;
C
C        M                NUMBER OF ROWS OF THE MATRICES A, B, AND C;
C
C        N                NUMBER OF COLUMNS OF THE MATRICES A, B, AND C;
C
C        A                AN M X N MATRIX;
C
C        B                AN M X N MATRIX.
C
C     ON OUTPUT:
C
C        C                AN M X N ARRAY CONTAINING A+B.
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
            C(I,J)=A(I,J)+B(I,J)
10       CONTINUE
20    CONTINUE
      RETURN
C
C     LAST LINE OF MADD
C
      END
