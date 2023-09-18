      SUBROUTINE TRNATA (NA,N,A)
C
C     *****PARAMETERS:
      INTEGER NA,N
      DOUBLE PRECISION A(NA,N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,NM1,JP1
      DOUBLE PRECISION TEMP
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE REPLACES THE N X N ARRAY A WITH THE TRANSPOSE
C     OF A.
C
C     *****PARAMETER DESCRIPTION:
C     ON INPUT:
C        NA               ROW DIMENSION OF THE ARRAY CONTAINING A AS
C                         DECLARED IN THE CALLING PROGRAM DIMENSION
C                         STATEMENT;
C
C        N                ORDER OF THE MATRIX A;
C
C        A                AN N X N MATRIX.
C
C     ON OUTPUT:
C
C        A                AN N X N ARRAY CONTAINING THE TRANSPOSE OF THE
C                         INPUT MATRIX A.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (ELEC. SYS. LAB., M.I.T., RM. 35-331,
C     CAMBRIDGE, MA 02139,  PH.: (617)-253-2125), SEPTEMBER 1977.
C     MOST RECENT VERSION: SEP. 21, 1977.
C
C     ------------------------------------------------------------------
C
      IF (N.EQ.1) RETURN
      NM1 = N- 1
      DO 20 J=1,NM1
         JP1=J+1
         DO 10 I=JP1,N
            TEMP=A(I,J)
            A(I,J)=A(J,I)
            A(J,I)=TEMP
10       CONTINUE
20    CONTINUE
      RETURN
C
C     LAST LINE OF TRNATA
C
      END
