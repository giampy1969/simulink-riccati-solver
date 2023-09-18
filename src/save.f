      SUBROUTINE SAVE (NA,NAS,M,N,A,ASAVE)
C
C     *****PARAMETERS:
      INTEGER NA,NAS,M,N
      DOUBLE PRECISION A(NA,N),ASAVE(NAS,N)
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
C     THIS SUBROUTINE COPIES THE CONTENTS OF THE M X N ARRAY A INTO
C     THE M X N ARRAY ASAVE.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C        NA,NAS           ROW DIMENSIONS OF THE ARRAYS CONTAINING A AND
C                         AS, RESPECTIVELY, AS DECLARED IN THE CALLING
C                         PROGRAM DIMENSION STATEMENT;
C
C        M                NUMBER OF ROWS OF THE MATRICES A AND ASAVE;
C
C        N                NUMBER OF COLUMNS OF THE MATRICES A AND ASAVE;
C
C        A                AN M X N MATRIX.
C
C     ON OUTPUT:
C
C        ASAVE            AN M X N ARRAY CONTAINING THE ARRAY A.
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
            ASAVE(I,J)=A(I,J)
10       CONTINUE
20    CONTINUE
      RETURN
C
C     LAST LINE OF SAVE
C
      END
