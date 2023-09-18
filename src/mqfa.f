      SUBROUTINE MQFA (NX,NY,M,N,X,Y,WORK)
C
C     *****PARAMETERS:
      INTEGER NX,NY,M,N
      DOUBLE PRECISION X(NX,N),Y(NY,N),WORK(M)
C
C     LOCAL VARIABLES:
      INTEGER I,J,JM1,K
C
C     *****FORTRAN FUNCTIONS:
C     NONE
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     -------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COMPUTES THE SYMMETRIC MATRIX PRODUCT
C         T
C        X * X = Y.  X IS M X N AND Y IS N X N.  THE ARRAY Y MUST
C     BE DISTINCT FROM X.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C        NX,NY            ROW DIMENSIONS OF THE ARRAYS CONTAINING X
C                         AND Y, RESPECTIVELY, AS DECLARED IN THE MAIN
C                         PROGRAM DIMENSION STATEMENT;
C
C        M                THE NUMBER OF ROWS OF X;
C
C        N                NUMBER OF COLUMNS OF X AND THE ORDER OF THE
C                         MATRIX Y;
C
C        X                AN M X N MATRIX.
C
C     ON OUTPUT:
C
C                                                              T
C        Y                A SYMMETRIC N X N MATRIX CONTAINING X *X;
C
C        WORK             A REAL SCRATCH VECTOR OF LENGTH M.
C
C     *****HISTORY:
C     WRITTEN BY W.F. ARNOLD (NAVAL WEAPONS CENTER, CODE 35104, CHINA LAKE,
C     CA.  93555, PH.:  (619)-939-3900), FEBRUARY 1984.
C     MOST RECENT VERSION:  FEBRUARY 1984.
C
C     -----------------------------------------------------------------------
C
C     COMPUTE THE LOWER TRIANGLE OF Y
C
      DO 40 K=1,N
         DO 10 I=1,M
            WORK(I) = X(I,K)
   10    CONTINUE
         DO 30 J=K,N
            Y(J,K) = 0.0D0
            DO 20 I=1,M
               Y(J,K) = Y(J,K) + WORK(I)*X(I,J)
   20       CONTINUE
   30    CONTINUE
   40 CONTINUE
C
C     DETERMINE THE STRICT UPPER TRIANGLE OF Y BY SYMMETRY
C
      DO 60 J=2,N
         JM1=J-1
         DO 50 I=1,JM1
            Y(I,J) = Y(J,I)
   50    CONTINUE
   60 CONTINUE
      RETURN
C
C     LAST LINE OF MQFA
C
      END
