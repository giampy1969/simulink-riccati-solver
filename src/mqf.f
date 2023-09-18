      SUBROUTINE MQF (NS,NX,NY,M,N,S,X,Y,WORK)
C
C     *****PARAMETERS:
      INTEGER NS,NX,NY,M,N
      DOUBLE PRECISION S(NS,M),X(NX,N),Y(NY,N),WORK(M)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,K,JM1
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COMPUTES THE SYMMETRIC MATRIX PRODUCT
C         T
C        X *S*X  WHERE S IS SYMMETRIC.  THE RESULT
C     IS STORED IN THE ARRAY Y.  X IS M X N, S IS M X M,
C     AND Y IS N X N.  THE ARRAY Y MUST BE DISTINCT FROM
C     BOTH X AND S.
C
C     *****PARAMETER DESCRIPTION:
C     ON INPUT:
C
C        NS,NX,NY         ROW DIMENSIONS OF THE ARRAYS CONTAINING S,X,
C                         AND Y, RESPECTIVELY, AS DECLARED IN THE
C                         PROGRAM DIMENSION STATEMENT;
C
C        M                NUMBER OF ROWS OF X AND ORDER OF THE MATRIX S;
C
C        N                NUMBER OF COLUMNS OF X AND ORDER OF THE
C                         MATRIX Y;
C
C        S                AN M X M SYMMETRIC MATRIX;
C
C        X                AN M X N MATRIX.
C
C     ON OUTPUT:
C
C                                                             T
C        Y                A SYMMETRIC N X N ARRAY CONTAINING X *S*X;
C
C        WORK             A REAL SCRATCH VECTOR OF LENGTH M.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (DEP'T. OF ELEC. AND COMP. ENGRG.,
C     UNIV. OF CALIF., SANTA BARBARA, CA 91306, PH.: (805) 961-3616),
C     SEPTEMBER 1977.
C     MOST RECENT VERSION: JUL. 29, 1983.
C
C     ------------------------------------------------------------------
C
C     COMPUTE THE LOWER TRIANGLE OF Y
C
      DO 60 K=1,N
         DO 10 I=1,M
            WORK(I) = 0.0D0
   10    CONTINUE
         DO 30 J=1,M
            DO 20 I=1,M
               WORK(I) = WORK(I) + S(I,J)*X(J,K)
   20       CONTINUE
   30    CONTINUE
         DO 50 I=K,N
            Y(I,K) = 0.0D0
            DO 40 J=1,M
               Y(I,K) = Y(I,K) + WORK(J)*X(J,I)
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
      IF (N.EQ.1) RETURN
C
C     DETERMINE THE STRICT UPPER TRIANGLE OF Y BY SYMMETRY
C
      DO 80 J=2,N
         JM1=J-1
         DO 70 I=1,JM1
            Y(I,J)=Y(J,I)
   70    CONTINUE
   80 CONTINUE
      RETURN
C
C     LAST LINE OF MQF
C
      END
