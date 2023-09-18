      SUBROUTINE MQFWO (NS,NX,N,S,X,WORK)
C
C     *****PARAMETERS:
      INTEGER NS,NX,N
      DOUBLE PRECISION S(NS,N),X(NX,N),WORK(N)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,K,JM1
C
C     *****SUBROUTINES CALLED:
C     MULWOA
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE COMPUTES THE SYMMETRIC MATRIX PRODUCT
C         T
C        X *S*X WHERE S IS SYMMETRIC AND OVERWRITES S WITH
C     THE RESULT.  BOTH S AND X ARE OF ORDER N.
C
C     *****PARAMETER DESCRIPTION:
C     ON INPUT:
C        NS,NX            ROW DIMENSIONS OF THE ARRAYS CONTAINING S
C                         AND A, RESPECTIVELY, AS DECLARED IN THE
C                         CALLING PROGRAM DIMENSION STATEMENT;
C
C        N                ORDER OF THE MATRICES S AND X;
C
C        S                AN N X N SYMMETRIC MATRIX;
C
C        X                AN N X N MATRIX.
C
C     ON OUTPUT:
C
C                                                             T
C        S                A SYMMETRIC N X N ARRAY CONTAINING X *S*X;
C
C        WORK             A REAL SCRATCH VECTOR OF LENGTH N.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (ELEC. SYS. LAB., M.I.T., RM. 35-331,
C     CAMBRIDGE, MA 02139,  PH.: (617)-253-2125), OCTOBER 1977.
C     MOST RECENT VERSION:  OCT. 12, 1977.
C
C     ------------------------------------------------------------------
C
C     COMPUTE S*X, OVERWRITING INTO S
C
      CALL MULWOA (NS,NX,N,S,X,WORK)
C
C                                    T
C     COMPUTE THE LOWER TRIANGLE OF X *S*X
C
      DO 50 J=1,N
         DO 10 I=J,N
            WORK(I)=0.0D0
10       CONTINUE
         DO 30 K=1,N
            DO 20 I=J,N
               WORK(I)=WORK(I)+X(K,I)*S(K,J)
20          CONTINUE
30       CONTINUE
         DO 40 I=J,N
            S(I,J)=WORK(I)
40       CONTINUE
50    CONTINUE
      IF (N.EQ.1) RETURN
C
C     DETERMINE THE STRICT UPPER TRIANGLE BY SYMMETRY
C
      DO 70 J=2,N
         JM1=J-1
         DO 60 I=1,JM1
            S(I,J)=S(J,I)
60       CONTINUE
70    CONTINUE
      RETURN
C
C     LAST LINE OF MQFWO
C
      END
