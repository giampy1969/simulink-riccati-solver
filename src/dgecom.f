      SUBROUTINE DGECOM (A,LDA,N,IPVT,RCOND,Z)
      INTEGER LDA,N,IPVT(N)
      DOUBLE PRECISION A(LDA,N),Z(N)
      DOUBLE PRECISION RCOND
C
C     DGECOM FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF RCOND IS NOT NEEDED, DGEFAM IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DGECOM BY DGESLM.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECOM BY DGESLM.
C
C     ON ENTRY:
C
C        A       DOUBLE PRECISION(LDA,N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY A AS DECLARED
C                IN THE MAIN CALLING PROGRAM.
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX A.
C
C     ON RETURN:
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND U IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C                   IPVT(K) = THE INDEX OF THE K-TH PIVOT ROW
C                THE DETERMINANT OF A CAN BE OBTAINED ON OUTPUT BY
C                   DET(A) = S*A(1,1)*A(2,2)* ... *A(N,N)
C                WHERE S = (-1)**(NUMBER OF TIMES IPVT(K) .NE. K)
C                BUT THIS RESULT SHOULD BE USED WITH CAUTION AS IT
C                MAY EASILY UNDERFLOW OR OVERFLOW.
C
C        RCOND   DOUBLE PRECISION
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF A.
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN A AND B OF SIZE EPSILON MAY CAUSE RELATIVE
C                PERTURBATIONS IN X OF SIZE  EPSILON/RCOND.
C                IF RCOND IS SO SMALL THAT THE LOGICAL EXPRESSION
C                             1.0D0 + RCOND  .EQ.  1.0D0
C                IS TRUE, THEN A MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR, RCOND IS ZERO IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       DOUBLE PRECISION(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF A IS CLOSE TO A SINGULAR MATRIX, THEN Z IS AN
C                APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                       NORM(A*Z) = RCOND*NORM(A)*NORM(Z).
C
C     THIS VERSION IS ADAPTED FROM THE LINPACK SUBROUTINE DGECO BY
C     REPLACEMENT OF CALLS TO THE BLAS WITH IN-LINE CODE.  MODIFICATION
C     DONE BY ALAN J. LAUB, UNIVERSITY OF SOUTHERN CALIFORNIA,
C     APRIL 1980.
C
C     SUBROUTINES AND FUNCTIONS CALLED:
C
C     DGEFAM
C
C     FORTRAN FUNCTIONS CALLED:
C
      DOUBLE PRECISION DABS,DSIGN
C
C     INTERNAL VARIABLES:
C
      INTEGER I,INFO,J,K,KB,KM1,KP1,L
      DOUBLE PRECISION ANORM,EK,S,SM,T,WK,WKM,YNORM
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0D0
      DO 20 J=1,N
         T = 0.0D0
         DO 10 I=1,N
            T = T+DABS(A(I,J))
   10    CONTINUE
         IF (T .GT. ANORM) ANORM = T
   20 CONTINUE
C
C     FACTOR
C
      CALL DGEFAM(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))).
C     ESTIMATE = NORM(Z)/NORM(Y)  WHERE  A*Z = Y AND  TRANS(A)*Y = E.
C     TRANS(A)  IS THE TRANSPOSE OF A.  THE COMPONENTS OF E ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE  TRANS(U)*W = E
C
      EK = 1.0D0
      DO 30 I=1,N
         Z(I) = 0.0D0
   30 CONTINUE
      DO 120 K=1,N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 50
         S = DABS(A(K,K))/DABS(EK-Z(K))
         DO 40 I=1,N
            Z(I) = S*Z(I)
   40    CONTINUE
         EK = S*EK
   50    CONTINUE
         WK = EK-Z(K)
         WKM = -EK-Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 60
         WK = WK/A(K,K)
         WKM = WKM/A(K,K)
         GO TO 70
   60    CONTINUE
         WK = 1.0D0
         WKM = 1.0D0
   70    CONTINUE
         KP1 = K+1
         IF (KP1 .GT. N) GO TO 110
         DO 80 J=KP1,N
            SM = SM+DABS(Z(J)+WKM*A(K,J))
            Z(J) = Z(J)+WK*A(K,J)
            S = S+DABS(Z(J))
   80    CONTINUE
         IF (S .GE. SM) GO TO 100
         T = WKM-WK
         WK = WKM
         DO 90 J=KP1,N
            Z(J) = Z(J)+T*A(K,J)
   90    CONTINUE
  100    CONTINUE
  110    CONTINUE
         Z(K) = WK
  120 CONTINUE
      T = 0.0D0
      DO 130 I=1,N
         T = T+DABS(Z(I))
  130 CONTINUE
      S = 1.0D0/T
      DO 140 I=1,N
         Z(I) = S*Z(I)
  140 CONTINUE
C
C     SOLVE  TRANS(L)*Y = W
C
      DO 190 KB=1,N
         K = N+1-KB
         IF (K .EQ. N) GO TO 160
         KP1 = K+1
         T = 0.0D0
         DO 150 I=KP1,N
            T = T+A(I,K)*Z(I)
  150    CONTINUE
         Z(K) = Z(K)-T
  160    CONTINUE
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 180
         S = 1.0D0/DABS(Z(K))
         DO 170 I=1,N
            Z(I) = S*Z(I)
  170    CONTINUE
  180    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  190 CONTINUE
      T = 0.0D0
      DO 200 I=1,N
         T = T+DABS(Z(I))
  200 CONTINUE
      S = 1.0D0/T
      DO 210 I=1,N
         Z(I) = S*Z(I)
  210 CONTINUE
C
      YNORM = 1.0D0
C
C     SOLVE  L*V = Y
C
      DO 260 K=1,N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .EQ. N) GO TO 230
         KP1 = K+1
         DO 220 I=KP1,N
            Z(I) = Z(I)-T*A(I,K)
  220    CONTINUE
  230    CONTINUE
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 250
         S = 1.0D0/DABS(Z(K))
         DO 240 I=1,N
            Z(I) = S*Z(I)
  240    CONTINUE
         YNORM = S*YNORM
  250    CONTINUE
  260 CONTINUE
      T = 0.0D0
      DO 270 I=1,N
         T = T+DABS(Z(I))
  270 CONTINUE
      S = 1.0D0/T
      DO 280 I=1,N
         Z(I) = S*Z(I)
  280 CONTINUE
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 330 KB=1,N
         K = N+1-KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 300
         S = DABS(A(K,K))/DABS(Z(K))
         DO 290 I=1,N
            Z(I) = S*Z(I)
  290    CONTINUE
         YNORM = S*YNORM
  300    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         IF (K .EQ. 1) GO TO 320
         KM1 = K-1
         DO 310 I=1,KM1
            Z(I) = Z(I)+T*A(I,K)
  310    CONTINUE
  320    CONTINUE
  330 CONTINUE
C
C     MAKE ZNORM = 1.0D0
C
      T = 0.0D0
      DO 340 I=1,N
         T = T+DABS(Z(I))
  340 CONTINUE
      S = 1.0D0/T
      DO 350 I=1,N
         Z(I) = S*Z(I)
  350 CONTINUE
      YNORM = S*YNORM
C
C     DETERMINE RCOND
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
