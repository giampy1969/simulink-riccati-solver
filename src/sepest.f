      SUBROUTINE SEPEST(NR,N,T,Q,SEP,INFO)
C
C     *****PARAMETERS:
      INTEGER NR,N,INFO
      DOUBLE PRECISION T(NR,N),Q(NR,N),SEP
C
C     *****LOCAL VARIABLES:
      INTEGER IPVT(4),I,II,IM1,J,JP1,K1,K2,K1M1,L1,L2,L1M1,L1MK1,
     X        M,ND,NM1
      DOUBLE PRECISION A(4,4),VEC(4),Z(4),A1N,RCOND,S,TEMP,T1N
C
C     *****FORTRAN FUNCTIONS:
      DOUBLE PRECISION DABS,DMAX1,DMIN1,DSIGN
C
C     *****FUNCTION SUBPROGRAMS:
      DOUBLE PRECISION D1NRM
C
C     *****SUBROUTINES CALLED:
C     DGECOM, DGESLM, MSCALE, SYMSLV
C
C     --------------------------------------------------------------
C
C     *****PURPOSE:
C     GIVEN A QUASITRIANGULAR MATRIX T AND A SYMMETRIC MATRIX Q,
C     THIS SUBROUTINE COMPUTES
C
C       SEP = MIN( 1-NORM( T-TRANSPOSE*Q + Q*T)/ 1-NORM(Q))
C
C     REF.:  BARTELS, R.H. AND G.W. STEWART, "SOLUTION OF THE MATRIX
C     EQUATION A*X + X*B = C," COMM. OF THE ACM, VOL. 15,
C     PP. 820-826, 1972.
C     CLINE, A.K., MOLER, C.B., STEWART, G.W. AND J.H. WILKINSON,
C     "AN ESTIMATE OF THE CONDITION NUMBER OF A MATRIX," SIAM J. OF
C     NUMERICAL ANALYSIS, VOL. 16, PP. 368-375, 1979.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C       NR      INTEGER
C               ROW DIMENSION OF THE ARRAYS CONTAINING THE MATRICES
C               T AND Q AS DECLARED IN THE MAIN CALLING PROGRAM
C               DIMENSION STATEMENT;
C
C       N       INTEGER
C               ORDER OF THE SQUARE MATRICES Q AND T;
C
C       T       REAL(NR,N)
C               QUASITRIANGULAR INPUT MATRIX.
C
C     ON OUTPUT:
C
C       SEP     REAL
C               AN ESTIMATE OF THE QUANTITY SPECIFIED ABOVE;
C
C       Q       REAL(NR,N)
C               THE MATRIX THAT MINIMIZES THE QUANTITY SEP;
C
C       INFO    INTEGER
C               ERROR FLAG WITH THE FOLLOWING MEANINGS
C               = 0 INDICATES NORMAL RETURN
C               = (N+1)*L+K INDICATES THE L-TH AND K-TH EIGENVALUES
C                 OF T FORM A +/- PAIR, SO SEP IS EQUAL TO ZERO.
C
C     *****ALGORITHM NOTES:
C                           T
C     LET PHI(Y) = T*Y + Y*T .  Q IS OBTAINED BY INVERSE ITERATION ON
C            T
C     PHI*PHI .  THE STARTING VALUE OF Q IS CHOSEN AS FOLLOWS:
C     PARTITION ALL MATRICES CONFORMALLY WITH T.  Q IS CHOSEN TO
C     SATISFY
C                 T           T
C        T*Y + Y*T  = B,     T *Q + Q*T = Y/1-NORM(Y).
C
C     *****HISTORY:
C     THIS SUBROUTINE IS A MODIFIED VERSION OF THE SUBROUTINE OF THE
C     SAME NAME WRITTEN BY RALPH BYERS, 2/82 REF.:  BYERS, R.,
C     "HAMILTONIAN AND SYMPLECTIC ALGORITHMS FOR THE ALGEBRAIC
C     RICCATI EQUATION," PHD THESIS, CORNELL UNIVERSITY, PP.289-295,
C     JANUARY 1983.  THE MODIFICATIONS WERE MADE BY W.F. ARNOLD,
C     NAVAL WEAPONS CENTER, CODE 35104, CHINA LAKE, CA 93555, AS
C     PART OF THE SOFTWARE PACKAGE RICPACK, SEPTEMBER 1983.
C
C     --------------------------------------------------------------
C
      INFO = 0
C
C     INITIALIZATION
C
      IF(N .NE. 1) GO TO 5
      SEP = 2.0D0*DABS(T(1,1))
      RETURN
    5 CONTINUE
      T1N = 0.0D0
      DO 10 I=1,N
          T1N = T1N + DABS(T(I,N))
          Q(I,N) = 0.0D0
   10 CONTINUE
      NM1 = N - 1
      DO 40 J=1,NM1
          JP1 = J + 1
          TEMP = 0.0D0
          DO 20 I=1,JP1
              TEMP = TEMP+DABS(T(I,J))
   20     CONTINUE
          DO 30 I=1,N
              Q(I,J) = 0.0D0
   30     CONTINUE
          T1N = DMAX1(T1N,TEMP)
   40 CONTINUE
      SEP = 2.0D0*T1N
C
C     START COLUMN LOOP (INDEX = L)
C
      L1 = N + 1
   50 CONTINUE
      L2 = L1 - 1
      IF(L2 .LT. 1) GO TO 300
      L1 = L2
      IF(L2 .EQ. 1) GO TO 60
      IF(T(L2,L2-1) .NE. 0.0D0) L1 = L2 - 1
   60 CONTINUE
      IF(L1 .NE. L2) GO TO 80
C
C     ONE BY ONE CASE
C
      SEP = DMIN1(SEP,2.0D0*DABS(T(L1,L1)))
      IF(T(L1,L1) .EQ. 0.0D0) INFO = (N+1)*L1+L2
      IF(INFO .NE. 0) RETURN
      Q(L1,L1) = (DSIGN(1.0D0,Q(L1,L1))+Q(L1,L1))/(2.0D0*T(L1,L1))
      L1M1 = L1 - 1
      DO 70 I = 1,L1M1
          Q(I,L1) = -Q(L1,L1)*T(I,L1)+Q(I,L1)
   70 CONTINUE
      GO TO 100
   80 CONTINUE
C
C     TWO BY TWO CASE
C
      A(1,1) = T(L1,L1)
      A(1,2) = T(L1,L2)
      A(1,3) = 0.0D0
      A(2,1) = T(L2,L1)
      A(2,2) = T(L1,L1) + T(L2,L2)
      A(2,3) = T(L1,L2)
      A(3,1) = 0.0D0
      A(3,2) = T(L2,L1)
      A(3,3) = T(L2,L2)
      A1N = DMAX1(DABS(A(1,1))+DABS(A(2,1)),
     X            DABS(A(1,2))+DABS(A(2,2))+DABS(A(3,2)),
     X            DABS(A(2,3))+DABS(A(3,3)))
      VEC(1) = Q(L1,L1)/2.0D0
      VEC(2) = Q(L2,L1)
      VEC(3) = Q(L2,L2)/2.0D0
      CALL DGECOM(A,4,3,IPVT,RCOND,Z)
      SEP = DMIN1(SEP,2.0D0*A1N*RCOND)
      IF(1.0D0+RCOND .EQ. 1.0D0) INFO = (N+1)*L1+L2
      IF(INFO .NE. 0) RETURN
      CALL DGESLM(A,4,3,IPVT,VEC)
      CALL DGESLM(A,4,3,IPVT,Z)
      S = -DSIGN(1.0D0,Z(1)*VEC(1)+Z(2)*VEC(2)+Z(3)*VEC(3))
      Q(L1,L1) = VEC(1) - S*Z(1)
      Q(L2,L1) = VEC(2) - S*Z(2)
      Q(L1,L2) = Q(L2,L1)
      Q(L2,L2) = VEC(3) - S*Z(3)
      IF(L1M1 .LT. 1) GO TO 100
      DO 90 I=1,L1M1
          Q(I,L1) = Q(I,L1)-T(I,L1)*Q(L1,L1)-T(I,L2)*Q(L2,L1)
          Q(I,L2) = Q(I,L2)-T(I,L1)*Q(L1,L2)-T(I,L2)*Q(L2,L2)
   90 CONTINUE
  100 CONTINUE
C
C     START ROW LOOP (INDEX = K)
C
      K1 = L1
  110 CONTINUE
      K2 = K1 - 1
      IF(K2 .LT. 1) GO TO 290
      K1 = K2
      IF(K2 .EQ. 1) GO TO 120
      IF(T(K2-1,K2) .NE. 0.0D0) K1 = K2-1
  120 CONTINUE
C
C     CASE K NOT EQUAL TO L
C
      IF(.NOT. (L1.EQ.L2 .AND. K1.EQ.K2)) GO TO 130
C
C     ONE BY ONE CASE
C
      ND = 1
      A(1,1) = T(L1,L1) + T(K1,K1)
      A1N = DABS(A(1,1))
      GO TO 160
  130 CONTINUE
      IF(.NOT. (L1.EQ.L2 .AND. K1.NE.K2)) GO TO 140
C
C     TWO BY ONE CASE
C
      ND = 2
      A(1,1) = T(L1,L1) + T(K1,K1)
      A(1,2) = T(K1,K2)
      A(2,1) = T(K2,K1)
      A(2,2) = T(L2,L2) + T(K2,K2)
      A1N = DMAX1(DABS(A(1,1))+DABS(A(2,1)),
     X            DABS(A(1,2))+DABS(A(2,2)))
      GO TO 160
  140 CONTINUE
      IF(.NOT. (L1.NE.L2 .AND. K1.EQ.K2)) GO TO 150
C
C     ONE BY TWO CASE
C
      ND = 2
      A(1,1) = T(L1,L1) + T(K1,K1)
      A(1,2) = T(L1,L2)
      A(2,1) = T(L2,L1)
      A(2,2) = T(L2,L2) + T(K2,K2)
      A1N = DMAX1(DABS(A(1,1))+DABS(A(2,1)),
     X            DABS(A(1,2))+DABS(A(2,2)))
      GO TO 160
  150 CONTINUE
C
C     TWO BY TWO CASE
C
      ND = 4
      A(1,1) = T(L1,L1) + T(K1,K1)
      A(1,2) = T(K1,K2)
      A(1,3) = T(L1,L2)
      A(1,4) = 0.0D0
      A(2,1) = T(K2,K1)
      A(2,2) = T(K2,K2) + T(L1,L1)
      A(2,3) = 0.0D0
      A(2,4) = T(L1,L2)
      A(3,1) = T(L2,L1)
      A(3,2) = 0.0D0
      A(3,3) = T(K1,K1) + T(L2,L2)
      A(3,4) = T(K1,K2)
      A(4,1) = 0.0D0
      A(4,2) = T(L2,L1)
      A(4,3) = T(K2,K1)
      A(4,4) = T(K2,K2) + T(L2,L2)
      A1N = DMAX1(DABS(A(1,1))+DABS(A(2,1))+DABS(A(3,1)),
     X            DABS(A(1,2))+DABS(A(2,2))+DABS(A(4,2)),
     X            DABS(A(1,3))+DABS(A(3,3))+DABS(A(4,3)),
     X            DABS(A(2,4))+DABS(A(3,4))+DABS(A(4,4)))
  160 CONTINUE
      M = 1
      DO 180 J=L1,L2
          DO 170 I=K1,K2
              VEC(M) = Q(I,J)
              M = M + 1
  170     CONTINUE
  180 CONTINUE
      CALL DGECOM(A,4,ND,IPVT,RCOND,Z)
      SEP = DMIN1(SEP,2.0D0*A1N*RCOND)
      IF(1.0D0+RCOND .EQ. 1.0D0) INFO = (N+1)*L1+K1
      IF(INFO .NE. 0) RETURN
      CALL DGESLM(A,4,ND,IPVT,VEC)
      CALL DGESLM(A,4,ND,IPVT,Z)
      S = 0.0D0
      DO 190 J=1,ND
          S = S + Z(J)*VEC(J)
  190 CONTINUE
      S = -DSIGN(1.0D0,S)
      M = 1
      DO 210 J=L1,L2
          DO 200 I = K1,K2
              VEC(M) = VEC(M)+S*Z(M)
              Q(I,J) = VEC(M)
              M = M + 1
  200     CONTINUE
  210 CONTINUE
      DO 280 J=L1,L2
          DO 270 I=K1,K2
              L1MK1 = L1 - K1
              K1M1 = K1 - 1
              IF(K1M1 .EQ. 0) GO TO 230
              DO 220 II = 1,K1M1
                  Q(II,J) = -Q(I,J)*T(II,I)+Q(II,J)
  220         CONTINUE
  230         CONTINUE
              DO 240 II=1,K2
                  Q(II,I) = -Q(I,J)*T(II,J)+Q(II,I)
  240         CONTINUE
              IF(L1MK1 .EQ. 0) GO TO 260
              M = K1
              DO 250 II=1,L1MK1
                  Q(I,M) = -Q(I,J)*T(K1M1+II,J)+Q(I,M)
                  M = M + 1
  250         CONTINUE
  260         CONTINUE
  270     CONTINUE
  280 CONTINUE
C
C     END OF ROW LOOP
C
      GO TO 110
  290 CONTINUE
C
C     END OF COLUMN LOOP
C
      GO TO 50
  300 CONTINUE
C
C     FILL IN LOWER TRIANGLE OF Y
C
      DO 320 I=2,N
          IM1 = I - 1
          DO 310 J = 1,IM1
              Q(I,J) = Q(J,I)
  310     CONTINUE
  320 CONTINUE
C
C     SOLVE FOR Q
C
      S = -1.0D0/D1NRM(NR,N,N,Q)
      CALL MSCALE(NR,N,N,S,Q)
      CALL SYMSLV(NR,NR,N,T,Q)
      S = 1.0D0/D1NRM(NR,N,N,Q)
      SEP = DMIN1(SEP,S)
      CALL MSCALE(NR,N,N,S,Q)
      RETURN
C
C     LAST LINE OF SEPEST
C
      END
