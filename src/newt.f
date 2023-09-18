      SUBROUTINE NEWT(NR,NRR,NRW,NRX,N,M,E,A,B,CQC,R,S,RI,RSDM,X,W1,
     X                W2,W3,ALFR,ALFI,WK1,IPVT,RTOL,EFLAG,RFLAG,
     X                RDFLG,SFLAG,SEP,CTIME,MAXIT,INFO,NOUT)
C
C     *****PARAMETERS:
      INTEGER NR,NRR,NRW,NRX,N,M,IPVT(N),MAXIT,INFO,NOUT
      CHARACTER EFLAG,RFLAG,RDFLG,SFLAG
      DOUBLE PRECISION E(NR,N),A(NR,N),B(NR,M),CQC(NR,N),R(NR,M),
     X                 S(NR,M),RI(NR,M),RSDM(NRR,N),X(NRX,N),
     X                 W1(NRW,N),W2(NRW,N),W3(NRW,N),ALFR(N),
     X                 ALFI(N),WK1(N),RTOL,SEP
      LOGICAL CTIME
C
C     *****LOCAL VARIABLES:
      INTEGER I,IER1,IER2,ITER,J
      DOUBLE PRECISION A1N,B1N,BT1N,DPN1,EPS,E1N,H1N,RI1N,TOL,
     X    T1,T2,T3,X1N
C
C     *****FORTRAN FUNCTIONS:
      DOUBLE PRECISION DMAX1,DSQRT
C
C     *****FUNCTION SUBPROGRAMS:
      DOUBLE PRECISION D1NRM
C
C     *****SUBROUTINES CALLED:
C     FBGAIN, LYPCND, LYPDSD, MADD, MMUL, MQF, MQFA, MSCALE, MSUB,
C     SAVE, SEPEST, TRNATA, TRNATB
C
C     --------------------------------------------------------------
C
C     *****PURPOSE:
C     GIVEN THE MODEL MATRICES FOR THE CONTINUOUS- OR DISCRETE-TIME
C     OPTIMAL REGULATOR PROBLEM THAT RESULTS IN A GENERALIZED
C     ALGEBRAIC RICCATI EQUATION (GARE), AND AN INITIAL GUESS FOR
C     THE SOLUTION TO THE GARE, THIS SUBROUTINE APPLIES A NEWTON
C     TYPE ITERATIVE REFINEMENT PROCEDURE.  TO GUARANTEE
C     CONVERGENCE, THE INITIAL GUESS MUST STABILIZE THE CLOSED-LOOP
C     SYSTEM MATRIX E**-1*(A - B*K), WHERE
C
C          K = (R**-1)*(BT*X*E + ST)              CONTINUOUS
C
C          K = ((R + BT*X*B)**-1)*(BT*X*A + ST)   DISCRETE
C
C     AND T DENOTES MATRIX TRANSPOSE.
C
C     REF.:  ARNOLD, W.F., "ON THE NUMERICAL SOLUTION OF
C     ALGEBRAIC MATRIX RICCATI EQUATIONS," PHD THESIS, USC,
C     DECEMBER 1983.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C       NR      INTEGER
C               ROW DIMENSION OF THE ARRAYS CONTAINING THE MATRICES
C               E, A, B, CQC, R, S AND RI AS DECLARED IN THE MAIN
C               CALLING PROGRAM DIMENSION STATEMENT;
C
C       NRR     INTEGER
C               ROW DIMENSION OF THE ARRAY CONTAINING THE MATRIX
C               RSDM AS DECLARED IN THE MAIN CALLING PROGRAM
C               DIMENSION STATEMENT;
C
C       NRW     INTEGER
C               ROW DIMENSION OF THE ARRAYS CONTAINING THE MATRICES
C               W1, W2 AND W3 AS DECLARED IN THE MAIN CALLING
C               PROGRAM DIMENSION STATEMENT;
C
C       NRX     INTEGER
C               ROW DIMENSION OF THE ARRAY CONTAINING THE MATRIX X
C               AS DECLARED IN THE MAIN CALLING PROGRAM DIMENSION
C               STATEMENT;
C
C       N       INTEGER
C               ORDER OF THE SQUARE MATRICES E, A, CQC, RSDM AND X
C               ROW DIMENSION OF THE MATRICES B AND S;
C
C       M       INTEGER
C               ORDER OF THE SQUARE MATRICES R AND RI
C               COLUMN DIMENSION OF THE MATRICES B AND S;
C
C       E       REAL(NR,N)
C               MODEL DESCRIPTOR MATRIX;
C
C       A       REAL(NR,N)
C               MODEL SYSTEM MATRIX;
C
C       B       REAL(NR,M)
C               MODEL INPUT MATRIX;
C
C       CQC     REAL(NR,N)
C               MATRIX PRODUCT CT*Q*C WHERE T DENOTES MATRIX
C               TRANSPOSE;
C
C       R       REAL(NR,M)
C               CONTROL WEIGHTING MATRIX;
C
C       S       REAL(NR,M)
C               STATE - INPUT CROSS-WEIGHTING MATRIX;
C
C       RI      REAL(NR,M)
C               INVERSE OF THE CONTROL WEIGHTING MATRIX;
C
C       X       REAL(NRX,N)
C               INITIAL GUESS FOR RICCATI SOLUTION THAT MUST
C               STABILIZE THE CLOSED-LOOP SYSTEM MATRIX;
C
C       W1,W2,W3
C               REAL(NRW,N)
C               SCRATCH ARRAYS OF SIZE AT LEAST N BY N;
C
C       WK1     REAL(N)
C               WORK VECTOR OF LENGTH AT LEAST N;
C
C       IPVT    INTEGER(N)
C               WORK VECTOR OF LENGTH AT LEAST N;
C
C       RTOL    REAL
C               TOLERANCE ON THE CONDITION ESTIMATE OF R+BT*X*B
C               WITH RESPECT TO INVERSION (DISCRETE PROBLEM).
C               ERROR RETURN IF THIS TOLERANCE IS EXCEEDED;
C
C       EFLAG   CHARACTER
C               FLAG SET TO 'Y' IF E IS OTHER THAN THE IDENTITY
C               MATRIX;
C
C       RFLAG   CHARACTER
C               FLAG SET TO 'Y' IF R IS OTHER THAN THE IDENTITY
C               MATRIX;
C
C       RDFLG   CHARACTER
C               FLAG SET TO 'Y' IF R IS A DIAGONAL MATRIX;
C
C       SFLAG   CHARACTER
C               FLAG SET TO 'Y' IF S IS OTHER THAN THE ZERO MATRIX;
C
C       CTIME    LOGICAL
C               = .TRUE.  FOR CONTINUOUS-TIME SYSTEM
C               = .FALSE. FOR DISCRETE-TIME SYSTEM;
C
C       MAXIT   INTEGER
C               MAXIMUM NUMBER OF NEWTON ITERATIONS PERMITTED;
C
C       NOUT    INTEGER
C               UNIT NUMBER OF OUTPUT DEVICE FOR ERROR WARNING
C               MESSAGES.
C
C     ON OUTPUT:
C
C       RSDM    REAL(NRR,N)
C               DIFFERENCE MATRIX BETWEEN SOLUTIONS FROM LAST TWO
C               ITERATIONS;
C
C       X       REFINED RICCATI SOLUTION MATRIX;
C
C       ALFR    REAL(N)
C               REAL PARTS OF THE CLOSED-LOOP SYSTEM EIGENVALUES
C               UNDER THE OPTIMAL FEEDBACK FOR THE RICCATI SOLUTION
C               AT THE LAST ITERATION;
C
C       ALFI    REAL(N)
C               IMAGINARY PARTS OF THE CLOSED-LOOP SYSTEM
C               EIGENVALUES UNDER THE OPTIMAL FEEDBACK FOR THE
C               RICCATI SOLUTION AT THE LAST ITERATION;
C
C       WK1(1)  1-NORM OF THE RSDM MATRIX;
C
C       IPVT(1) CONVERGENCE CRITERION INDICATOR;
C
C       SEP     REAL
C               ESTIMATE OF THE SEPARATION OF THE CLOSED-LOOP
C               SPECTRUM AT THE LAST ITERATION (CONTINUOUS PROBLEM);
C
C       MAXIT   NUMBER OF ITERATIONS PERFORMED;
C
C       INFO    INTEGER
C               ERROR FLAG WITH THE FOLLOWING MEANINGS
C               = -1  NO CONVERGENCE
C               = -2  INITIAL GUESS NOT STABILIZING
C               = -3  E MATRIX NOT IDENTITY
C               = -4  INDICATES A FAILURE OF THE QR ALGORITHM TO
C                     DETERMINE THE EIGENVALUES IN SOLVING THE
C                     LYAPUNOV EQUATION
C               = -5  CONDITION ESTIMATE OF R+BT*X*B WITH RESPECT TO
C                     INVERSION EXCEEDS THE TOLERANCE VALUE RTOL.
C
C     *****ALGORITHM NOTES:
C     THE ALGORITHM CURRENTLY EMPLOYED IS BASED ON THE BARTELS-
C     STEWART ALGORITHM FOR LYAPUNOV EQUATIONS AND IS VALID FOR THE
C     CASE E = IDENTITY ONLY AT THIS TIME.
C
C     *****HISTORY:
C     THIS SUBROUTINE WAS WRITTEN BY W.F. ARNOLD, NAVAL WEAPONS
C     CENTER, CODE 35104, CHINA LAKE, CA  93555, AS PART OF THE
C     SOFTWARE PACKAGE RICPACK, SEPTEMBER 1983.
C     MODIFIED BY ALAN J. LAUB (UCSB): 01/06/85
C
C     --------------------------------------------------------------
C
      IF(EFLAG .EQ. 'Y' .OR. EFLAG .EQ. 'y') GO TO 200
      INFO = 0
      IER1 = 0
      IER2 = 0
      DPN1 = -1.0D0
C
C     CALCULATE MACHINE PRECISION (EPS)
C
      EPS = 1.0D0
   10 CONTINUE
      EPS = EPS/2.0D0
      IF(EPS+1.0D0 .GT. 1.0D0) GO TO 10
      EPS = EPS * 2.0D0
C
C     START ITERATION LOOP
C
      A1N = D1NRM(NR,N,N,A)
      B1N = D1NRM(NR,N,M,B)
      CALL TRNATB(NR,NRW,N,M,B,W1)
      BT1N = D1NRM(NRW,M,N,W1)
      RI1N = D1NRM(NR,M,M,RI)
      H1N = D1NRM(NR,N,N,CQC)
      DO 180 ITER=1,MAXIT
C
C     CALCULATE FEEDBACK GAIN K AND STORE IT IN W1
C
      CALL FBGAIN(NR,NRX,NRW,N,M,A,B,E,R,RI,S,X,W1,W2,WK1,IPVT,EFLAG,
     X            RDFLG,RFLAG,SFLAG,CTIME)
      IF(WK1(1) .GT. RTOL) GO TO 220
C
C     SET UP LYAPUNOV EQUATION
C
      CALL MMUL(NR,NRW,NRW,N,N,M,B,W1,W2)
      CALL MSUB(NR,NRW,NRW,N,N,A,W2,W2)
      CALL SAVE(NRX,NRR,N,N,X,RSDM)
      IF (SFLAG.EQ.'Y' .OR. SFLAG.EQ.'y')
     X               CALL MMUL(NR,NRW,NRW,N,N,M,S,W1,W3)
      IF (RFLAG.EQ.'Y' .OR. RFLAG.EQ.'y') GO TO 30
   20 CONTINUE
      CALL MQFA(NRW,NRX,M,N,W1,X,WK1)
      GO TO 70
   30 CONTINUE
      IF (RDFLG .EQ. 'Y' .OR. RDFLG .EQ. 'y') GO TO 40
      CALL MQF(NR,NRW,NRX,M,N,R,W1,X,WK1)
      GO TO 70
   40 CONTINUE
      DO 60 I=1,M
         T1 = DSQRT(R(I,I))
         DO 50 J=1,N
            W1(I,J) = T1*W1(I,J)
   50    CONTINUE
   60 CONTINUE
      GO TO 20
   70 CONTINUE
      T1 = D1NRM(NRW,N,N,W1)
      CALL MADD(NRX,NR,NRX,N,N,X,CQC,X)
      IF (SFLAG .NE. 'Y' .AND. SFLAG .NE. 'y') GO TO 80
      CALL MSUB(NRX,NRW,NRX,N,N,X,W3,X)
      CALL TRNATA(NRW,N,W3)
      CALL MSUB(NRX,NRW,NRX,N,N,X,W3,X)
   80 CONTINUE
      IF (.NOT. CTIME) CALL MSCALE(NRX,N,N,DPN1,X)
C
C     SOLVE LYAPUNOV EQUATION
C
      IF (CTIME) GO TO 100
C
C     SAVE SOME PARAMETERS FOR USE IN CONVERGENCE TESTS
C
      T1 = D1NRM(NRW,N,N,W2)
      CALL TRNATB(NRW,NRW,N,N,W2,W3)
      T2 = D1NRM(NRW,N,N,W3)
      T3 = 1.0D0 - T1*T2
C
C     DISCRETE SOLUTION
C
      CALL LYPDSD(NRW,NRX,N,W2,X,W1,ALFR,ALFI,WK1,W3,IPVT,IER1,IER2)
      IF(IER1 .NE. 0) GO TO 210
C
C     CHECK FOR STABILIZING SOLUTION ON FIRST ITERATION
C
      IF(ITER .GT. 1) GO TO 120
      DO 90 I=1,N
          IF(DSQRT(ALFR(I)*ALFR(I)+ALFI(I)*ALFI(I)) .GT. 1.0D0)
     X    GO TO 190
   90 CONTINUE
      GO TO 120
  100 CONTINUE
C
C     CONTINUOUS SOLUTION
C
      CALL LYPCND(NRW,NRX,N,W2,X,W1,ALFR,ALFI,WK1,IER1,IER2)
      IF(IER1 .NE. 0) GO TO 210
C
C     CHECK FOR STABILIZING SOLUTION ON FIRST ITERATION
C
      IF(ITER .GT. 1) GO TO 120
      DO 110 I=1,N
          IF(ALFR(I) .GE. 0.0D0) GO TO 190
  110 CONTINUE
C
C     CONVERGENCE CHECKS
C
  120 CONTINUE
C
C     MAIN CHECK, CONTINUOUS OR DISCRETE
C
      CALL MSUB(NRX,NRR,NRR,N,N,X,RSDM,RSDM)
      E1N = D1NRM(NRR,N,N,RSDM)
      WK1(1) = E1N
      X1N = D1NRM(NRX,N,N,X)
      IF(E1N/X1N .LE. 10.0D0*N*EPS) GO TO 230
      IF (CTIME) GO TO 150
C
C     SPECIAL CHECKS - DISCRETE
C
      IF (T3 .LE. 0.0D0) GO TO 130
C
C     DISCRETE CHECK 2
C
      IF (E1N .LE. EPS*H1N/T3) GO TO 240
C
C     DISCRETE CHECK 3
C
      TOL = EPS*(H1N/X1N+2.0D0*T1*T2)/T3
      IF (E1N/X1N .LE. TOL) GO TO 250
      GO TO 170
C
C     DISCRETE CHECK 4
C
  130 CONTINUE
      DO 140 I=1,N
          T3 = DMAX1(T3, ALFR(I)*ALFR(I)+ALFI(I)*ALFI(I))
  140 CONTINUE
      T3 = 1.0D0 - T3
      IF (E1N .LE. EPS*H1N/T3) GO TO 260
C
C     DISCRETE CHECK 5
C
      TOL = EPS*(H1N/X1N+2.0D0*T1*T2)/T3
      IF (E1N/X1N .LE. TOL) GO TO 270
      GO TO 170
C
C     SPECIAL CHECKS - CONTINUOUS
C
  150 CONTINUE
C
C     CONTINUOUS CHECK 2
C
      CALL SEPEST(NRW,N,W2,W1,SEP,INFO)
      IF(INFO .EQ. 0) GO TO 160
C      WRITE(NOUT,900) ITER
C  900 FORMAT(//' UNSUCCESSFUL ATTEMPT TO ESTIMATE SEPARATION ON',
C    X       I3,' -TH ITERATION OF NEWTON''S METHOD'//)
C      GO TO 170
  160 CONTINUE
      IF(E1N .LE. EPS*H1N/SEP) GO TO 240
C
C     CONTINUOUS CHECK 3
C
      TOL = EPS*(H1N+2.0D0*A1N*X1N+T1)
      IF (E1N*SEP .LE. TOL) GO TO 250
C
C     CONTINUOUS CHECK 4
      TOL = EPS*(H1N+2.0D0*A1N*X1N+B1N*RI1N*BT1N*X1N*X1N)
      IF (E1N*SEP .LE. TOL) GO TO 260
  170 CONTINUE
  180 CONTINUE
C
C     END OF ITERATION LOOP
C
      IF(INFO .NE. 0) RETURN
      INFO = -1
      RETURN
  190 CONTINUE
      INFO = -2
      RETURN
  200 CONTINUE
      INFO = -3
      RETURN
  210 CONTINUE
      INFO = -4
      RETURN
  220 CONTINUE
      INFO = -5
      RETURN
  230 CONTINUE
      IPVT(1) = 1
      MAXIT = ITER
      RETURN
  240 CONTINUE
      IPVT(1) = 2
      MAXIT = ITER
      RETURN
  250 CONTINUE
      IPVT(1) = 3
      MAXIT = ITER
      RETURN
  260 CONTINUE
      IPVT(1) = 4
      MAXIT = ITER
      RETURN
  270 CONTINUE
      IPVT(1) = 5
      MAXIT = ITER
      RETURN
      END