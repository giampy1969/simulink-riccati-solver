      SUBROUTINE SYMSLV (NA,NC,N,A,C)
C
C     *****PARAMETERS:
      INTEGER NA,NC,N
      DOUBLE PRECISION A(NA,N),C(NC,N)
C
C     *****LOCAL VARIABLES:
      INTEGER K,KK,DK,KM1,L,LL,DL,LDL,I,IA,J,NSYS,NT,IPVT(4)
      DOUBLE PRECISION ALFA,T(4,4),P(4),COND,WORK(4),S,EPS,EPSP1
C
C     *****SUBROUTINES CALLED:
C     LINEQ(DGECOM,DGESLM),MSCALE
C
C     ------------------------------------------------------------------
C
C     ******PURPOSE:
C     THIS DOUBLE PRECISION SUBROUTINE SOLVES THE REAL MATRIX EQUATION
C      T
C     A *X + X*A + C = 0, WHERE C IS SYMMETRIC AND A IS IN UPPER
C     REAL SCHUR FORM.
C
C     *****PARAMETER DESCRIPTION:
C     ON INPUT:
C        NA,NC            ROW DIMENSIONS OF THE ARRAYS CONTAINING A AND
C                         C AS DECLARED IN THE CALLING PROGRAM DIMENSION
C                         STATEMENT;
C
C        N                ORDER OF THE MATRICES A AND C;
C
C        A                AN N X N (REAL) MATRIX IN UPPER SCHUR FORM;
C
C        C                AN N X N (REAL) ARRAY.
C
C     ON OUTPUT:
C
C        C                AN N X N (REAL) MATRIX CONTAINING THE
C                         SOLUTION;
C
C     *****HISTORY:
C     THIS IS A MODIFICATION OF THE BARTELS-STEWART SUBROUTINE SYMSLV
C     WRITTEN BY ALAN J. LAUB (DEP'T. OF EE-SYSTEMS, UNIV. OF SOUTHERN
C     CALIF., LOS ANGELES, CA 90089, PH.: (213) 743-5535) SEP. 1977
C     MOST RECENT VERSION: JUN. 29, 1982.
C
C     ------------------------------------------------------------------
C
C     EPS IS AN INTERNALLY GENERATED MACHINE DEPENDENT PARAMETER
C     SPECIFYING THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C     FOR EXAMPLE, EPS = 16.0D0**(-13) FOR DOUBLE PRECISION ARITHMETIC
C     ON IBM S360/S370.
C
      EPS=1.0D0
    1 EPS=EPS/2.0D0
      EPSP1=EPS+1.0D0
      IF(EPSP1.GT.1.0D0) GO TO 1
      EPS=2.0D0*EPS
C
      ALFA = -1.0D0
      CALL MSCALE (NC,N,N,ALFA,C)
      NT = 4
      L = 1
   10 DL = 1
      IF(L.EQ.N) GO TO 20
      S=DABS(A(L,L))+DABS(A(L+1,L+1))
      IF(DABS(A(L+1,L)).GT.EPS*S) DL=2
   20 LL = L+DL-1
      K = L
   30 KM1 = K-1
      DK = 1
      IF(K.EQ.N) GO TO 35
      S=DABS(A(K,K))+DABS(A(K+1,K+1))
      IF(DABS(A(K+1,K)).GT.EPS*S) DK=2
   35 KK= K+DK-1
      IF(K.EQ.L) GO TO 45
      DO 40 I=K,KK
         DO 40 J=L,LL
            DO 40 IA=L,KM1
               C(I,J) = C(I,J) - A(IA,I)*C(IA,J)
   40 CONTINUE
   45 IF(DL.EQ.2) GO TO 60
      IF(DK.EQ.2) GO TO 50
      T(1,1) = A(K,K) + A(L,L)
      IF (T(1,1).EQ.0.0D0) RETURN
      C(K,L) = C(K,L)/T(1,1)
      GO TO 90
   50 T(1,1) = A(K,K) + A(L,L)
      T(1,2) = A(KK,K)
      T(2,1) = A(K,KK)
      T(2,2) = A(KK,KK) + A(L,L)
      P(1) = C(K,L)
      P(2) = C(KK,L)
      NSYS = 2
      CALL LINEQ (NT,NSYS,T,P,COND,IPVT,WORK)
      C(K,L) = P(1)
      C(KK,L) = P(2)
      GO TO 90
   60 IF(DK.EQ.2) GO TO 70
      T(1,1) = A(K,K) + A(L,L)
      T(1,2) = A(LL,L)
      T(2,1) = A(L,LL)
      T(2,2) = A(K,K) + A(LL,LL)
      P(1) = C(K,L)
      P(2) = C(K,LL)
      NSYS= 2
      CALL LINEQ (NT,NSYS,T,P,COND,IPVT,WORK)
      C(K,L) = P(1)
      C(K,LL) = P(2)
      GO TO 90
   70 IF(K.NE.L) GO TO 80
      T(1,1) = A(L,L)
      T(1,2) = A(LL,L)
      T(1,3) = 0.D0
      T(2,1) = A(L,LL)
      T(2,2) = A(L,L) + A(LL,LL)
      T(2,3) = T(1,2)
      T(3,1) = 0.D0
      T(3,2) = T(2,1)
      T(3,3) = A(LL,LL)
      P(1) = C(L,L)/2.D0
      P(2) =  C(LL,L)
      P(3) = C(LL,LL)/2.D0
      NSYS = 3
      CALL LINEQ (NT,NSYS,T,P,COND,IPVT,WORK)
      C(LL,L) = P(2)
      C(L,L) = P(1)
      C(L,LL) = P(2)
      C(LL,LL) = P(3)
      GO TO 90
   80 T(1,1) = A(K,K) + A(L,L)
      T(1,2) = A(KK,K)
      T(1,3) = A(LL,L)
      T(1,4) = 0.D0
      T(2,1) = A(K,KK)
      T(2,2) = A(KK,KK) + A(L,L)
      T(2,3) = 0.D0
      T(2,4) = T(1,3)
      T(3,1) = A(L,LL)
      T(3,2) = 0.D0
      T(3,3) = A(K,K) + A(LL,LL)
      T(3,4) = T(1,2)
      T(4,1) = 0.D0
      T(4,2) = T(3,1)
      T(4,3) = T(2,1)
      T(4,4) = A(KK,KK) + A(LL,LL)
      P(1) = C(K,L)
      P(2) = C(KK,L)
      P(3) = C(K,LL)
      P(4) = C(KK,LL)
      NSYS = 4
      CALL LINEQ (NT,NSYS,T,P,COND,IPVT,WORK)
      C(K,L) = P(1)
      C(KK,L) = P(2)
      C(K,LL) = P(3)
      C(KK,LL) = P(4)
   90 K = K + DK
      IF(K.LE.N) GO TO 30
      LDL = L + DL
      IF(LDL.GT.N) RETURN
      DO 120 J=LDL,N
         DO 100 I=L,LL
            C(I,J) = C(J,I)
  100    CONTINUE
         DO 120 I=J,N
            DO 110 K=L,LL
               C(I,J) = C(I,J) - C(I,K)*A(K,J) - A(K,I)*C(K,J)
  110       CONTINUE
            C(J,I) = C(I,J)
  120 CONTINUE
      L = LDL
      GO TO 10
C
C     LAST LINE OF SYMSLV
C
      END
