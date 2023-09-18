      SUBROUTINE DSTSLV (NA,NC,N,A,C,U,IDIM)
C
C     *****PARAMETERS:
      INTEGER NA,NC,N,IDIM(N)
      DOUBLE PRECISION A(NA,N),U(NA,N),C(NC,N)
C
C     *****LOCAL VARIABLES:
      INTEGER K,KK,DK,KM1,L,LL,DL,LDL,I,IA,J,NSYS,DM,IBLK,
     +        NT,IPVT(4),IBM1,IC,JJ,IDM,IDK,IB
      DOUBLE PRECISION T(4,4),P(4),COND,WORK(4),S,EPS,EPSP1
C
C     *****SUBROUTINES CALLED:
C     LINEQ
C
C     ------------------------------------------------------------------
C
C     ******PURPOSE:
C     THIS SUBROUTINE SOLVES THE REAL MATRIX EQUATION
C      T
C     A *X*A - X = C, WHERE C IS SYMMETRIC, A IS IN UPPER
C     REAL SCHUR FORM.
C
C     *****PARAMETER DESCRIPTION:
C     ON INPUT:
C        NA,NC            ROW DIMENSIONS OF THE ARRAYS CONTAINING A
C                         (AND U), AND C, RESPECTIVELY, AS DECLARED IN
C                         THE CALLING PROGRAM DIMENSION STATEMENT;
C
C        N                ORDER OF THE MATRICES A AND C;
C
C        A                N X N (REAL) MATRIX IN UPPER SCHUR FORM;
C
C        C                N X N (REAL) SYMMETRIC MATRIX.
C
C     ON OUTPUT:
C
C        C                N X N (REAL) MATRIX CONTAINING THE
C                         SOLUTION;
C
C        U                N X N (REAL) SCRATCH ARRAY;
C
C        IDIM             INTEGER WORK VECTOR OF LENGTH N.
C
C     *****HISTORY:
C     WRITTEN BY J.A.K. CARRIG (ELEC. SYS. LAB., M.I.T., RM. 35-427,
C     CAMBRIDGE, MA 02139, PH.: (617) - 253-7263, AUGUST 1978.
C     MOST RECENT VERSION: AUG. 29, 1978.
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
      DO 5 J=1,N
        DO 5 I=1,N
           U(I,J)=0.0D0
    5 CONTINUE
      NT = 4
      L = 1
      IBLK=0
   10 DL = 1
      IBLK=IBLK+1
      IF(L.EQ.N) GO TO 20
      S=DABS(A(L,L))+DABS(A(L+1,L+1))
      IF(DABS(A(L+1,L)).GT.EPS*S) DL=2
      IDIM(IBLK)=DL
   20 LL = L+DL-1
      K = L
   30 KM1 = K-1
      DK = 1
      IF(K.EQ.N) GO TO 35
      S=DABS(A(K,K))+DABS(A(K+1,K+1))
      IF(DABS(A(K+1,K)).GT.EPS*S) DK=2
   35 KK= K+DK-1
      IF(K.EQ.1) GO TO 45
      IF(L.NE.1) GO TO 37
      DO 36 J=L,LL
         DO 36 I=K,KK
            U(I,J)= 0.0D0
   36 CONTINUE
   37 CONTINUE
      DO 38 J=L,LL
         DO 38 I=K,KK
            DO 38 IA=1,KM1
               U(I,J)=U(I,J)+A(IA,I)*C(IA,J)
   38 CONTINUE
      DO 39 J=L,LL
         DO 39 I=K,KK
            DO 39 IB=L,LL
               C(I,J)=C(I,J)-U(I,IB)*A(IB,J)
   39 CONTINUE
      IF(IBLK.EQ.1) GO TO 45
      IBM1=IBLK-1
      JJ=1
      DO 44 IC=1,IBM1
         DM=IDIM(IC)
         IDM=DM+JJ-1
         IDK=DK+K-1
         IF(IC.NE.IBM1) GO TO 42
         DO 41 J=JJ,IDM
            DO 41 I=K,KK
               DO 41 IA=K,IDK
                  U(I,J)=U(I,J)+A(IA,I)*C(IA,J)
  41    CONTINUE
  42    CONTINUE
        DO 43 J=L,LL
           DO 43 I=K,KK
              DO 43 IA=JJ,IDM
                 C(I,J)=C(I,J)-U(I,IA)*A(IA,J)
   43   CONTINUE
        JJ=JJ+DM
   44 CONTINUE
   45 IF(DL.EQ.2) GO TO 60
      IF(DK.EQ.2) GO TO 50
      T(1,1) = A(K,K)*A(L,L)-1.0D0
      IF (T(1,1).EQ.0.0D0) RETURN
      C(K,L) = C(K,L)/T(1,1)
      GO TO 90
   50 T(1,1) = A(K,K)*A(L,L)-1.0D0
      T(1,2) = A(KK,K)*A(L,L)
      T(2,1) = A(K,KK)*A(L,L)
      T(2,2) = A(KK,KK)*A(L,L)-1.0D0
      P(1) = C(K,L)
      P(2) = C(KK,L)
      NSYS = 2
      CALL LINEQ (NT,NSYS,T,P,COND,IPVT,WORK)
      C(K,L) = P(1)
      C(KK,L) = P(2)
      GO TO 90
   60 IF(DK.EQ.2) GO TO 70
      T(1,1) = A(K,K)*A(L,L)-1.0D0
      T(1,2) = A(LL,L)*A(K,K)
      T(2,1) = A(L,LL)*A(K,K)
      T(2,2) = A(K,K)*A(LL,LL)-1.0D0
      P(1) = C(K,L)
      P(2) = C(K,LL)
      NSYS= 2
      CALL LINEQ (NT,NSYS,T,P,COND,IPVT,WORK)
      C(K,L) = P(1)
      C(K,LL) = P(2)
      GO TO 90
   70 IF(K.NE.L) GO TO 80
      T(1,1) = A(L,L)*A(L,L)-1.0D0
      T(1,2) = 2.0D0*A(LL,L)*A(L,L)
      T(1,3) = A(LL,L)*A(LL,L)
      T(2,1) = A(L,L)*A(L,LL)
      T(2,2) = A(L,L)*A(LL,LL)+A(L,LL)*A(LL,L)-1.0D0
      T(2,3) = A(LL,LL)*A(LL,L)
      T(3,1) = A(L,LL)*A(L,LL)
      T(3,2) = 2.0D0*A(L,LL)*A(LL,LL)
      T(3,3) = A(LL,LL)*A(LL,LL)-1.0D0
      P(1) = C(L,L)
      P(2) =  C(LL,L)
      P(3) = C(LL,LL)
      NSYS = 3
      CALL LINEQ (NT,NSYS,T,P,COND,IPVT,WORK)
      C(LL,L) = P(2)
      C(L,L) = P(1)
      C(L,LL) = P(2)
      C(LL,LL) = P(3)
      GO TO 90
   80 T(1,1) = A(K,K)*A(L,L)-1.0D0
      T(1,2) = A(KK,K)*A(L,L)
      T(1,3) = A(LL,L)*A(K,K)
      T(1,4) = A(KK,K)*A(LL,L)
      T(2,1) = A(L,L)*A(K,KK)
      T(2,2) = A(KK,KK)*A(L,L)-1.0D0
      T(2,3) = A(K,KK)*A(LL,L)
      T(2,4) = A(KK,KK)*A(LL,L)
      T(3,1) = A(K,K)*A(L,LL)
      T(3,2) = A(KK,K)*A(L,LL)
      T(3,3) = A(K,K)*A(LL,LL)-1.0D0
      T(3,4) = A(LL,LL)*A(KK,K)
      T(4,1) = A(K,KK)*A(L,LL)
      T(4,2) = A(L,LL)*A(KK,KK)
      T(4,3) = A(K,KK)*A(LL,LL)
      T(4,4) = A(KK,KK)*A(LL,LL)-1.0D0
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
   90 K= K + DK
      IF(K.LE.N) GO TO 30
      LDL = L + DL
      IF(LDL.GT.N) RETURN
      DO 100 J=LDL,N
         DO 100 I=L,LL
            C(I,J) = C(J,I)
  100 CONTINUE
      L = LDL
      GO TO 10
C
C     LAST LINE OF DSTSLV
C
      END
