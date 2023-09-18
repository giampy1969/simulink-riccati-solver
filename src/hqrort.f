      SUBROUTINE HQRORT (NM,N,LOW,IGH,H,WR,WI,Z,IERR)
C
C     *****PARAMETERS:
      INTEGER NM,N,LOW,IGH,IERR
      DOUBLE PRECISION H(NM,N),WR(N),WI(N),Z(NM,N)
C
C     ****LOCAL VARIABLES:
      LOGICAL NOTLAS
      INTEGER I,J,K,L,M,EN,LL,MM,NA,ITS,MP2,ENM2
      DOUBLE PRECISION P,Q,R,S,T,W,X,Y,ZZ,NORM,EPS,EPSP1
C
C     *****FUNCTIONS:
      INTEGER MIN0
      DOUBLE PRECISION DSQRT,DABS,DSIGN
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE IS A MODIFICATION OF THE ALGOL PROCEDURE HQR2,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C     THE FORTRAN CODE IS A MINOR MODIFICATION OF THE EISPACK
C     SUBROUTINE HQR2.
C
C     THIS SUBROUTINE COMPUTES A REAL SCHUR FORM OF A REAL
C     HESSENBERG MATRIX BY THE QR METHOD.  THE ORTHOGONAL
C     TRANSFORMATION WHICH REDUCES A MATRIX TO REAL SCHUR FORM
C     IS PRODUCED BY ACCUMULATING THE SIMILARITY TRANSFORMATIONS
C     COMPUTED BY ORTHES, ORTRAN, AND HQRORT.
C
C     *****PARAMETER DESCRIPTION:
C     ON INPUT:
C        NM               ROW DIMENSION OF THE ARRAYS CONTAINING H
C                         (AND A) AS DECLARED IN THE CALLING PROGRAM
C                         DIMENSION STATEMENT;
C
C        N                ORDER OF THE MATRIX H;
C
C        LOW AND IGH      INTEGERS DETERMINED BY THE BALANCING
C                         SUBROUTINE BALANC.  IF BALANC HAS NOT BEEN
C                         USED SET LOW=1, IGH=N;
C
C        H                N X N ARRAY CONTAINING THE UPPER HESSENBERG
C                         MATRIX;
C
C        Z                N X N ARRAY CONTAINING THE TRANSFORMATION
C                         MATRIX PRODUCED BY ORTRAN AFTER THE
C                         REDUCTION BY ORTHES.
C
C     ON OUTPUT:
C
C        H                HAS BEEN DESTROYED;
C
C        WR,WI            REAL SCRATCH VECTORS OF LENGTH N CONTAINING
C                         THE REAL AND IMAGINARY PARTS, RESPECTIVELY,
C                         OF THE EIGENVALUES.  THE EIGENVALUES
C                         ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE
C                         PAIRS OF VALUES APPEAR CONSECUTIVELY WITH THE
C                         EIGENVALUES HAVING THE POSITIVE IMAGINARY PART
C                         FIRST; IF AN ERROR EXIT IS MADE, THE
C                         EIGENVALUES SHOULD BE CORRECT FOR INDICES
C                         IERR+1,. . .,N;
C
C        Z                N X N ARRAY CONTAINING THE ORTHOGONAL MATRIX
C                         WHICH REDUCES (TO REAL SCHUR FORM) THE MATRIX
C                         ORIGINALLY USED AS INPUT TO ORTHES;
C
C        IERR             AN INTEGER SET TO
C                         ZERO       FOR NORMAL RETURN,
C                         J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                                    DETERMINED AFTER 30 ITERATIONS.
C
C     *****HISTORY:
C     PROGRAM BASED ON EISPACK SUBROUTINE HQR2 AS MODIFIED BY
C     ALAN J. LAUB (ELEC. SYS. LAB., M.I.T., RM. 35-331,
C     CAMBRIDGE, MA 02139, PH.: (617)-253-2125), OCTOBER 1977.
C     MOST RECENT VERSION: OCT. 12, 1977.
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
      IERR = 0
      NORM = 0.0D0
      K = 1
C     ---------- STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ----------
      DO 50 I = 1, N
C
         DO 40 J = K, N
            NORM = NORM + DABS(H(I,J))
   40    CONTINUE
C
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0D0
   50 CONTINUE
C
      EN = IGH
      T = 0.0D0
C     ---------- SEARCH FOR NEXT EIGENVALUES ----------
   60 IF (EN .LT. LOW) GO TO 1001
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     ---------- LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ----------
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = DABS(H(L-1,L-1)) + DABS(H(L,L))
         IF (S .EQ. 0.0D0) S = NORM
         IF (DABS(H(L,L-1)) .LE. EPS * S) GO TO 100
   80 CONTINUE
C     ---------- FORM SHIFT ----------
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITS .EQ. 30) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     ---------- FORM EXCEPTIONAL SHIFT ----------
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = DABS(H(EN,NA)) + DABS(H(NA,ENM2))
      X = 0.75D0 * S
      Y = X
      W = -0.4375D0 * S * S
  130 ITS = ITS + 1
C     ---------- LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ----------
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = DABS(P) + DABS(Q) + DABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         IF (DABS(H(M,M-1)) * (DABS(Q) + DABS(R)) .LE. EPS * DABS(P)
     X    * (DABS(H(M-1,M-1)) + DABS(ZZ) + DABS(H(M+1,M+1)))) GO TO 150
  140 CONTINUE
C
  150 MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0D0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0D0
  160 CONTINUE
C     ---------- DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ----------
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0D0
         IF (NOTLAS) R = H(K+2,K-1)
         X = DABS(P) + DABS(Q) + DABS(R)
         IF (X .EQ. 0.0D0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = DSIGN(DSQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
C     ---------- ROW MODIFICATION ----------
         DO 210 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            IF (.NOT. NOTLAS) GO TO 200
            P = P + R * H(K+2,J)
            H(K+2,J) = H(K+2,J) - P * ZZ
  200       H(K+1,J) = H(K+1,J) - P * Y
            H(K,J) = H(K,J) - P * X
  210    CONTINUE
C
         J = MIN0(EN,K+3)
C     ---------- COLUMN MODIFICATION ----------
         DO 230 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            IF (.NOT. NOTLAS) GO TO 220
            P = P + ZZ * H(I,K+2)
            H(I,K+2) = H(I,K+2) - P * R
  220       H(I,K+1) = H(I,K+1) - P * Q
            H(I,K) = H(I,K) - P
  230    CONTINUE
C     ---------- ACCUMULATE TRANSFORMATIONS ----------
         DO 250 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1)
            IF (.NOT. NOTLAS) GO TO 240
            P = P + ZZ * Z(I,K+2)
            Z(I,K+2) = Z(I,K+2) - P * R
  240       Z(I,K+1) = Z(I,K+1) - P * Q
            Z(I,K) = Z(I,K) - P
  250    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     ---------- ONE ROOT FOUND ----------
  270 H(EN,EN) = X + T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.0D0
      EN = NA
      GO TO 60
C     ---------- TWO ROOTS FOUND ----------
  280 P = (Y - X) / 2.0D0
      Q = P * P + W
      ZZ = DSQRT(DABS(Q))
      H(EN,EN) = X + T
      X = H(EN,EN)
      H(NA,NA) = Y + T
      IF (Q .LT. 0.0D0) GO TO 320
C     ---------- REAL PAIR ----------
      ZZ = P + DSIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0D0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0D0
      WI(EN) = 0.0D0
      X = H(EN,NA)
      S = DABS(X) + DABS(ZZ)
      P = X / S
      Q = ZZ / S
      R = DSQRT(P*P+Q*Q)
      P = P / R
      Q = Q / R
C     ---------- ROW MODIFICATION ----------
      DO 290 J = NA, N
         ZZ = H(NA,J)
         H(NA,J) = Q * ZZ + P * H(EN,J)
         H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
C     ---------- COLUMN MODIFICATION ----------
      DO 300 I = 1, EN
         ZZ = H(I,NA)
         H(I,NA) = Q * ZZ + P * H(I,EN)
         H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
C     ---------- ACCUMULATE TRANSFORMATIONS ----------
      DO 310 I = LOW, IGH
         ZZ = Z(I,NA)
         Z(I,NA) = Q * ZZ + P * Z(I,EN)
         Z(I,EN) = Q * Z(I,EN) - P * ZZ
  310 CONTINUE
C
      GO TO 330
C     ---------- COMPLEX PAIR ----------
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
1000  IERR=EN
1001  RETURN
C
C     LAST LINE OF HQRORT
C
      END
