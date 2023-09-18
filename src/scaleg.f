      SUBROUTINE SCALEG (N,MA,A,MB,B,LOW,IGH,CSCALE,CPERM,WK)
C
C     *****PARAMETERS:
      INTEGER IGH,LOW,MA,MB,N
      DOUBLE PRECISION A(MA,N),B(MB,N),CPERM(N),CSCALE(N),WK(N,6)
C
C     *****LOCAL VARIABLES:
      INTEGER I,IR,IT,J,JC,KOUNT,NR,NRP2
      DOUBLE PRECISION ALPHA,BASL,BETA,CMAX,COEF,COEF2,COEF5,COR,
     *                 EW,EWC,FI,FJ,GAMMA,PGAMMA,SUM,T,TA,TB,TC
C
C     *****FORTRAN FUNCTIONS:
      DOUBLE PRECISION DABS, DLOG10, DSIGN
C     FLOAT
C
C     *****SUBROUTINES CALLED:
C     NONE
C
C     ---------------------------------------------------------------
C
C     *****PURPOSE:
C     SCALES THE MATRICES A AND B IN THE GENERALIZED EIGENVALUE
C     PROBLEM A*X = (LAMBDA)*B*X SUCH THAT THE MAGNITUDES OF THE
C     ELEMENTS OF THE SUBMATRICES OF A AND B (AS SPECIFIED BY LOW
C     AND IGH) ARE CLOSE TO UNITY IN THE LEAST SQUARES SENSE.
C     REF.:  WARD, R. C., BALANCING THE GENERALIZED EIGENVALUE
C     PROBLEM, SIAM J. SCI. STAT. COMPUT., VOL. 2, NO. 2, JUNE 1981,
C     141-152.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C       MA,MB   INTEGER
C               ROW DIMENSIONS OF THE ARRAYS CONTAINING MATRICES
C               A AND B RESPECTIVELY, AS DECLARED IN THE MAIN CALLING
C               PROGRAM DIMENSION STATEMENT;
C
C       N       INTEGER
C               ORDER OF THE MATRICES A AND B;
C
C       A       REAL(MA,N)
C               CONTAINS THE A MATRIX OF THE GENERALIZED EIGENPROBLEM
C               DEFINED ABOVE;
C
C       B       REAL(MB,N)
C               CONTAINS THE B MATRIX OF THE GENERALIZED EIGENPROBLEM
C               DEFINED ABOVE;
C
C       LOW     INTEGER
C               SPECIFIES THE BEGINNING INDEX FOR THE ROWS AND
C               COLUMNS OF A AND B TO BE SCALED;
C
C       IGH     INTEGER
C               SPECIFIES THE ENDING INDEX FOR THE ROWS AND COLUMNS
C               OF A AND B TO BE SCALED;
C
C       CPERM   REAL(N)
C               WORK ARRAY.  ONLY LOCATIONS LOW THROUGH IGH ARE
C               REFERENCED AND ALTERED BY THIS SUBROUTINE;
C
C       WK      REAL(N,6)
C               WORK ARRAY THAT MUST CONTAIN AT LEAST 6*N LOCATIONS.
C               ONLY LOCATIONS LOW THROUGH IGH, N+LOW THROUGH N+IGH,
C               ..., 5*N+LOW THROUGH 5*N+IGH ARE REFERENCED AND
C               ALTERED BY THIS SUBROUTINE.
C
C     ON OUTPUT:
C
C       A,B     CONTAIN THE SCALED A AND B MATRICES;
C
C       CSCALE  REAL(N)
C               CONTAINS IN ITS LOW THROUGH IGH LOCATIONS THE INTEGER
C               EXPONENTS OF 2 USED FOR THE COLUMN SCALING FACTORS.
C               THE OTHER LOCATIONS ARE NOT REFERENCED;
C
C       WK      CONTAINS IN ITS LOW THROUGH IGH LOCATIONS THE INTEGER
C               EXPONENTS OF 2 USED FOR THE ROW SCALING FACTORS.
C
C     *****ALGORITHM NOTES:
C     NONE.
C
C     *****HISTORY:
C     WRITTEN BY R. C. WARD.......
C     Modified 8/86 by Bobby Bodenheimer so that if
C       SUM = 0 (corresponding to the case where the matrix
C       doesn't need to be scaled) the routine returns.
C
C     ---------------------------------------------------------------
C
      IF (LOW .EQ. IGH) GO TO 410
      DO 210 I = LOW,IGH
         WK(I,1) = 0.0D0
         WK(I,2) = 0.0D0
         WK(I,3) = 0.0D0
         WK(I,4) = 0.0D0
         WK(I,5) = 0.0D0
         WK(I,6) = 0.0D0
         CSCALE(I) = 0.0D0
         CPERM(I) = 0.0D0
  210 CONTINUE
C
C     COMPUTE RIGHT SIDE VECTOR IN RESULTING LINEAR EQUATIONS
C
      BASL = DLOG10(2.0D0)
      DO 240 I = LOW,IGH
         DO 240 J = LOW,IGH
            TB = B(I,J)
            TA = A(I,J)
            IF (TA .EQ. 0.0D0) GO TO 220
            TA = DLOG10(DABS(TA)) / BASL
  220       CONTINUE
            IF (TB .EQ. 0.0D0) GO TO 230
            TB = DLOG10(DABS(TB)) / BASL
  230       CONTINUE
            WK(I,5) = WK(I,5) - TA - TB
            WK(J,6) = WK(J,6) - TA - TB
  240 CONTINUE
      NR = IGH-LOW+1
      COEF = 1.0D0/FLOAT(2*NR)
      COEF2 = COEF*COEF
      COEF5 = 0.5D0*COEF2
      NRP2 = NR+2
      BETA = 0.0D0
      IT = 1
C
C     START GENERALIZED CONJUGATE GRADIENT ITERATION
C
  250 CONTINUE
      EW = 0.0D0
      EWC = 0.0D0
      GAMMA = 0.0D0
      DO 260 I = LOW,IGH
         GAMMA = GAMMA + WK(I,5)*WK(I,5) + WK(I,6)*WK(I,6)
         EW = EW + WK(I,5)
         EWC = EWC + WK(I,6)
  260 CONTINUE
      GAMMA = COEF*GAMMA - COEF2*(EW**2 + EWC**2)
     +        - COEF5*(EW - EWC)**2.0D0
      IF (IT .NE. 1) BETA = GAMMA / PGAMMA
      T = COEF5*(EWC - 3.0D0*EW)
      TC = COEF5*(EW - 3.0D0*EWC)
      DO 270 I = LOW,IGH
         WK(I,2) = BETA*WK(I,2) + COEF*WK(I,5) + T
         CPERM(I) = BETA*CPERM(I) + COEF*WK(I,6) + TC
  270 CONTINUE
C
C     APPLY MATRIX TO VECTOR
C
      DO 300 I = LOW,IGH
         KOUNT = 0
         SUM = 0.0D0
         DO 290 J = LOW,IGH
            IF (A(I,J) .EQ. 0.0D0) GO TO 280
            KOUNT = KOUNT+1
            SUM = SUM + CPERM(J)
  280       CONTINUE
            IF (B(I,J) .EQ. 0.0D0) GO TO 290
            KOUNT = KOUNT+1
            SUM = SUM + CPERM(J)
  290    CONTINUE
         WK(I,3) = FLOAT(KOUNT)*WK(I,2) + SUM
  300 CONTINUE
      DO 330 J = LOW,IGH
         KOUNT = 0
         SUM = 0.0D0
         DO 320 I = LOW,IGH
            IF (A(I,J) .EQ. 0.0D0) GO TO 310
            KOUNT = KOUNT+1
            SUM = SUM + WK(I,2)
  310       CONTINUE
            IF (B(I,J) .EQ. 0.0D0) GO TO 320
            KOUNT = KOUNT+1
            SUM = SUM + WK(I,2)
  320    CONTINUE
         WK(J,4) = FLOAT(KOUNT)*CPERM(J) + SUM
  330 CONTINUE
      SUM = 0.0D0
      DO 340 I = LOW,IGH
         SUM = SUM + WK(I,2)*WK(I,3) + CPERM(I)*WK(I,4)
  340 CONTINUE
      IF(SUM.EQ.0.0D0) RETURN
      ALPHA = GAMMA / SUM
C
C     DETERMINE CORRECTION TO CURRENT ITERATE
C
      CMAX = 0.0D0
      DO 350 I = LOW,IGH
         COR = ALPHA * WK(I,2)
         IF (DABS(COR) .GT. CMAX) CMAX = DABS(COR)
         WK(I,1) = WK(I,1) + COR
         COR = ALPHA * CPERM(I)
         IF (DABS(COR) .GT. CMAX) CMAX = DABS(COR)
         CSCALE(I) = CSCALE(I) + COR
  350 CONTINUE
      IF (CMAX .LT. 0.5D0) GO TO 370
      DO 360 I = LOW,IGH
         WK(I,5) = WK(I,5) - ALPHA*WK(I,3)
         WK(I,6) = WK(I,6) - ALPHA*WK(I,4)
  360 CONTINUE
      PGAMMA = GAMMA
      IT = IT+1
      IF (IT .LE. NRP2) GO TO 250
C
C     END GENERALIZED CONJUGATE GRADIENT ITERATION
C
  370 CONTINUE
      DO 380 I = LOW,IGH
         IR = WK(I,1) + DSIGN(0.5D0,WK(I,1))
         WK(I,1) = IR
         JC = CSCALE(I) + DSIGN(0.5D0,CSCALE(I))
         CSCALE(I) = JC
  380 CONTINUE
C
C     SCALE A AND B
C
      DO 400 I = 1,IGH
         IR = WK(I,1)
         FI = 2**IR
         IF (I .LT. LOW) FI = 1.0D0
         DO 400 J =LOW,N
            JC = CSCALE(J)
            FJ = 2**JC
            IF (J .LE. IGH) GO TO 390
            IF (I .LT. LOW) GO TO 400
            FJ = 1.0D0
  390       CONTINUE
            A(I,J) = A(I,J)*FI*FJ
            B(I,J) = B(I,J)*FI*FJ
  400 CONTINUE
  410 CONTINUE
      RETURN
C
C     LAST LINE OF SCALEG
C
      END
