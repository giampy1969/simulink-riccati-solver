      SUBROUTINE RINV(NR,NRD,N,NN,M,E,A,B,CQC,RI,G,F,WK1,WRK,RDFLG,
     X                RFLAG,EFLAG,CTIME)
C
C     *****PARAMETERS:
      INTEGER NR,NRD,N,NN,M
      CHARACTER RDFLG,RFLAG,EFLAG
      DOUBLE PRECISION E(NR,N),A(NR,N),B(NR,M),CQC(NR,N),RI(NR,M),
     X    G(NRD,NN),F(NRD,NN),WK1(NR,N),WRK(N)
      LOGICAL CTIME
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,K,NPI,NPJ
C
C     *****FORTRAN FUNCTIONS:
C     NONE
C
C     *****SUBROUTINES CALLED:
C     MULB, TRNATB
C
C     --------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE FORMS THE MATRIX PENCIL:
C
C                 ( E   0 )     (   A   -B*RI*BT )
C          LAMBDA*(       )  -  (                )
C                 ( 0  ET )     ( -CQC      -AT  )
C
C      =:  LAMBDA * F  -  G
C
C     WHERE T DENOTES MATRIX TRANSPOSE.
C     THIS SUBROUTINE IS USEFUL IN SOLVING THE CONTINUOUS-TIME GARE.
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
C               E, A, B, CQC, RI AND WK1 AS DECLARED IN THE MAIN
C               CALLING PROGRAM DIMENSION STATEMENT;
C
C       NRD     INTEGER
C               ROW DIMENSION OF THE ARRAYS CONTAINING THE MATRICES
C               F AND G AS DECLARED IN THE MAIN CALLING PROGRAM
C               DIMENSION STATEMENT;
C
C       N       INTEGER
C               ORDER OF THE SQUARE MATRICES E, A AND CQC
C               ROW DIMENSION OF THE MATRIX B;
C
C       NN      INTEGER
C               SIZE OF THE MATRIX PENCIL;
C
C       M       INTEGER
C               ORDER OF THE SQUARE MATRIX RI
C               COLUMN DIMENSION OF THE MATRIX B;
C
C       E       REAL(NR,N)
C               MODEL DESCRIPTOR MATRIX;
C
C       A       REAL(NR,N)
C               = A - B*RI*ST IN THE GENERALIZED CASE;
C
C       B       REAL(NR,M)
C               MODEL INPUT MATRIX;
C
C       CQC     REAL(NR,N)
C               =  CT*Q*C - S*RI*ST IN THE GENERALIZED CASE;
C
C       RI      REAL(NR,M)
C               INVERSE OF THE CONTROL WEIGHTING MATRIX;
C
C       WK1     REAL(NR,N)
C               SCRATCH ARRAY OF SIZE AT LEAST M BY N;
C
C       WRK     REAL(N)
C               WORKING VECTOR OF SIZE AT LEAST N;
C
C       RDFLG   CHARACTER
C               FLAG SET TO 'Y' IF RI IS A DIAGONAL MATRIX;
C
C       RFLAG   CHARACTER
C               FLAG SET TO 'Y' IF RI IS OTHER THAN THE IDENTITY
C               MATRIX;
C
C       EFLAG   CHARACTER
C               FLAG SET TO 'Y' IF E IS OTHER THAN THE IDENTITY
C               MATRIX;
C
C       CTIME    LOGICAL
C               = .TRUE.  FOR CONTINUOUS-TIME SYSTEM
C               = .FALSE. FOR DISCRETE-TIME SYSTEM.
C
C     ON OUTPUT:
C
C       G       REAL(NRD,NN)
C               PENCIL MATRIX AS DEFINED ABOVE;
C
C       F       REAL(NRD,NN)
C               PENCIL MATRIX AS DEFINED ABOVE.
C
C     *****ALGORITHM NOTES:
C     NONE
C
C     *****HISTORY:
C     THIS SUBROUTINE WAS WRITTEN BY W.F. ARNOLD, NAVAL WEAPONS
C     CENTER, CODE 35104, CHINA LAKE, CA  93555, AS PART OF THE
C     SOFTWARE PACKAGE RICPACK, SEPTEMBER 1983.
C
C     --------------------------------------------------------------
C
      DO 20 I=1,NN
          DO 10 J=1,NN
              F(J,I) = 0.0D0
              G(J,I) = 0.0D0
   10     CONTINUE
   20 CONTINUE
C
C     FORM G MATRIX
C
      DO 40 J=1,N
          NPJ = N + J
          DO 30 I=1,N
              NPI = N + I
              G(I,J) = A(I,J)
              G(NPI,J) = -CQC(I,J)
              G(NPI,NPJ) = -A(J,I)
   30     CONTINUE
          F(J,J) = 1.0D0
          F(NPJ,NPJ) = 1.0D0
   40 CONTINUE
      IF(M .EQ. 0) GO TO 120
      IF (RDFLG . NE. 'Y' .AND. RDFLG .NE. 'y') GO TO 70
      DO 60 J=1,N
          DO 50 I=1,M
              WK1(I,J) = RI(I,I)*B(J,I)
   50     CONTINUE
   60 CONTINUE
      GO TO 80
   70 CONTINUE
      CALL TRNATB(NR,NR,N,M,B,WK1)
      IF(RFLAG .NE. 'Y' .AND. RFLAG .NE. 'y') GO TO 80
      CALL MULB(NR,NR,M,M,N,RI,WK1,WRK)
   80 CONTINUE
      DO 110 J=1,N
          NPJ = N + J
          DO 100 K=1,N
              DO 90 I=1,M
                  G(K,NPJ) = G(K,NPJ) - B(K,I)*WK1(I,J)
   90         CONTINUE
  100     CONTINUE
  110 CONTINUE
  120 CONTINUE
      IF(EFLAG .NE. 'Y' .AND. EFLAG .NE. 'y') RETURN
C
C     FORM F MATRIX FOR CASE E NOT IDENTITY
C
      DO 140 J=1,N
          NPJ = N + J
          DO 130 I=1,N
              NPI = N + I
              F(I,J) = E(I,J)
              F(NPI,NPJ) = E(J,I)
  130     CONTINUE
  140 CONTINUE
      RETURN
C
C     LAST LINE OF RINV
C
      END
