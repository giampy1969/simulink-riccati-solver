      SUBROUTINE CMPRS(NR,NRD,NRT,N,NN,NNPM,M,E,A,B,CQC,R,S,G,F,U,
     X                 WK,WK1,WK2,WK3,EFLAG,SFLAG,INFO)
C
C     *****PARAMETERS:
      INTEGER NR,NRD,NRT,N,NN,NNPM,M,INFO
      CHARACTER EFLAG,SFLAG
      DOUBLE PRECISION E(NR,N),A(NR,N),B(NR,M),CQC(NR,N),
     X    R(NR,M),S(NR,M),G(NRD,NN),F(NRD,NN),U(NRT,NNPM),
     X    WK(NRT,M),WK1(M),WK2(M),WK3(NNPM)
C
C     *****LOCAL VARIABLES:
      INTEGER I,J,JOB,K,NNPI,NPI,NPJ,NPK
C
C     *****FORTRAN FUNCTIONS:
C     NONE
C
C     *****SUBROUTINES CALLED:
C     DSVDC
C
C     --------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE EMPLOYS THE SINGULAR VALUE DECOMPOSITION TO
C     DETERMINE AN ORTHOGONAL MATRIX U, (2N+M) BY (2N+M), SUCH THAT
C
C
C               (     (     )  ( B )     ( 0 )
C               ( U11 ( U12 )  (   )     (   )
C               (     (     )  (-S )  =  ( 0 )
C               (-----(-----)  (---)     (---)
C               ( U21 ( U22 )  ( R )     ( RB)
C
C     AND THEN FORMS THE MATRIX PENCIL
C
C                   ( E  0 )     (     (  A    0 )               )
C        LAMBDA*U11*(      )  -  ( U11*(         ) + U12*(ST BT) )
C                   ( 0 ET )     (     (-CQC -AT )               )
C
C          =:  LAMBDA * F - G
C
C     WHERE T DENOTES MATRIX TRANSPOSE.
C     THIS PENCIL CAN THEN BE USED FOR SOLVING THE CONTINUOUS-TIME
C     GARE.
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
C               E, A, B, CQC, R AND S AS DECLARED IN THE MAIN
C               CALLING PROGRAM DIMENSION STATEMENT;
C
C       NRD     INTEGER
C               ROW DIMENSION OF THE ARRAYS CONTAINING THE MATRICES
C               G AND F AS DECLARED IN THE MAIN CALLING PROGRAM
C               DIMENSION STATEMENT;
C
C       NRT     INTEGER
C               ROW DIMENSION OF THE ARRAYS CONTAINING THE MATRICES
C               U AND WK AS DECLARED IN THE MAIN CALLING PROGRAM
C               DIMENSION STATEMENT;
C
C       N       INTEGER
C               ORDER OF THE SQUARE MATRICES E, A AND CQC
C               ROW DIMENSION OF THE MATRICES B AND S;
C
C       NN      INTEGER
C               ORDER OF THE SQUARE MATRICES G AND F;
C
C       NNPM    INTEGER
C               = NN + M;
C
C       M       INTEGER
C               ORDER OF THE SQUARE MATRIX R
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
C       WK      REAL(NRT,M)
C               SCRATCH ARRAY OF SIZE AT LEAST (NN+M) BY M;
C
C       WK1,WK2 REAL(M)
C               WORKING VECTORS OF SIZE AT LEAST M;
C
C       WK3     REAL(NNPM)
C               WORKING VECTOR OF SIZE AT LEAST NN+M;
C
C       EFLAG   CHARACTER
C               FLAG SET TO 'Y' IF E IS OTHER THAN THE IDENTITY
C               MATRIX;
C
C       SFLAG   CHARACTER
C               FLAG SET TO 'Y' IF S IS OTHER THAN THE ZERO MATRIX;
C
C     ON OUTPUT:
C
C       G       REAL(NRD,NN)
C               MATRIX OF THE COMPRESSED PENCIL AS DEFINED ABOVE;
C
C       F       REAL(NRD,NN)
C               MATRIX OF THE COMPRESSED PENCIL AS DEFINED ABOVE;
C
C       U       REAL(NRT,NNPM)
C               ORTHOGONAL COMPRESSION MATRIX AS DEFINED ABOVE;
C
C       INFO    INTEGER
C               ERROR FLAG WITH MEANING AS FOLLOWS
C               = 0  NORMAL RETURN
C               = NONZERO  IF SINGULAR VALUE DECOMPOSITION COULD
C                          NOT BE CALCULATED.
C
C     *****ALGORITHM NOTES:
C     NONE
C
C     *****HISTORY:
C     THIS SUBROUTINE WAS WRITTEN BY W.F. ARNOLD, NAVAL WEAPONS
C     CENTER, CODE 35104, CHINA LAKE, CA  93555, AS PART OF THE
C     SOFTWARE PACKAGE RICPACK, SEPTEMBER 1983.
C     MODIFIED BY ALAN J. LAUB (UCSB): 01/06/85
C
C     --------------------------------------------------------------
C
C     SET UP MATRIX TO BE COMPRESSED
C
      JOB = 10
      DO 30 J=1,M
          DO 10 I=1,N
              NPI = N + I
              WK(I,J) = B(I,J)
              WK(NPI,J) = 0.0D0
   10     CONTINUE
          DO 20 I=1,M
              NNPI = NN + I
              WK(NNPI,J) = R(I,J)
   20     CONTINUE
   30 CONTINUE
      IF(SFLAG .NE. 'Y' .AND. SFLAG .NE. 'y') GO TO 55
      DO 50 J=1,M
          DO 40 I=1,N
              NPI = N + I
              WK(NPI,J) = -S(I,J)
   40     CONTINUE
   50 CONTINUE
   55 CONTINUE
C
C     CALCULATE COMPRESSION MATRIX U
C
      CALL DSVDC(WK,NRT,NNPM,M,WK1,WK2,U,NRT,G,NRD,WK3,JOB,INFO)
      IF(INFO .NE. 0) RETURN
C
C     FORM MATRIX PENCIL
C
      DO 80 J=1,NN
          MPJ = M + J
          DO 60 I=1,NNPM
              U(I,J) = U(I,MPJ)
   60     CONTINUE
          DO 70 I=1,NN
              G(I,J) = 0.0D0
              F(I,J) = 0.0D0
   70     CONTINUE
   80 CONTINUE
      IF(EFLAG .EQ. 'Y' .OR. EFLAG .EQ. 'y') GO TO 110
      DO 100 I=1,NN
          DO 90 J=1,NN
              F(I,J) = U(J,I)
   90     CONTINUE
  100 CONTINUE
      GO TO 145
  110 CONTINUE
      DO 140 J=1,N
          NPJ = N + J
          DO 130 K=1,N
              NPK = N + K
              DO 120 I=1,N
                  NPI = N + I
                  F(I,J) = F(I,J) + U(K,I)*E(K,J)
                  F(NPI,J) = F(NPI,J) + U(K,NPI)*E(K,J)
                  F(I,NPJ) = F(I,NPJ) + U(NPK,I)*E(J,K)
                  F(NPI,NPJ) = F(NPI,NPJ) + U(NPK,NPI)*E(J,K)
  120         CONTINUE
  130     CONTINUE
  140 CONTINUE
  145 CONTINUE
      DO 180 J=1,N
          NPJ = N + J
          DO 170 K=1,N
              NPK = N + K
              DO 150 I=1,N
                  NPI = N + I
                  G(I,J) = G(I,J) + U(K,I)*A(K,J)
                  G(I,J) = G(I,J) - U(NPK,I)*CQC(K,J)
                  G(NPI,J) = G(NPI,J) + U(K,NPI)*A(K,J)
                  G(NPI,J) = G(NPI,J) - U(NPK,NPI)*CQC(K,J)
                  G(I,NPJ) = G(I,NPJ) - U(NPK,I)*A(J,K)
                  G(NPI,NPJ) = G(NPI,NPJ) - U(NPK,NPI)*A(J,K)
  150         CONTINUE
              DO 160 I=1,M
                  NNPI = NN + I
                  G(K,NPJ) = G(K,NPJ) + U(NNPI,K)*B(J,I)
                  G(NPK,NPJ) = G(NPK,NPJ) + U(NNPI,NPK)*B(J,I)
  160         CONTINUE
  170     CONTINUE
  180 CONTINUE
      IF(SFLAG .NE. 'Y' .AND. SFLAG .NE. 'y') RETURN
      DO 210 J=1,N
          NPJ = N + J
          DO 200 K=1,N
              NPK = N + K
              DO 190 I=1,M
                  NNPI = NN + I
                  G(K,J) = G(K,J) + U(NNPI,K)*S(J,I)
                  G(NPK,J) = G(NPK,J) + U(NNPI,NPK)*S(J,I)
  190         CONTINUE
  200     CONTINUE
  210 CONTINUE
      RETURN
C
C     LAST LINE OF CMPRS
C
      END
