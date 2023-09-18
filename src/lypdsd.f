      SUBROUTINE LYPDSD(NF,NH,N,F,H,Z,WR,WI,WK,U,IDIM,IER1,IER2)
C
C     *****PARAMETERS:
      INTEGER NF,NH,N,IER1,IER2,IDIM(N)
      DOUBLE PRECISION F(NF,N),U(NF,N),H(NH,N),Z(NF,N),WR(N),WI(N),
     X                 WK(N)
C
C     *****LOCAL VARIABLES:
      INTEGER LOW,IGH,KIN,KOUT
C
C     *****FORTRAN FUNCTIONS:
C     NONE
C
C     *****SUBROUTINES CALLED:
C     ORTHES,ORTRAN,HQRORT,MQFWO(MULWOA),DSTSLV(LINEQ),TRNATA
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE SOLVES THE DISCRETE TIME LYAPUNOV EQUATION
C
C           T
C          F *X*F- X = H.
C
C     BY A MODIFICATION OF THE BARTELS-STEWART ALGORITHM (SEE REFS.(1)
C     AND (2)).
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C        NF,NH            ROW DIMENSIONS OF THE ARRAYS CONTAINING F
C                         (AND Z,U),AND H ,RESPECTIVELY, AS DECLARED IN
C                         THE CALLING PROGRAM DIMENSION STATEMENT;
C
C        N                ORDER OF THE MATRICES F AND C;
C
C        F                N X N (REAL) MATRIX;
C
C        H                N X N SYMMETRIC MATRIX;
C
C        IER1             AN INTEGER VARIABLE; NORMALLY SET IER1 TO 0;
C                         IF IER1 IS SET TO A NON-ZERO INTEGER, THE
C                         REDUCTION OF F TO REAL SCHUR FORM IS SKIPPED
C                         AND THE ARRAYS F AND Z ARE ASSUMED TO CONTAIN
C                         THE REAL SCHUR FORM AND ACCOMPANYING
C                         ORTHOGONAL MATRIX THUS PERMITTING THE
C                         EFFICIENT SOLUTION OF SEVERAL EQUATIONS WITH
C                         DIFFERENT CONSTANT TERMS H;
C
C        IER2             AN INTEGER VARIABLE; NORMALLY SET IER2 TO 0;
C                         IF ONLY A REAL SCHUR FORM OF F AND ASSOCIATED
C                         ORTHOGONAL SIMILARITY MATRIX Z ARE DESIRED,
C                         SET IER2 TO A NON-ZERO INTEGER.
C
C     ON OUTPUT:
C
C        H                N X N ARRAY CONTAINING THE (SYMMETRIC)
C                         SOLUTION X OF THE LYAPUNOV EQUATION;
C
C        F                N X N ARRAY CONTAINING IN ITS UPPER
C                         TRIANGLE AND FIRST SUBDIAGONAL A REAL SCHUR
C                         FORM OF F;
C
C        Z                N X N ARRAY CONTAINING, ON OUTPUT, THE
C                         ORTHOGONAL MATRIX THAT REDUCES F
C                         TO REAL SCHUR FORM;
C
C        WR               REAL SCRATCH VECTOR OF LENGTH N; ON OUTPUT
C                         (WR(I),I=1,N) CONTAINS THE REAL PARTS OF THE
C                         EIGENVALUES OF F AND THUS CAN BE USED TO TEST
C                         THE STABILITY OF F;
C
C        WI               REAL SCRATCH VECTOR OF LENGTH N; ON OUTPUT
C                         CONTAINS THE IMAGINARY PART OF THE
C                         EIGENVALUES OF F;
C
C        WK               REAL SCRATCH VECTOR OF LENGTH N;
C
C        U                N X N REAL SCRATCH ARRAY;
C
C        IDIM             INTEGER SCRATCH VECTOR OF LENGTH N;
C
C        IER1             =0 FOR NORMAL RETURN (IF =0 ON INPUT),
C                         =J IF THE J-TH EIGENVALUE HAS NOT BEEN
C                         DETERMINED IN THE QR ALGORITHM (IF =0 ON
C                         INPUT).
C
C     *****ALGORITHM NOTES:
C     IT IS ASSUMED THAT F HAS NO EIGENVALUES WITH PRODUCT EQUAL TO ZERO
C     (THIS CAN BE CHECKED FROM THE ARRAY WR). THIS IS SUFFICIENT TO
C     GUARANTEE A UNIQUE SOLUTION.
C     IF, MOREOVER, F IS STABLE THEN X IS NONNEGATIVE DEFINITE.
C
C     *****REFERENCES:
C        (1)   BARTELS, R.H., AND G.W. STEWART, SOLUTION
C              OF THE MATRIX EQUATION AX + XB = C,
C              ALGORITHM 432, COMM. ACM, 15(1972),820-826.
C                                                            T
C        (2)   BARRAUD,A.Y., A NUMERICAL ALGORITHM TO SOLVE A XA-X=Q,
C              IEEE TRANSACTIONS ON AUTOMATIC CONTROL, VOL. AC-22, NO.5,
C              OCTOBER 1977,883-885.
C
C     *****HISTORY:
C     WRITTEN BY J.A.K. CARRIG (ELEC. SYS. LAB., M.I.T., RM. 35-427,
C     CAMBRIDGE, MA 02139, PH.: (617) 653-7263, SEPTEMBER 1978.
C     MOST RECENT VERSION: SEPT. 20, 1978.
C
C     modified 8/18/86:bb: deleted references to common block
C     and write statement on error.
C
C     ------------------------------------------------------------------
C
      IF(IER1.NE.0) GO TO 30
      LOW=1
      IGH=N
C
C     COMPUTE A REAL SCHUR FORM OF F
C
      CALL ORTHES(NF,N,LOW,IGH,F,WR)
      CALL ORTRAN(NF,N,LOW,IGH,F,WR,Z)
      CALL HQRORT(NF,N,LOW,IGH,F,WR,WI,Z,IER1)
      IF(IER1.EQ.0) GO TO 20
      RETURN
   20 IF(IER2.NE.0) RETURN
C
C     DO THE BACK SUBSTITUTIONS
C
   30 CALL MQFWO(NH,NF,N,H,Z,WK)
      CALL DSTSLV(NF,NH,N,F,H,U,IDIM)
      CALL TRNATA(NF,N,Z)
      CALL MQFWO(NH,NF,N,H,Z,WK)
      CALL TRNATA(NF,N,Z)
      RETURN
C
C     LAST LINE OF LYPDSD
C
      END
