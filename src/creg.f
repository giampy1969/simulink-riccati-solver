      SUBROUTINE CREG(NR,N,M,L,A,B,C,R,Q,G,F,Z,ACL,
     1                 CQC,E,RI,RS,S,U,WK,ALFI,ALFR,BETA,
     2                 CPERM,CSCALE,RSD,RTOL,MAXIT,
     3                 RSTRUC,IBAL,IND,IERR)
C
C       FUNCTION:
CF
CF       The subroutine CREG designs a steady state optimal linear
CF       regulator for the continous-time linear system:
CF
CF                   X'(t) = A*X(t) + B*U(t)
CF                    Y(t) = C*X(t)
CF
CF       which minimizes 
CF
CF                      infinity  T              T
CF          J = Integral(       (Y (t)*Q*Y(t) + U (t)*R*U(t))dt
CF                      0
CF
CF       The optimal gain matrix Z is of the form
CF
CF                     -1   T
CF                Z = R  * B  *  F
CF
CF       where F is the solution to the continuous time 
CF       Algebraic Riccati Equation (ARE)
CF
CF         T                -1  T      T
CF        A *F + F*A - F*B*R  *B *F + C *Q*C = 0.
CF
CF      The closed loop gain matrix ACL of the system is then given
CF       by:
CF                ACL = A - B*Z.
CF
C       USAGE:
CU
CU          The subroutine CREG is used to design a regulator for a system.
CU      It can be invoked by the call
CU
CU      SUBROUTINE CREG(NR,N,M,L,A,B,C,R,Q,G,F,Z,ACL,
CU     1                CQC,E,RI,RS,S,U,WK,ALFI,ALFR,BETA,
CU     2                CPERM,CSCALE,RSD,RTOL,MAXIT,
CU     3                RSTRUC,IBAL,IND,IERR,ERRMSG)
CU
CU      where the variables are described below. When the subroutine is
CU      called, it is only necessary to have NR, N, M, L, A, B, C, R, Q,
CU      RTOL, MAXIT, RSTRUC, and IBAL initialized.
CU
C       INPUTS:
CI
CI         NR - INTEGER variable which is the row dimension of
CI              all the arrays as declared in the calling program.
CI
CI         N - INTEGER variable; the order the system.
CI
CI         M - INTEGER variable; the number of system inputs.
CI
CI         L - INTEGER variable; the number of system outputs.
CI
CI         A - DP variable; the N by N `A' (state) matrix of the system.
CI
CI         B - DP variable; the N by M `B' (input) matrix of the system.
CI
CI         C - DP variable; the L by N `C' (output) matrix of the system.
CI
CI         R - DP variable; the M by M control weighting matrix.
CI 
CI         Q - DP variable; the L by L output weighting matrix. 
CI
CI         G - DP variable; a scratch array.
CI
CI         CQC - DP variable; a scratch array, it contains the matrix
CI               product CT*Q*C, where T is the transpose.
CI
CI         E - DP variable; a scratch array.
CI 
CI         S - DP variable; a scratch array.
CI
CI         U - DP variable; a scratch array.
CI
CI         WK - DP variable; a scratch array.
CI
CI         RTOL - DP variable; set this to the worst condition R can be 
CI                and still be considered invertible. The condition RTOL
CI                is compared to is 1/RCOND where RCOND is the condition
CI                estimate used in LINPACK. If RTOL=0, then a default 
CI                value of 1.0D+6 is used.
CI
CI         MAXIT - the maximum number of iterations to use for improving
CI                 the accuracy of the solution by Newton's method.
CI
CI         RSTRUC - INTEGER variable indicating the structure of the 
CI                  R matrix:
CI                        RSTRUC = 0 if R is the identity;
CI                        RSTRUC = 1 if R is diagonal but
CI                                   not the identity;
CI                        RSTRUC = 2 otherwise.
CI
CI         IBAL - INTEGER variable determining if the equation is
CI                balanced:
CI                         IBAL = 1 for WARD (pre-QZ) balancing.
CI                         IBAL = 0 for no balancing.
CI
C       OUTPUTS:
CO
CO         F - DP variable; F contains the N by N solution to the ARE.
CO
CO         Z - DP variable; Z contains the M by N optimal gain matrix.
CO
CO         ACL - DP variable; ACL contains the N by N closed loop 
CO               matrix for the system, which is A-B*Z.
CO
CO         RI - DP variable; contains the M by M matrix R-inverse; 
CO              however, the program may not have have used 
CO              R-inverse in the calculation of the ARE solution.
CO
CO         RS - DP variable; contains the N by N residual matrix of
CO              of the ARE with F as the solution.
CO
CO         RSD - DP variable; the 1-Norm of the residual matrix RS 
CO               divided by the 1-Norm of the solution matrix F.
CO
CO         MAXIT - INTEGER variable; on output, MAXIT contains the
CO                 number Newton iterations actually performed by
CO                 NEWT. A number less than the original MAXIT is
CO                 not necessarily bad since the solution may have
CO                 converged.
CO
CO         ALFR - DP variable; this matrix contains the real part of
CO                the eigenvalues of ACL.
CO
CO         ALFI - DP variable; this matrix contains the imaginary part of
CO                the eigenvalues of ACL.
CO
CO         IERR - An INTEGER variable indicating whether an error or
CO                warning has occured during execution.
CO                    IERR = 0 for a normal return.
CO                    IERR < 0 is a warning, indicating abnormal
CO                           but nonfatal conditions. NOTE that if
CO                           a warning occurs, then it is the last
CO                           warning condition which occured, and
CO                           therefore may be misleading, since an
CO                           earlier warning may have caused other
CO                           warnings and then been masked.
CO                    IERR > 0 indicates a fatal error.
CO
CO         ERRMSG = a CHARACTER variable which is only meaningful if
CO                  IERR .NE. 0, it contains a description of the 
CO                  error or warning.
CO
CO
C       ALGORITHM:
CA
CA          The algorithm and much of the code is adapted (copied where
CA          possible) from the interactive driver program to RICPACK.
CA          See:
CA              1. Arnold, W. F., and A. J. Laub, "Generalized
CA                 eigenproblem algorithms and software for
CA                 algebraic Riccati equations," Proceedings of
CA                 the IEEE, vol. 72, no. 12, December 1984,
CA                 pp. 1746-1754.
CA
CA              2. The RICPACK code.
CA
CA            Note that this program currently only solves the ARE for
CA            the case where S=0 and E=I.
CA
C       MACHINE DEPENDENCIES:
CM
CM          None.
CM
C       HISTORY:
CH
CH      written by:             Bobby Bodenheimer
CH      date:                   September 1986
CH      current version:        1.2
CH      modifications:          zero'ed out cross term S - jdb - 5/11/88
CH                              defined undefined variables - jdb - 5/13/88
CH
CH      modified by:            Giampiero Campa
CH      date:                   October 2002
CH      current version:        1.2. (simulink compatible version)
CH      modifications:          heavily commented to avoid:
CH                              1) any possible complaint on standart output 
CH                              2) initializing S to 0 since it is passed as 0
CH                                 and it is not changed by any routine
CH                              3) initializing E to identity since all the 
CH                                 subroutines are called with EFLAG='N'
CH                              4) calculation of CQC since C is identity 
CH                              5) execution of newt and resid (see below)
CH                              6) computation of closed loop matrix
CH
C       ROUTINES CALLED:
CC
CC        RINV, RICSOL, CMPRS, SAVE, MADD, MMUL, TRNATB,
CC        MLINEQ, FBGAIN, NEWT, RESID, MSUB, XTY, ICNVRT
CC
CC
C----------------------------------------------------------------------
C       written for:    The CASCADE Project
C                       Oak Ridge National Laboratory
C                       U.S. Department of Energy
C                       contract number DE-AC05-840R21400
C                       subcontract number 37B-07685C S13
C                       organization:   The University of Tennessee
C----------------------------------------------------------------------
C       THIS SOFTWARE IS IN THE PUBLIC DOMAIN
C       NO RESTRICTIONS ON ITS USE ARE IMPLIED
C----------------------------------------------------------------------
C
C  Global Variables.
C
      INTEGER               NR
      INTEGER               RSTRUC
      INTEGER               IERR
      INTEGER               IND(NR)
      INTEGER               IBAL
      INTEGER               INFO
      INTEGER               L
      INTEGER               M
      INTEGER               MAXIT
      INTEGER               N

      DOUBLE PRECISION      A(NR,NR)
      DOUBLE PRECISION      B(NR,NR)
      DOUBLE PRECISION      C(NR,NR)
      DOUBLE PRECISION      CQC(NR,NR)
      DOUBLE PRECISION      E(NR,NR)
      DOUBLE PRECISION      F(NR,NR)
      DOUBLE PRECISION      G(NR,NR)
      DOUBLE PRECISION      Q(NR,NR)
      DOUBLE PRECISION      R(NR,NR)
      DOUBLE PRECISION      RI(NR,NR)
      DOUBLE PRECISION      RS(NR,NR)
      DOUBLE PRECISION      S(NR,NR)
      DOUBLE PRECISION      U(NR,NR)
      DOUBLE PRECISION      WK(NR,NR)
      DOUBLE PRECISION      Z(NR,NR)
      DOUBLE PRECISION      ACL(NR,NR)
C
      DOUBLE PRECISION      ALFI(NR)
      DOUBLE PRECISION      ALFR(NR)
      DOUBLE PRECISION      BETA(NR)
      DOUBLE PRECISION      CPERM(NR)
      DOUBLE PRECISION      CSCALE(NR)
C
      DOUBLE PRECISION      RSD
      DOUBLE PRECISION      RTOL
C
C   Local Variables.
C
      INTEGER               I
      INTEGER               J
      INTEGER               LENGTH
      INTEGER               MAXNL
      INTEGER               MERR
      INTEGER               NN
      INTEGER               NNPM
      INTEGER               NOUT
C
      PARAMETER(NOUT=6)
C
      DOUBLE PRECISION      COND
      DOUBLE PRECISION      SEP
C
      CHARACTER             RDFLG
      CHARACTER             RFLAG
C      CHARACTER*30          STRING
C
C   Begin
C
      IERR = 0
      IF (RTOL.EQ.0.0D0) RTOL = 1.0D+6
      NN = 2*N
      NNPM = NN + M
C
C      Check system size.
C
C      IF (NNPM.GT.NR) THEN
C         CALL ICNVRT(0,NNPM,STRING,LENGTH,MERR)
C         IERR = 1
C         ERRMSG = 'CREG: Fatal Error: System too large.'//
C     1            ' For this system, NR should be at '//
C     2            ' at least '//STRING(1:LENGTH)
C         RETURN
C      END IF
C
C    Set the structure flags for R.
C
      IF (RSTRUC.EQ.0) THEN
         RFLAG = 'N'
         RDFLG = 'Y'
      ELSE IF (RSTRUC.EQ.1) THEN
         RFLAG = 'Y'
         RDFLG = 'Y'
      ELSE IF (RSTRUC.EQ.2) THEN
         RFLAG = 'Y'
         RDFLG = 'N'
      ELSE
         IERR = 2
C         ERRMSG = ' CREG: Fatal Error: RSTRUC must be '//
C     1            '0, 1, or 2.'
         RETURN
      END IF
C
C  Zero out the cross term S
C
C      DO 4 I = 1, M
C         DO 4 J = 1, N
C            S(J,I) = 0.D0
C 4    CONTINUE
C
C   Make E the identity.
C
C      MAXNL = MAX(N,L)
C      DO 6 J=1,MAXNL
C         DO 5 I=1,MAXNL
C            E(I,J) = 0.0D0
C 5       CONTINUE
C         E(J,J) = 1.0D0
C 6    CONTINUE
C
C                           T
C  Form the matrix product C *Q*C, storing it in CQC. Use G as
C   temporary storage.
C
C      CALL SYMPRD (NR,NR,NR,NR,L,N,C,Q,CQC,G)
C
      IF (RSTRUC.EQ.0) THEN
C
C        R = Identity.
C
         CALL SAVE(NR,NR,M,M,R,RI)
         CALL RINV(NR,NR,N,NN,M,
     1             E,A,B,CQC,RI,G,F,RS,CPERM,
     2             RDFLG,RFLAG,'N',.TRUE.)
C
C        R is diagonal.
C
      ELSE IF (RSTRUC.EQ.1) THEN
            DO 10 I=1,M
               RI(I,I) = 1.0D0/R(I,I)
 10         CONTINUE
            CALL RINV(NR,NR,N,NN,M,
     1                E,A,B,CQC,RI,G,F,RS,CPERM,
     2                RDFLG,RFLAG,'N',.TRUE.)
C
      ELSE IF (RSTRUC.EQ.2) THEN
C
C   No special structure in R to exploit.
C
            DO 30 J=1,M
               DO 20 I=1,M
                  RI(I,J) = 0.0D0
 20            CONTINUE
               RI(J,J) = 1.0D0
 30         CONTINUE
C
C        Check the condition of R with respect to inversion.
C          RI contains R-inverse, but if R is ill-conditioned,
C          we won't use it.
C
            CALL MLINEQ(NR,NR,M,M,R,RI,
     1                  COND,IND,CPERM)
C
C          R ill-conditioned.
C
            IF (COND.GE.RTOL) THEN
                CALL CMPRS(NR,NR,NR,N,NN,
     1                     NNPM,M,E,A,B,CQC,R,
     2                     S,G,F,U,WK,CPERM,CSCALE,
     3                     BETA,'N','N',INFO)
                 IF (INFO.NE.0) THEN
                    IERR = 100 + INFO
C                    ERRMSG = 'CREG: Fatal Error: Attempt to '//
C     1                       'compress R failed. IERR = 100 +'//
C     2                       'INFO, INFO returned from DSVDC.'
                    RETURN
                 END IF
            ELSE
C
C           R not ill-conditioned.
C
               CALL RINV(NR,NR,N,NN,M,
     1                   E,A,B,CQC,RI,G,F,RS,CPERM,
     2                   RDFLG,RFLAG,'N',.TRUE.)
            END IF
      END IF
C
C  Compute the Riccati solution.
C
      CALL RICSOL(NR,NR,NN,N,G,F,E,Z,
     1            ALFR,ALFI,BETA,CPERM,CSCALE,
     2            IND,-1,IBAL,.TRUE.,'N')
C
C     Compute the generalized eigenvalues in case the routine
C     bombs later.
C
C      DO 40, I=1,N
C         ALFR(I) = ALFR(I)/BETA(I)
C         ALFI(I) = ALFI(I)/BETA(I)
C40    CONTINUE
C
C     Check the RICSOL computation for errors.
C
      IF (IND(1).NE.0) THEN
         IERR = 200 + IND(1)
C         ERRMSG = 'CREG: Fatal Error: More than 50 '//
C     1            'iterations required by QZITW. '//
C     2            'IERR = 200 + IND(1), where IND(1) '//
C     3            ' is the ierr from QZITW.'
         RETURN
      END IF
      IF (CPERM(1).EQ.1.0D+20) THEN
C          CALL SAVE(NR,NR,N,N,Z,F)
          IERR = 300
C          ERRMSG = ' CREG: Fatal Error: Schur Vector '//
C     1             'singular to working precision, '//
C     2             'probably because the original system '//
C     3             'is not stabilizable. Hence solution '//
C     4             'by this method is not possible. The '//
C     5             'Schur matrix is returned in F.'
          RETURN
      END IF
      IF (IND(2).NE.0) THEN
         IERR = 400
C         ERRMSG = 'CREG: Fatal Error: Convergence '//
C     1             ' failure in ordering routine '//
C     2             'ORDER. '
         RETURN
      END IF
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  The part below has been almost completely cut off since it does 
C  not really improve upon the solution when R is diagonal, and it 
C  often screws everything up when L=N and R is nondiagonal.
C  In any case it considerably slows down the routine.
C
C  Giampiero Campa, October 2002
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Attempt to improve the solution by iterating
C  with Newton's method.
C
C       CALL SAVE(NR,NR,N,N,F,U)
C       CALL NEWT(NR,NR,NR,NR,N,M,E,A,B,CQC,R,S,RI,RS,U,
C     1           G,F,Z,ALFR,ALFI,BETA,IND,RTOL,'N',RFLAG,
C     2           RDFLG,'N',SEP,.TRUE.,MAXIT,INFO,NOUT)
C       IF (INFO.NE.0) THEN
C          IERR = -200 - INFO
CC          ERRMSG = 'CREG: Warning: Newton iteration '//
CC     1             'failed. Solution after last '//
CC     2             'iteration given. IERR = -200 - INFO '//
CC     3             'where INFO is returned from NEWT.'
C      END IF
C
C   Compute the Residual of the refined ARE solution.
C
C       CALL RESID(NR,NR,NR,NR,N,M,E,A,B,CQC,R,S,RI,RS,
C     1            U,G,Z,CPERM,IND,RTOL,'N',RFLAG,RDFLG,
C     2            'N',RSD,.TRUE.,NOUT)
C       IF (IND(1).NE.0) THEN
C          IERR = -300
CC         ERRMSG = 'CREG: Warning: Residual estimation '//
CC    1             'may be invalid. IPVT(1) from RESID '//
CC    2             'is non-zero.'
C       END IF
C       CALL SAVE(NR,NR,N,N,U,F)
C
C   Compute the optimal gain matrix Z.
C
       CALL FBGAIN(NR,NR,NR,N,M,A,B,E,R,RI,S,F,Z,G,CPERM,
     1             IND,'N',RDFLG,RFLAG,'N',.TRUE.)
C
C  Calculate the Closed Loop Matrix for the system.
C              ACL = A - B*Z
C
C       CALL MMUL(NR,NR,NR,N,N,M,B,Z,G)
C       CALL MSUB(NR,NR,NR,N,N,A,G,ACL)
C
C  Last Lines of CREG.
C
       RETURN
       END
