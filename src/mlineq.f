      SUBROUTINE MLINEQ (NA,NB,N,M,A,B,COND,IPVT,WORK)
C
C     *****PARAMETERS:
      INTEGER NA,NB,N,M,IPVT(N)
      DOUBLE PRECISION A(NA,N),B(NB,M),COND,WORK(N)
C
C     *****LOCAL VARIABLES:
      INTEGER J
C
C     *****SUBROUTINES CALLED:
C     DGECOM,DGESLM
C
C     ------------------------------------------------------------------
C
C     *****PURPOSE:
C     THIS SUBROUTINE SOLVES THE MATRIX LINEAR EQUATION
C                             A*X = B
C     WHERE A IS AN N X N (INVERTIBLE) MATRIX AND B IS AN N X M
C     MATRIX.  SUBROUTINE DGECOM IS CALLED ONCE FOR THE LU-FACTOR-
C     IZATION OF A AND SUBROUTINE DGESLM IS CALLED M TIMES FOR
C     FORWARD ELIMINATION AND BACK SUBSTITUTION TO PRODUCE THE
C     M COLUMNS OF THE SOLUTION MATRIX X = (A-INVERSE)*B.  AN
C     ESTIMATE OF THE CONDITION OF A IS RETURNED.  SHOULD A BE
C     SINGULAR TO WORKING ACCURACY, COND IS SET TO 1.0D+20.
C
C     *****PARAMETER DESCRIPTION:
C
C     ON INPUT:
C
C        NA,NB            ROW DIMENSIONS OF THE ARRAYS CONTAINING A AND
C                         B, RESPECTIVELY, AS DECLARED IN THE CALLING
C                         PROGRAM DIMENSION STATEMENT;
C
C        N                ORDER OF THE MATRIX A AND NUMBER OF ROWS OF
C                         THE MATRIX B;
C
C        M                NUMBER OF COLUMNS OF THE MATRIX B;
C
C        A                N X N COEFFICIENT MATRIX;
C
C        B                N X M RIGHT HAND SIDE MATRIX.
C
C     ON OUTPUT:
C
C        B                SOLUTION MATRIX  X = (A-INVERSE)*B;
C
C        COND             AN ESTIMATE OF THE CONDITION OF A.
C                         = 1/RCOND  WHERE RCOND IS THE INVERSE
C                         OF THE CONDITION ESTIMATE (SEE THE LINPACK
C                         USER'S GUIDE FOR DETAILS);
C
C        IPVT             PIVOT VECTOR OF LENGTH N (SEE DGECOM
C                         DOCUMENTATION);
C
C        WORK             A REAL SCRATCH VECTOR OF LENGTH N.
C
C     *****APPLICATIONS AND USAGE RESTRICTIONS:
C     (1)THE VALUE OF COND SHOULD ALWAYS BE CHECKED BY THE CALLING
C        PROGRAM.  SHOULD A BE NEAR-SINGULAR (OR SINGULAR TO WORKING
C        ACCURACY) THE DATA SHOULD BE INVESTIGATED FOR POSSIBLE
C        ERRORS.  IF THERE ARE NONE AND THE PROBLEM IS APPARENTLY
C        WELL-POSED AND/OR MEANINGFUL, SINGULAR VALUE ANALYSIS MAY
C        THEN BE A MORE RELIABLE SOLUTION TECHNIQUE (E.G., EISPACK
C        SUBROUTINES  SVD  AND  MINFIT).
C     (2)MLINEQ CAN BE USED TO COMPUTE THE INVERSE OF A:  SIMPLY SOLVE
C        A*X = I WHERE  I  IS THE N X N IDENTITY MATRIX.
C     (3)IF THE SOLUTION TO X*A = B  (X = B*(A-INVERSE)) IS DESIRED,
C        SIMPLY TRANSPOSE THE SOLUTION OF
C         T      T
C        A *X = B .
C
C     *****ALGORITHM NOTES:
C     THE CONTENTS OF A ARE MODIFIED BY THIS SUBROUTINE.  SHOULD THE
C     ORIGINAL COEFFICIENTS OF A BE NEEDED SUBSEQUENTLY, THE
C     CONTENTS OF A SHOULD BE SAVED PRIOR TO THE CALL TO MLINEQ.
C
C     *****HISTORY:
C     WRITTEN BY ALAN J. LAUB (DEP'T. OF ELEC. ENGRG. - SYSTEMS,
C     UNIVERSITY OF SOUTHERN CALIFORNIA, LOS ANGELES, CA 90007,
C     PH.: (213)-743-5535), MAY 1980.
C     MOST RECENT VERSION:  MAY 6, 1980.
C
C     ------------------------------------------------------------------
C
      CALL DGECOM (A,NA,N,IPVT,COND,WORK)
      IF ((1.0D0 + COND) .GT. 1.0D0) GO TO 20
      COND = 1.0D+20
      RETURN
   20 CONTINUE
      COND = 1.0D0/COND
      DO 30 J=1,M
C
C        COMPUTE  (J-TH COLUMN OF X) = (A-INVERSE)*(J-TH COLUMN OF B)
C
         CALL DGESLM (A,NA,N,IPVT,B(1,J))
   30 CONTINUE
      RETURN
      END