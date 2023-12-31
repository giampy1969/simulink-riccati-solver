      SUBROUTINE QZHESW(NM,N,A,B,LOW,IGH,MATZ,Z)
C
      INTEGER I,J,K,L,N,LB,L1,NM,NK1,LOW,IGH,IGHM1,IGHM2
      DOUBLE PRECISION A(NM,N),B(NM,N),Z(NM,N)
      DOUBLE PRECISION R,S,T,U1,U2,V1,V2,RHO
      DOUBLE PRECISION DSQRT,DABS,DSIGN
      LOGICAL MATZ
C
C     THIS SUBROUTINE IS THE FIRST STEP OF THE QZ ALGORITHM
C     FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C     THE CODE IS A MODIFICATION OF SUBROUTINE QZHES FOUND IN EISPACK
C     SO AS TO ACCEPT PARTITIONED MATRICES AS PROVIDED BY THE
C     BALANCING ALGORITHM OF WARD IN SIAM J. SCIEN. STAT. COMPU. 2,
C     141-152 (1981).
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF REAL GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM AND THE OTHER
C     TO UPPER TRIANGULAR FORM USING ORTHOGONAL TRANSFORMATIONS.
C     IT IS USUALLY FOLLOWED BY  QZITW, QZVAL  AND, POSSIBLY,  QZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A CONTAINS A REAL GENERAL MATRIX,
C
C        B CONTAINS A REAL GENERAL MATRIX,
C
C        LOW IS THE STARTING INDEX FOR THE SCALED SUBMATRICES OF
C          A AND B OBTAINED FROM THE BALANCING ALGORITHM.
C          IF THE MATRICES HAVE NOT BEEN BALANCED, SET LOW TO 1.
C
C        IGH IS THE ENDING INDEX FOR THE SCALED SUBMATRICES OF
C          A AND B OBTAINED FROM THE BALANCING ALGORITHM.
C          IF THE MATRICES HAVE NOT BEEN BALANCED, SET IGH TO N.
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO,
C
C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS IF
C          MATZ HAS BEEN SET TO .TRUE.  OTHERWISE, Z IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO R. C. WARD,
C     UNION CARBIDE CORP.-NUCLEAR DIVISION, P.O. BOX Y, BLDG.
C     9704-1, OAK RIDGE, TENNESSEE 37830.
C
C     ------------------------------------------------------------------
C
C     ********** INITIALIZE Z **********
      IF (.NOT. MATZ) GO TO 10
C
      DO 3 I = 1, N
C
         DO 2 J = 1, N
            Z(I,J) = 0.0D0
    2    CONTINUE
C
         Z(I,I) = 1.0D0
    3 CONTINUE
C     ********** REDUCE B TO UPPER TRIANGULAR FORM **********
   10 IF (LOW .EQ. IGH) GO TO 170
      IGHM1 = IGH-1
C
      DO 100 L =LOW, IGHM1
         L1 = L + 1
         S = 0.0D0
C
         DO 20 I = L1, IGH
            S = S + DABS(B(I,L))
   20    CONTINUE
C
         IF (S .EQ. 0.0D0) GO TO 100
         S = S + DABS(B(L,L))
         R = 0.0D0
C
         DO 25 I = L, IGH
            B(I,L) = B(I,L) / S
            R = R + B(I,L)**2
   25    CONTINUE
C
         R = DSIGN(DSQRT(R),B(L,L))
         B(L,L) = B(L,L) + R
         RHO = R * B(L,L)
C
         DO 50 J = L1, N
            T = 0.0D0
C
            DO 30 I = L, IGH
               T = T + B(I,L) * B(I,J)
   30       CONTINUE
C
            T = -T / RHO
C
            DO 40 I = L, IGH
               B(I,J) = B(I,J) + T * B(I,L)
   40       CONTINUE
C
   50    CONTINUE
C
         DO 80 J = LOW, N
            T = 0.0D0
C
            DO 60 I = L, IGH
               T = T + B(I,L) * A(I,J)
   60       CONTINUE
C
            T = -T / RHO
C
            DO 70 I = L, IGH
               A(I,J) = A(I,J) + T * B(I,L)
   70       CONTINUE
C
   80    CONTINUE
C
         B(L,L) = -S * R
C
         DO 90 I = L1, IGH
            B(I,L) = 0.0D0
   90    CONTINUE
C
  100 CONTINUE
C     ********** REDUCE A TO UPPER HESSENBERG FORM, WHILE
C                KEEPING B TRIANGULAR **********
      IF (LOW .EQ. IGHM1) GO TO 170
      IGHM2 = IGH-2
C
      DO 160 K = LOW, IGHM2
         NK1 = IGHM1 - K
C     ********** FOR L=IGH-1 STEP -1 UNTIL K+1 DO -- **********
         DO 150 LB = 1, NK1
            L = IGH - LB
            L1 = L + 1
C     ********** ZERO A(L+1,K) **********
            S = DABS(A(L,K)) + DABS(A(L1,K))
            IF (S .EQ. 0.0D0) GO TO 150
            U1 = A(L,K) / S
            U2 = A(L1,K) / S
            R = DSIGN(DSQRT(U1*U1+U2*U2),U1)
            V1 =  -(U1 + R) / R
            V2 = -U2 / R
            U2 = V2 / V1
C
            DO 110 J = K, N
               T = A(L,J) + U2 * A(L1,J)
               A(L,J) = A(L,J) + T * V1
               A(L1,J) = A(L1,J) + T * V2
  110       CONTINUE
C
            A(L1,K) = 0.0D0
C
            DO 120 J = L, N
               T = B(L,J) + U2 * B(L1,J)
               B(L,J) = B(L,J) + T * V1
               B(L1,J) = B(L1,J) + T * V2
  120       CONTINUE
C     ********** ZERO B(L+1,L) **********
            S = DABS(B(L1,L1)) + DABS(B(L1,L))
            IF (S .EQ. 0.0D0) GO TO 150
            U1 = B(L1,L1) / S
            U2 = B(L1,L) / S
            R = DSIGN(DSQRT(U1*U1+U2*U2),U1)
            V1 =  -(U1 + R) / R
            V2 = -U2 / R
            U2 = V2 / V1
C
            DO 130 I = 1, L1
               T = B(I,L1) + U2 * B(I,L)
               B(I,L1) = B(I,L1) + T * V1
               B(I,L) = B(I,L) + T * V2
  130       CONTINUE
C
            B(L1,L) = 0.0D0
C
            DO 140 I = 1, IGH
               T = A(I,L1) + U2 * A(I,L)
               A(I,L1) = A(I,L1) + T * V1
               A(I,L) = A(I,L) + T * V2
  140       CONTINUE
C
            IF (.NOT. MATZ) GO TO 150
C
            DO 145 I = LOW, IGH
               T = Z(I,L1) + U2 * Z(I,L)
               Z(I,L1) = Z(I,L1) + T * V1
               Z(I,L) = Z(I,L) + T * V2
  145       CONTINUE
C
  150    CONTINUE
C
  160 CONTINUE
C
  170 RETURN
C     ********** LAST CARD OF QZHESW **********
      END
