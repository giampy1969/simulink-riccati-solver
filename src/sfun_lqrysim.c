 /*************************************************************************
 *                                                                        *
 *  ALGEBRAIC RICCATI SOLUTION IN SIMULINK VIA C + FORTRAN                *
 *  Giampiero Campa                                                       *
 *  October 2002, January 2009                                            *
 *                                                                        *
 *  Jason Hall, Riccardo Bevilacqua                                       *
 *  Spacecraft Robotics Laboratory                                        *
 *  Naval Postgraduate School, Monterey, CA                               *
 *  September 2008                                                        *
 *                                                                        *
 *  The s-function sfun_lqrysim solves Riccati equations in Simulink      *
 *  without calling the Matlab interpreter.                               *
 *  It is based on the "ricpack" FORTRAN code written by Arnold and Laub  *
 *  in the early eighties and publicly available under the cascade        *
 *  sublibrary of netlib:  http://www.netlib.org/ieeecss/cascade/         *
 *                                                                        *
 *  In particular, the s-function is a level-2 gateway that performs a    *
 *  coordinate transformation to rearrange a general LQR problem (see     *
 *  the Matlab function lqry) into a simplified form, calls a simplified  *
 *  version of the creg.f routine to solve the Riccati equation, and uses *
 *  the solution of the Riccati equation to find the feedback gain matrix *
 *  that solves the original LQR problem.                                 *
 *                                                                        *
 *  The parameters to be given to the s-function are the number of states,*
 *  the number of inputs, and the number of outputs, in that order.       *
 *                                                                        *
 *  The inputs to the s-function are, in order, the matrices A,B,C,D,Q,R  *
 *                                                                        *
 *  The outputs are, in the following order, the solution of the Riccati  *
 *  equation, the feedback gain matrix, and a 2D vector.  The first       *
 *  element in the 2D vector indicates an error when greater than zero or *
 *  a somewhat "unreliable" result when negative. The second element is   *
 *  simply the condition number of the R matrix.                          *
 *                                                                        *
 *  A self-explanatory example (the Simulink scheme sfcndemo_lqry.mdl) is *
 *  presented along with appropriate mex function to allow it to be       *
 *  executed in both a linux and win32 environment.  Additionally, the    *
 *  fortran filesnare included to allow re-mex of this .c file and use in *
 *  generating code for use in a real-time environments such as RTAI.     *
 *                                                                        *
 *  The file creg.f has been simplified by commenting-out several sections*
 *  (see the source code to know why), and that a few other routines have *
 *  been modified to avoid complaints on standard output.                 *
 *  The original code is available at the address given above under the   *
 *  link "creg.f plus dependencies"                                       *
 *                                                                        *
 *  If creating the mex file in the win32 environment, ensure that MinGW  *
 *  is installed with at least the g77 compiler selected (the latest      *
 *  version can be downloaded and installed from www.sourceforge.net).    *
 *  Copy the g77.exe file that will be created using the default install  *
 *  from C:\MinGW\bin to the working directory in MATLAB and verify that  *
 *  the Fortran files and the sfun_lqrysim.c files are also in the working*
 *  directory and then call the following in the MATLAB Command Window:   *
 *                                                                        *
 *  >> eval('! g77 -c *.f')                                               *
 *  >> mex sfun_lqrysim.c *.o                                             *
 *                                                                        *
 *  If creating the mex file in a Linux enviroment, ensure g77 is         *
 *  installed following these instructions from a root (ubuntu) terminal: *
 *                                                                        *
 *  #apt-get install g77-3.4                                              *
 *  #ln -s /usr/bin/g77-3.4 /usr/bin/g77                                  *
 *                                                                        *
 *  Again, verify that the Fortran files and the sfun_lqrysim.c files are *
 *  in the working directory (you do not need the g77.exe file as in the  *
 *  win32 case) and then call the following in the Command Window:        *
 *                                                                        *
 *  >> eval('! g77 -c *.f')                                               *
 *  >> mex sfun_lqrysim.c *.o                                             *
 *                                                                        *
 *  To compile using Real Time Workshop, choose the appropriate .tlc file *
 *  under Configuration Parameters>Real-Time Workshop, then go to         *
 *  Configurations Parameters>Real-Time Workshop>Custom Code and put in   *
 *  the full path to where the source code exists (*.f, *.o files, and    *
 *  sfun_lqrysim.c) in the Include list of the additional:>Include        *
 *  directories: of the RTW GUI window                                    *
 *                                                                        *
 *  For example:  For NPS SRL - ../toolbox/nps_srl/src/lqry/              *
 *                                                                        *
 *  Next, type all the *.o files into Include list of additional:>Source  *
 *  files: GUI window.  Go back to Configurations Parameters> Real-Time   *
 *  Workshop> and select Generate Code.                                   *
 *                                                                        *
 *  If the target is such that the executable can be built on the machine *
 *  where you run the simulation, then at this point you can deselect     *
 *  the Generate Code option under Configurations Parameters >            *
 *  Real-Time Workshop and build the code.                                *
 *                                                                        *
 *  Assumed that the target is RTAI, copy the all the Fortran files,      *
 *  their object files and the sfun_lqrysim.c to the created              *
 *  (mdl name)_rtai folder. If you are not on a windows machine then you  *
 *  neet to move the whole folder to the linux machine (essentially       *
 *  following point 17 of the RTAI_TARGET_HOWTO.TXT file, available on    *
 *  the matlabcentral web site, assuming that the previous things have    *
 *  already been done), then on linux, open a root terminal, browse to    *
 *  the (mdl name)_rtai folder, and then type:                            *
 *                                                                        *
 *  # make -f *.mk                                                        *
 *                                                                        *
 *  The executable will be generated and placed in the folder one up from *
 *  the (mdl name)_rtai folder.                                           *
 *                                                                        *
 *************************************************************************/
 
 
#define S_FUNCTION_LEVEL    2
#undef  S_FUNCTION_NAME
#define S_FUNCTION_NAME     sfun_lqrysim

#include "simstruc.h"

extern void CREG(int *,int *,int *,int *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,double *,int *,int *,int *,int *,int *);
extern void SYMPRD(int *,int *,int *,int *,int *,int *,double *,double *,double *,double *);
extern void MADD(int *,int *,int *,int *,int *,double *,double *,double *);
extern void MSUB(int *,int *,int *,int *,int *,double *,double *,double *);
extern void MMUL(int *,int *,int *,int *,int *,int *,double *,double *,double *);
extern void XTY(int *,int *,int *,int *,int *,int *,double *,double *,double *);
extern void SAVE(int *,int *,int *,int *,double *,double *);
extern void MLINEQ(int *,int *,int *,int *,double *,double *,double *,int *,double *);

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

#define MAX(a,b)  ((a) > (b)  ? (a) : (b))

/* Input Arguments */
#define NUMBER_OF_ARGS        	(3)
#define NUMSTATE                ssGetSFcnParam(S,0)
#define NUMCONTR                ssGetSFcnParam(S,1)
#define NUMOUT                  ssGetSFcnParam(S,2)

#define NO_OUTPUTS              (3) /* number of outputs */
#define NO_INPUTS               (6) /* number of inputs */

static char_T msg[256];

static void mdlInitializeSizes(SimStruct *S)
{
    int     n  = mxGetPr(NUMSTATE)[0];
    int     m  = mxGetPr(NUMCONTR)[0];
    int     p  = mxGetPr(NUMOUT)[0];
    int     dims[2];
    int     nr; 
   
    DECL_AND_INIT_DIMSINFO(di); /* Initializes structure */
    
    ssSetNumSFcnParams(S, NUMBER_OF_ARGS);  /* Number of expected params */
    if (ssGetNumSFcnParams(S) != ssGetSFcnParamsCount(S)) {
        sprintf(msg,"Wrong number of input arguments passed.\n"
        "%d arguments are expected\n", NUMBER_OF_ARGS);
        ssSetErrorStatus(S,msg);
        return;
        }
    
    /* Allow signal dimensions greater than 2 */
    ssAllowSignalsWithMoreThan2D(S);
    
    /* Set number of input and output ports */
    ssSetNumInputPorts(S, NO_INPUTS);
    ssSetNumOutputPorts(S, NO_OUTPUTS);
    
    /*
     * Non-scalar parameter: output dimensions are the same as
     * the parameter dimensions. To support n-D signals, must
     * use a dimsInfo structure to specify dimensions.
    */
    
    /* Assign A Matrix Input */
    dims[0] = n;
    dims[1] = n;
    di.dims = (int *)&dims;
    di.width = dims[0]*dims[1];
    di.numDims = 2;
    if(!ssSetInputPortDimensionInfo(S,  0, &di)) return;
    ssSetInputPortDirectFeedThrough(S, 0, 1);
    
    /* Assign B Matrix Input */
    dims[0] = n;
    dims[1] = m;
    di.dims = (int *)&dims;
    di.width = dims[0]*dims[1];
    di.numDims = 2;
    if(!ssSetInputPortDimensionInfo(S,  1, &di)) return;
    ssSetInputPortDirectFeedThrough(S, 1, 1);
    
    /* Assign C Matrix Input */
    dims[0] = p;
    dims[1] = n;
    di.dims = (int *)&dims;
    di.width = dims[0]*dims[1];
    di.numDims = 2;
    if(!ssSetInputPortDimensionInfo(S,  2, &di)) return;
    ssSetInputPortDirectFeedThrough(S, 2, 1);
    
    /* Assign D Matrix Input */
    dims[0] = p;
    dims[1] = m;
    di.dims = (int *)&dims;
    di.width = dims[0]*dims[1];
    di.numDims = 2;
    if(!ssSetInputPortDimensionInfo(S,  3, &di)) return;
    ssSetInputPortDirectFeedThrough(S, 3, 1);
    
    /* Assign Q Matrix Input */
    dims[0] = p;
    dims[1] = p;
    di.dims = (int *)&dims;
    di.width = dims[0]*dims[1];
    di.numDims = 2;
    if(!ssSetInputPortDimensionInfo(S,  4, &di)) return;
    ssSetInputPortDirectFeedThrough(S, 4, 1);
    
    /* Assign R Matrix Input */
    dims[0] = m;
    dims[1] = m;
    di.dims = (int *)&dims;
    di.width = dims[0]*dims[1];
    di.numDims = 2;
    if(!ssSetInputPortDimensionInfo(S,  5, &di)) return;
    ssSetInputPortDirectFeedThrough(S, 5, 1);
        
    /* Assign S Matrix Output */
    dims[0] = n;
    dims[1] = n;
    di.dims = (int *)&dims;
    di.width = dims[0]*dims[1];
    di.numDims = 2;
    if(!ssSetOutputPortDimensionInfo(S, 0, &di)) return;
    
    /* Assign K Matrix Output */
    dims[0] = m;
    dims[1] = n;
    di.dims = (int *)&dims;
    di.width = dims[0]*dims[1];
    di.numDims = 2;
    if(!ssSetOutputPortDimensionInfo(S, 1, &di)) return;
    
    /* Assign E Matrix Output */
    dims[0] = 2;
    dims[1] = 1;
    di.dims = (int *)&dims;
    di.width = dims[0]*dims[1];
    di.numDims = 2;
    if(!ssSetOutputPortDimensionInfo(S, 2, &di)) return;
    
    ssSetNumSampleTimes(S, 0);
    ssSetOptions(S, SS_OPTION_WORKS_WITH_CODE_REUSE |SS_OPTION_EXCEPTION_FREE_CODE);
        
	nr=(int) (2*n+MAX(m,p));
    ssSetNumRWork(S,25*nr*nr+6*nr); /* number real work vector elements  */

    ssSetNumIWork(S,nr);         /* number integer work vector elements  */

    ssSetNumPWork(S,0);          /* number ptr work vector elements      */
    ssSetNumModes(S,0);          /* number mode work vector elements     */
    ssSetNumNonsampledZCs(S,0);  /* number of nonsampled zero crossing   */
}
 
/* Function to initialize sample times */
static void mdlInitializeSampleTimes(SimStruct *S)
{
    
    ssSetSampleTime(S, 0, INHERITED_SAMPLE_TIME);
    ssSetOffsetTime(S, 0, 0.0);
    ssSetModelReferenceSampleTimeDefaultInheritance(S);
}

/* mdlStart - initialize work vectors */
#define MDL_START
static void mdlStart(SimStruct *S)
{
    int         n  = mxGetPr(NUMSTATE)[0];
    int         m  = mxGetPr(NUMCONTR)[0];
    int         p  = mxGetPr(NUMOUT)[0];
    double     *rv = ssGetRWork(S);
    int        *iv = ssGetIWork(S), nr, i;
    
    nr=(int) (2*n+MAX(m,p));
    
    for (i=0;i<25*nr*nr+6*nr;i++) rv[i]=0;
    for (i=0;i<nr;i++) iv[i]=0;
    
}

/* mdlOutputs - compute the outputs */
static void mdlOutputs(SimStruct *S, int_T tid)
{
    int                i, j, NR, NS, MAXIT, IERR, IBAL, RSTRUC;
    int               *iv  = ssGetIWork(S);
    
    int                 n  = mxGetPr(NUMSTATE)[0];
    int                 m  = mxGetPr(NUMCONTR)[0];
    int                 p  = mxGetPr(NUMOUT)[0];
    double             *rv  = ssGetRWork(S), RTOL, RSD, COND;
    
    double             *P  = ssGetOutputPortRealSignal(S,0);
    double             *K  = ssGetOutputPortRealSignal(S,1);
    double             *E  = ssGetOutputPortRealSignal(S,2);
    
    InputRealPtrsType   a    = ssGetInputPortRealSignalPtrs(S,0);
    InputRealPtrsType   b    = ssGetInputPortRealSignalPtrs(S,1);
    InputRealPtrsType   c    = ssGetInputPortRealSignalPtrs(S,2);
    InputRealPtrsType   d    = ssGetInputPortRealSignalPtrs(S,3);
    InputRealPtrsType   q    = ssGetInputPortRealSignalPtrs(S,4);
    InputRealPtrsType   r    = ssGetInputPortRealSignalPtrs(S,5);

    
     /* Note1: Matrix signals are stored in column major order.
      *  Note2: Access each matrix element by one index not two indices.
     *         For example, if the output signal is a [2x2] matrix signal,
     *        -          -
     *       | y[0]  y[2] |
     *       | y[1]  y[3] |
     *       -           -
     *       Output elements are stored as follows:
     *           y[0] --> row = 0, col = 0
     *           y[1] --> row = 1, col = 0
     *           y[2] --> row = 0, col = 1
     *           y[3] --> row = 1, col = 1 */
    
NR=2*n+MAX(m,p);NS=NR*NR;

/* fill A (#0)*/
for (j=0;j<n;j++) {
	for (i=0;i<n;i++) {
		rv[0*NS +i+j*NR]=*a[i+j*n];
	}
}

/* fill B (#1) */
for (j=0;j<m;j++) {
	for (i=0;i<n;i++) {
		rv[1*NS +i+j*NR]=*b[i+j*n];
	}
}

/* fill C (#2) */
for (j=0;j<n;j++) {
	for (i=0;i<p;i++) {
		rv[2*NS +i+j*NR]=*c[i+j*p];
	}
}

/* fill D (#16.6) */
for (j=0;j<m;j++) {
	for (i=0;i<p;i++) {
		rv[16*NS+6*NR +i+j*NR]=*d[i+j*p];
	}
}

/* fill Q (#4) */
for (j=0;j<p;j++) {
	for (i=0;i<p;i++) {
		rv[4*NS +i+j*NR]=*q[i+j*p];
	}
}

/* fill R (#3) */
for (j=0;j<m;j++) {
	for (i=0;i<m;i++) {
		rv[3*NS +i+j*NR]=*r[i+j*m];
	}
}

/* call symprd to form C'*Q*C (#9) using G (#5)
 SYMPRD( NR, NR, NR, NR, p, n,       C ,       Q ,     CQC ,       G )*/
   symprd_(&NR,&NR,&NR,&NR,&p,&n,&rv[2*NS],&rv[4*NS],&rv[9*NS],&rv[5*NS]);

/* call mmul to form Z (#7) = Q*D
 MMUL( NR, NR, NR, m, p, p,       Q ,             D ,       Z )*/
   mmul_(&NR,&NR,&NR,&m,&p,&p,&rv[4*NS],&rv[16*NS+6*NR],&rv[7*NS]);

/* call xty to form ST (#19.6) = D'*Q'*C
 XTY( NR, NR, NR, p, m, n,       Z ,       C ,            ST )*/
   xty_(&NR,&NR,&NR,&p,&m,&n,&rv[7*NS],&rv[2*NS],&rv[19*NS+6*NR]);

/* call xty to form DQD (#17.6) = D'*Q*D
 XTY( NR, NR, NR, p, m, m,             D ,       Z ,           DQD )*/
   xty_(&NR,&NR,&NR,&p,&m,&m,&rv[16*NS+6*NR],&rv[7*NS],&rv[17*NS+6*NR]);

/* call madd to form RR (#18.6) = R+D'*Q*D
 MADD( NR, NR, NR, m, m,       R ,           DQD ,            RR )*/
   madd_(&NR,&NR,&NR,&m,&m,&rv[3*NS],&rv[17*NS+6*NR],&rv[18*NS+6*NR]);


/* fill RX (#20.6) = (RR+RR')/2 */
for (j=0;j<m;j++) {
	for (i=0;i<m;i++) {
		rv[20*NS+6*NR +i+j*NR]=(rv[18*NS+6*NR +i+j*NR]+rv[18*NS+6*NR +j+i*NR])/2;
	}
}

/* form RI (#11) = eye(m) */
for (j=0;j<m;j++) {
	for (i=0;i<m;i++) {
		rv[11*NS +i+j*NR]=0;
	}
	rv[11*NS +j+j*NR]=1;
}

/* save RX (#20.6) into RR (#18.6)
 SAVE( NR, NR, m, m,            RX ,            RR )*/
   save_(&NR,&NR,&m,&m,&rv[20*NS+6*NR],&rv[18*NS+6*NR]);

/* call mlineq to form RI=inv(RX) 
 MLINEQ( NR, NR, m, m,            RX ,       RI , COND,   IND,         CPERM )*/
   mlineq_(&NR,&NR,&m,&m,&rv[20*NS+6*NR],&rv[11*NS],&COND,&iv[0],&rv[15*NS+4*NR]);

/* call mmul to form RI*S'
 MMUL( NR, NR, NR, n, m, m,       RI ,            ST ,       G )*/
   mmul_(&NR,&NR,&NR,&n,&m,&m,&rv[11*NS],&rv[19*NS+6*NR],&rv[5*NS]);

/* call mmul to form B*RI*S'
 MMUL( NR, NR, NR, n, n, m,       B ,       G ,       F )*/
   mmul_(&NR,&NR,&NR,&n,&n,&m,&rv[1*NS],&rv[5*NS],&rv[6*NS]);

/* call msub to form AX=A-B*RI*S' (#24.6)
 MSUB( NR, NR, NR, n, n,       A ,       F ,            AX )*/
   msub_(&NR,&NR,&NR,&n,&n,&rv[0*NS],&rv[6*NS],&rv[24*NS+6*NR]);

/* call xty to form S*RI*S' (#21.6)
 XTY( NR, NR, NR, m, n, n,            ST ,       G ,           SRS )*/
   xty_(&NR,&NR,&NR,&m,&n,&n,&rv[19*NS+6*NR],&rv[5*NS],&rv[21*NS+6*NR]);

/* form QX (#22.6) and CX (#23.6) = eye(n)*/
for (j=0;j<n;j++) {
	for (i=0;i<n;i++) {
		rv[22*NS+6*NR+i+j*NR]=(rv[9*NS+i+j*NR]+rv[9*NS+j+i*NR])/2-rv[21*NS+6*NR+i+j*NR];
		rv[23*NS+6*NR+i+j*NR]=0;
	}
	rv[23*NS+6*NR+j+j*NR]=1;
}


IERR=0; MAXIT=5; RSTRUC=2; IBAL=0; 
RTOL=0.0; RSD=0.0;

/* call creg to compute riccati solution G (#5) and F (#6) 
 CREG( NR,  n,  m,  n,             AX ,         B ,             CX ,             RR ,             QX ,        G ,        F ,        K ,      ACL ,            CQC ,         E ,        RI ,        RS ,         S ,        UU ,             WK ,           ALFI ,           ALFR ,           BETA ,          CPERM ,         CSCALE ,  RSD,  RTOL,  MAXIT,  RSTRUC,  IBAL,    IND,  IERR)*/
   creg_(&NR, &n, &m, &n, &rv[24*NS+6*NR], &rv[ 1*NS], &rv[23*NS+6*NR], &rv[18*NS+6*NR], &rv[22*NS+6*NR], &rv[5*NS], &rv[6*NS], &rv[7*NS], &rv[8*NS], &rv[22*NS+6*NR], &rv[10*NS], &rv[11*NS], &rv[12*NS], &rv[13*NS], &rv[14*NS], &rv[15*NS+0*NR], &rv[15*NS+1*NR], &rv[15*NS+2*NR], &rv[15*NS+3*NR], &rv[15*NS+4*NR], &rv[15*NS+5*NR], &RSD, &RTOL, &MAXIT, &RSTRUC, &IBAL, &iv[0], &IERR);

/* call xty to form B'* F
 XTY( NR, NR, NR, n, m, n,       B ,       F ,       G )*/
   xty_(&NR,&NR,&NR,&n,&m,&n,&rv[1*NS],&rv[6*NS],&rv[5*NS]);

/* call madd to form ST + B'* F
 MADD( NR, NR, NR, m, n,            ST ,       G ,       UU )*/
   madd_(&NR,&NR,&NR,&m,&n,&rv[19*NS+6*NR],&rv[5*NS],&rv[14*NS]);

/* call mmul to form final state feedback matrix Z = RI*(ST + B'* F)
 MMUL( NR, NR, NR, n, m, m,       RI ,       UU ,       Z )*/
   mmul_(&NR,&NR,&NR,&n,&m,&m,&rv[11*NS],&rv[14*NS],&rv[7*NS]);

/* fill P */
for (j=0;j<n;j++) {
	for (i=0;i<n;i++) {
        P[i+j*n]=rv[6*NS+i+j*NR];
	}
}

/* fill K */
for (j=0;j<n;j++) {
	for (i=0;i<m;i++) {
		K[i+j*m]=rv[7*NS+i+j*NR];
	}
}

E[0]= (double) IERR; 
E[1]=          COND;

}

/* mdlTerminate - called when the simulation is terminated */
static void mdlTerminate(SimStruct *S) {}

/* Trailer information to set everything up for simulink usage */
#ifdef  MATLAB_MEX_FILE  /* Is this file being compiled as a MEX-file?   */
#include "simulink.c"    /* MEX-file interface mechanism                 */
#else
#include "cg_sfun.h"     /* Code generation registration function        */
#endif
