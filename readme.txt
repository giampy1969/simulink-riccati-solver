 /*************************************************************************
 *                                                                        *
 *  ALGEBRAIC RICCATI SOLUTION IN SIMULINK VIA C + FORTRAN                *
 *  Giampiero Campa                                                       *
 *  October 2002                                                          *
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
 