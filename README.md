# simulink-riccati-solver
Algebraic Riccati equation solution in Simulink via C and FORTRAN

[![View simulink-riccati-solver on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/2651-simulink-riccati-solver)

This s-function solves Riccati equations in Simulink&reg; without calling MATLAB&reg; it is fast and supports code generation.

Specifically, the s-function is a level-2 gateway that calls several FORTRAN routines written by Arnold and Laub in the early eighties and publicly available under the cascade sublibrary of netlib. FORTRAN code is included in the Readme file, along with further instructions on how to compile and link the s-function.

Giampiero Campa, October 2002 and January 2009
Riccardo Bevilacqua & Jason Hall, NPS Spacecraft Robotics Lab, October 2008
