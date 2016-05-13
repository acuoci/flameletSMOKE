flameletSMOKE
============
CFD solver for turbulent non-premixed flames based on the steady laminar flamelets method

Compulsory libraries
--------------------
- Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- RapidXML (http://rapidxml.sourceforge.net/)
- Boost C++ (http://www.boost.org/)
- OpenSMOKE++ (alberto.cuoci@polimi.it)

Optional libraries
------------------
- Intel MKL (https://software.intel.com/en-us/intel-mkl)

Compilation
-----------
Three different options are available to compile the code, according to the level of support for the solution of ODE systems
1. Minimalist: no external, optional libraries are required
2. Minimalist + Intel MKL: only the native OpenSMOKE++ ODE solver can be used, but linear algebra operations are managed by the Intel MKL libraries

1. Instructions to compile the Minimalist version
-------------------------------------------------
1. Open the `mybashrc.minimalist` and adjust the paths to the compulsory external libraries
2. Type: `source mybashrc.minimalist`
3. Compile the steady-state solver: from the `solver/laminarSimpleSMOKE` folder type `wmake`
4. Compile the unsteady solver: from the `solver/laminarPimpleSMOKE` folder type `wmake`

2. Instructions to compile the Minimalist+MKL version
-----------------------------------------------------
1. Open the `mybashrc.minimalist.mkl` and adjust the paths to the compulsory external libraries and the paths to the Intel MKL library
2. Type: `source mybashrc.minimalist.mkl`
3. Compile the steady-state solver: from the `solver/laminarSimpleSMOKE` folder type `wmake`
4. Compile the unsteady solver: from the `solver/laminarPimpleSMOKE` folder type `wmake`

