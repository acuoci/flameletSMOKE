flameletSMOKE
============
CFD solver for turbulent non-premixed flames based on the Steady Laminar Flamelets method

If you use flameletSMOKE for your publications, we kindly ask you to cite the following paper:

> Cuoci, A., Frassoldati, A., Faravelli, T., Ranzi, E., 
> OpenSMOKE++: An object-oriented framework for the numerical modeling of reactive systems with detailed kinetic mechanisms 
> (2015) Computer Physics Communications, 192, pp. 237-264, DOI: 10.1016/j.cpc.2015.02.014


Compulsory libraries
--------------------
- Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- RapidXML (http://rapidxml.sourceforge.net/)
- Boost C++ (http://www.boost.org/)

Optional libraries (strongly recommended)
-----------------------------------------
- Intel MKL (https://software.intel.com/en-us/intel-mkl)

External software
-----------------
In order to generate the lookup tables, the [`OpenSMOKE++Suite`][1] framework is required (alberto.cuoci@polimi.it).

Compilation
-----------
The flameletSMOKE solver is available for OpenFOAM versions 2.2, 2.3, 2.4, 4.x and dev
We strongly recommend to use OpenSMOKE-4.x or OpenSMOKE-dev. Please note that since April 2017
OpenFOAM versions lower than 4.x are not longer supported.

Two different options are available to compile the code:
1. Minimalist + Intel MKL (recommended): the Intel MKL libraries are used to carry out the most CPU expensive operations
2. Minimalist: the most expensive CPU operations are carried out by the internal functions and classes directly available in flameletSMOKE

1. Instructions to compile the Minimalist+MKL version (recommended)
-------------------------------------------------------------------
1. Open the `mybashrc.minimalist.mkl` and adjust the paths to the Intel MKL library
2. Type: `source mybashrc.minimalist.mkl`
3. Compile the flameletSMOKE library: from the `libs/thermophysicalModels/basic` folder type `wmake`
4. Compile the steady-state solver: from the `solvers/flameletSimpleSMOKE` folder type `wmake`
5. Compile the unsteady solvers:
  - from the `solver/flameletPimpleSMOKE` folder type `wmake`
  - from the `solver/flameletPisoSMOKE` folder type `wmake` (available only for OpenFOAM 2.2, 2.3, and 2.4)

2. Instructions to compile the Minimalist version
-------------------------------------------------
1. Open the `mybashrc.minimalist` and adjust the paths to the compulsory external libraries
2. Type: `source mybashrc.minimalist`
3. Compile the flameletSMOKE library: from the `libs/thermophysicalModels/basic` folder type `wmake`
4. Compile the steady-state solver: from the `solvers/flameletSimpleSMOKE` folder type `wmake`
5. Compile the unsteady solvers:
   - from the `solver/flameletPimpleSMOKE` folder type `wmake`
   - from the `solver/flameletPisoSMOKE` folder type `wmake` (available only for OpenFOAM 2.2, 2.3, and 2.4)

Run your first case
-------------------
The `cases` folder contains simple test cases (Sandia CO/H2/N2 turbulent jet flame).

1. The first operation to carry out is the generation of the lookup-table using the `OpenSMOKE++Suite` framework. This can be done in 2 steps:
  i. generation of a database of steady state, non-adiabatic laminar flamelets: from the `cases/lookUpTableGeneration/Sandia_COH2N2/flamelets` folder, after modifying the paths in the `input.dic` file, type `Run.sh`. In the `Output` folder the generated flamelets are available as XML files.
  ii. generation of the look-up table from the flamelets generated in 1a: from the `cases/lookUpTableGeneration/Sandia_COH2N2/library` folder, after modifying the paths in the `input.dic` file, type `Run.sh`. In the `OutputXML` folder the generated look-up tables are available as XML files.

2. Unsteady simulation: open the `cases/flameletPimpleSMOKE/Sandia_COH2N2` folder, build the mesh using the `blockMesh` utility, and run the case using the `flameletPimpleSMOKE` solver. 
   Even if you are interested in steady state conditions, we strongly suggest to always start with unsteady calculations to create a reasonable first-guess solution for the application of the steady state solver. 

3. Steady state simuation: open the `cases/flameletSimpleSMOKE/Sandia_COH2N2` folder, build the mesh using the `blockMesh` utility, and run the case using the `flameletSimpleSMOKE` solver. 

[1]: https://www.opensmokepp.polimi.it/
