# passiveScalarPimpleFoam

The repository contains the source code for the addition of a passive scalar transport to the incompressible pimpleFoam solver of OpenFoam. 

The passive scalar is indicated as T (for temperature) since heat transfer modelling is probably going to be the most common use for the solver. Neverthless, it can be used to track the evolution of any passive scalar which follows the following transport equation

```cpp
alphat = turbulence->nut()/Prt;
alphat.correctBoundaryConditions();

volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

fvScalarMatrix TEqn
(
    fvm::ddt(T)
     + fvm::div(phi, T)
     - fvm::laplacian(alphaEff, T)
);
```
As it can be inferred from the code above, a Simple Gradiend Diffusion Hypothesis is used to evaluate the turbulent flux of T. The turbulent diffusivity of T is evaluated using a constat turbulent Schmidt number assumption (the turbulent Schmidt number is actually indicadet as Prt in the code, since it corresponds to the turbulent Prandtl number when T represents the temperature). Similarly, the molecular Schmidt number is indicated with Pr within the code.

I have tested the solver with different Pr and Prt values for the PitzDaily case ($FOAM_TUTORIALS/incompressible/pimpleFoam/RAS/pitzDaily) and it gave reasonable (qualitative) results for different values of Pr and Prt; nevertheless, any independent assessment and any suggesion/feedback are greatly appreciated.

The repository also includes:
* A steady-state simpleFoam-based version of the solver 
* A test case based on the pitzDaily tutorial
