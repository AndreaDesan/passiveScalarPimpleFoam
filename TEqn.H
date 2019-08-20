{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
//     ==
//      + fvOptions(T)
    );

    TEqn.relax();

    //fvOptions.constrain(TEqn);

    TEqn.solve();

    //fvOptions.correct(T);


    //bound T between MaxT and MinT
    /*
    dimensionedScalar MaxT ("MaxT", T.dimensions(), 1.0);
    dimensionedScalar MinT ("MinT", T.dimensions(), 1.0e-15);
    T = max(min(T,MaxT),MinT);
    */
}
