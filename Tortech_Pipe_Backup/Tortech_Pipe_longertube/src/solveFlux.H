{

    Info << nl << "Solving Flux" << endl;

    fvScalarMatrix NEqn
    (
         fvm::div(phiAir, N, "div(phiAir,N)")
      - fvm::laplacian(rho*D, N, "laplacian(rho*D,N)")

    );

    NEqn.solve();


    Info << "N min mean max     = "
         << Foam::gMin(N.internalField()) << "    " 
         << Foam::gAverage(N.internalField()) << "    "
         << Foam::gMax(N.internalField())
         << nl <<endl;

    Info<< " sum NEqn.source() =   " << gSum(NEqn.source()) << nl;
    Info<< endl;
}
