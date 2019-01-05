{
   double f,tem,cc;

    Info << nl << "Solving Flux" << endl;

   Info << nl << "Diffusion Coefficient" << endl;

	forAll(Nair, cellI)
{
        Nair[cellI] = pAir[cellI] / (BoltzConst.value() * T[cellI]);
}




        forAll(lambda, cellI)
{
        lambda[cellI] =1 /( std::sqrt(2) * M_PI * pow(dAir.value(),2) * Nair[cellI]);
//	std::cout << lambda[cellI] << std::endl;
}


        forAll(kn, cellI)
{
        kn[cellI] = 2 * lambda[cellI] * pow(dp.value(), -1);
//        std::cout << kn[cellI] << std::endl;
}


        forAll(D, cellI)
{

	if(kn[cellI] < 0.1) {
        	f = 3 * M_PI * muAir[cellI] * dp.value();
	} else if(kn[cellI] > 10 ) {        
//		tem = Foam::sqrt( 2 * M_PI * BoltzConst.value() * T[cellI] * Na.value() / 58.5);
		tem = Foam::sqrt( 2 * M_PI * BoltzConst.value() * T[cellI]  * Na.value() / 2 );
//		std::cout << tem << std::endl;  
		f =(2.0 / 3) * pow(dp.value(),2) * rho[cellI] * tem  * (1 + M_PI * 0.9 / 8 );
//                std::cout << f << std::endl;
	} else {
		cc = 1 + kn[cellI] * (1.257 + 0.400 * Foam::exp(- 1.1 / kn[cellI])); 
		f = 3 * M_PI * muAir[cellI] * dp.value() / cc;
	}
        D[cellI] = BoltzConst.value() * T[cellI] / f;
// 	std::cout << D[cellI] << std::endl;


}










    fvScalarMatrix NEqn
    (
         fvm::div(phiAir, N, "div(phiAir,N)")
      - fvm::laplacian(rho * D, N, "laplacian(rho*D,N)")

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
