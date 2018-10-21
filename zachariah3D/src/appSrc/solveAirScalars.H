{
    Info << nl << "Solving air scalars" << endl;

	rhoTFlux = mesh.magSf()*fvc::snGrad(T);//-3*mu*rhoAir*rhoTFlux/4/T

// equations for dispersed phase represented by size nodes

	volScalarField& N1 = Nair[0];
	sigma = (ACoeff-BCoeff*T)/1000;
	pSat = pAir*exp(CCoeff-DCoeff/T);
	S = N1*kb*T/pSat;
	theta = s1*sigma/(T*kb);

	Info<< nl << "N1 " << N1[0]
	<< nl<< "N1 " << N1[100] 
	<< nl<< "N1 " << N1[1000] 
	<< nl<< "S " << S[0] 
	<< nl << "S " << S[100] 
	<< nl << "S " << S[1000] 
	<< nl << "pSat " << pSat[1] 
	<< nl << "pSat " << pSat[100] 
	<< nl << "pSat " << pSat[1000] << nl <<endl;

	volScalarField term1 = Foam::pow((pSat/(kb*T)),2);
	volScalarField term2 = Foam::pow((2*sigma)/(pi*m1),0.5);
	volScalarField term3 = Foam::log(S+0.000000001);
	volScalarField term4 = theta - (4*Foam::pow(theta,3))/(27*Foam::pow(term3,2));
	volScalarField term5 = Foam::exp(term4);
	volScalarField term6 = 2*theta/(3*term3);
	Jk = S*v1*term1*term2*term5;
	kstar = Foam::pow(term6,3);
	dpstar = 4*sigma*v1/(kb*T*term3);
	vpstar = pi/6*Foam::pow(dpstar,3);
	volScalarField t1 = 0.0*Jk;
	volScalarField t2 = 0.0*Jk;
	volScalarField t3 = 0.0*Jk;
	volScalarField t4 = 0.0*Jk;


	// complicated calculation of equation [5] in Prakash et al. (2003) because pow(..,0.3333) does not give integer dimension exponents
	volScalarField term7 = 0.0*S; //sigma*m1*Foam::pow(d1,2.5);
    forAll(term7, cellI)
    {
       term7[cellI] = Foam::pow((3/(pi*4)),0.16666666)*Foam::pow((6*kb.value()*T[cellI]/rhoDisp.value()),0.5);
    }

    for (label i = 0; i < Nair.size(); i++)
    {
        Info << " Nair[i].name() = " << Nair[i].name() << nl;

            volScalarField& Ni = Nair[i];
			volScalarField N1sati = 0.0*Nair[i];
			volScalarField N1satiPlusOne = 0.0*Nair[i];
			volScalarField N1satiMinusOne = 0.0*Nair[i];
	    	volScalarField NiPlusOne = 0.0*Nair[i];
            volScalarField NiMinusOne = 0.0*Nair[i];
			volScalarField gammai = Ngamma[i];
			volScalarField gammaiPlusOne = 0.0*Ngamma[i];
			volScalarField gammaiMinusOne = 0.0*Ngamma[i];
			volScalarField zetai = Nzeta[i];
			volScalarField addterm = 0.0*Jk;
			volScalarField subterm = 0.0*Jk2;
            dimensionedScalar vi = vp[i];
            dimensionedScalar viMinusOne = 0.0*vp[i];
            dimensionedScalar viPlusOne = 0.0*vp[i];
			dimensionedScalar di = dp[i];
			dimensionedScalar diPlusOne = 0.0*dp[i];
			dimensionedScalar diMinusOne = 0.0*dp[i];
		
	    if (i < Nair.size()-1)
	    {	
			NiPlusOne = Nair[i+1];
			viPlusOne = vp[i+1];
			diPlusOne = dp[i+1];
			//saturation concentration
			N1sati = pSat/(kb*T)*Foam::exp(4*sigma*MWDisp/(R*T*rhoDisp*di));
			N1satiPlusOne = pSat/(kb*T)*Foam::exp(4*sigma*MWDisp/(R*T*rhoDisp*diPlusOne));
			//collision kernels
			scalar gammaExpTerm = Foam::pow(v1.value(),0.3333333)+Foam::pow(vi.value(),0.3333333);
			scalar gammaExpTerm1 = Foam::pow(v1.value(),0.3333333)+Foam::pow(viPlusOne.value(),0.3333333);

			gammai = dimCorr*term7*Foam::pow(1/vi.value()+1/v1.value(),0.5)*Foam::pow(gammaExpTerm,2);
			gammaiPlusOne = dimCorr*term7*Foam::pow(1/vi.value()+1/v1.value(),0.5)*Foam::pow(gammaExpTerm1,2);
	    }


	    if (i > 1)
	    {	
		    NiMinusOne = Nair[i-1];
			viMinusOne = vp[i-1];
			diMinusOne = dp[i-1];
			N1satiMinusOne = pSat/(kb*T)*Foam::exp(4*sigma*MWDisp/(R*T*rhoDisp*diMinusOne));
			scalar gammaExpTerm2 = Foam::pow(v1.value(),0.3333333)+Foam::pow(viMinusOne.value(),0.3333333);
			gammaiMinusOne = dimCorr*term7*Foam::pow(1/vi.value()+1/v1.value(),0.5)*Foam::pow(gammaExpTerm2,2);
		}
		//growth source terms
		if (i < Nair.size()-1) 
		{
			//growth of node k due to condensation of monomers on node k-1
			if (N1>N1satiMinusOne) 
			{
				if (i>1) 
				{
					addterm = (v1/(vi-viMinusOne))*gammaiMinusOne*(N1-N1satiMinusOne)*NiMinusOne;
					t1 = t1 + gammaiMinusOne*(N1-N1satiMinusOne)*NiMinusOne;
				}
			}
			//growth of node k due to evaporation of monomers from node k+1
			if (N1<N1satiPlusOne) 
			{
					addterm = -(v1/(viPlusOne-vi))*gammaiPlusOne*(N1-N1satiPlusOne)*NiPlusOne;
					t2 = t2 + gammaiPlusOne*(-N1+N1satiPlusOne)*NiPlusOne;
			}
			//evaporation of monomers from node k --> increase of Nk
			if (N1<N1sati) 
			{
					subterm = -(v1/(vi-viMinusOne))*gammai*(N1-N1sati);
					t3 = t3 + gammai*(-N1+N1sati)*Ni;
			}
			//condensation of monomers on node k --> decrease of Nk
			if (N1>N1sati) 
			{
					subterm = (v1/(viPlusOne-vi))*gammai*(N1-N1sati);
					t4 = t4 + gammai*(N1-N1sati)*Ni;
			}
		}

		if (addterm>N1)
		{
			addterm = 0;
		}	

		if (subterm>N1)
		{
			subterm = 0;
		}	

	    if (i > 0)
	    {
       	solve
		(
        	fvm::div(phiAir, Ni,"div(phi,n)")
//              - fvm::laplacian(rhoAir*diffPart, Ni, "laplacian(diff,n)")
//              + fvm::div(-kEPart*rhoAir*rhoEFlux*CPPPart*indx, Ni,"div(phi,n)")
			  + fvm::Sp(rhoAir*subterm, Ni)
     	      - rhoAir*addterm,
                mesh.solver("nPart")
		);

        Ni.max(Foam::doubleScalarSMALL);
		}
	}

//	//Monomer Balance
//	if (monomerBal)
//	{
//       	solve
//		(
//		fvm::div(phiAir, N1,"div(phi,n)")
//     	      - rhoAir*(Jk*kstar - t1 - t4 + t3 + t2),
//                mesh.solver("nPart")
//		);
//	}

}

