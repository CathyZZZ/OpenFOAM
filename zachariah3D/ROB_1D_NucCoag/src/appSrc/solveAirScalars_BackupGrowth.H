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
	<< nl<< "N1 " << N1[4665] 
	<< nl<< "S " << S[0] 
	<< nl << "S " << S[4665] << nl << endl;

	volScalarField term1 = Foam::pow((pSat/(kb*T)),2);
	volScalarField term2 = Foam::pow((2*sigma)/(pi*m1),0.5);
	volScalarField term3 = Foam::log(S+0.00001);
	volScalarField term4 = theta - (4*Foam::pow(theta,3))/(27*Foam::pow(term3,2));
	volScalarField term5 = Foam::exp(term4);
	volScalarField term6 = 2*theta/(3*term3);

	Jk = S*v1*term1*term2*term5;
	kstar = Foam::pow(term6,3);
	dpstar = 4*sigma*v1/(kb*T*term3);
	vpstar = pi/6*Foam::pow(dpstar,3);
	volScalarField t1Plus = 0.0*Jk/N1;
	volScalarField t1Minus = 0.0*Jk;
	volScalarField t2Plus = 0.0*Jk/N1;
	volScalarField t2Minus = 0.0*Jk;
	volScalarField t3Plus = 0.0*Jk/N1;
	volScalarField t3Minus = 0.0*Jk;
	volScalarField t4Plus = 0.0*Jk/N1;
	volScalarField t4Minus = 0.0*Jk;
	volScalarField srcMonoPlus = 0.0*Jk/N1;
	volScalarField srcMonoMinus = 0.0*Jk;
	volScalarField srcMonoNuc = 0.0*Jk;

	// complicated calculation of equation [5] in Prakash et al. (2003) because pow(..,0.3333) does not give integer dimension exponents
	volScalarField term7 = 0.0*S; //sigma*m1*Foam::pow(d1,2.5);
    forAll(term7, cellI)
    {
       term7[cellI] = Foam::pow((3/(pi*4)),0.16666666)*Foam::pow((6*kb.value()*T[cellI]/rhoDisp.value()),0.5);
    }

    for (label i = 0; i < Nair.size(); i++)
    {
			//saturation concentration
			NSati[i] = pSat/(kb*T)*Foam::exp(4*sigma*MWDisp/(R*T*rhoDisp*dp[i]));
			//collision kernels
			scalar gammaExpTerm = Foam::pow(v1.value(),0.3333333)+Foam::pow(vp[i].value(),0.3333333);
			Ngamma[i] = dimCorr*term7*Foam::pow(1/vp[i].value()+1/v1.value(),0.5)*Foam::pow(gammaExpTerm,2);
	}

    for (label i = 1; i < Nair.size(); i++)
    {

            volScalarField& Ni = Nair[i];
			volScalarField addterm = 0.0*Jk;
			volScalarField subterm = 0.0*Jk2;


		
		//growth source terms
		if (i < Nair.size()-1) 
		{
			    forAll(N1, cellI)
    			{
					//condensation of monomers on node k-1 --> increase of Nk, decrease of N1				(decrease of N_(k-1))
					if (i>1) 
					{
						if (N1[cellI]>NSati[i-1][cellI]) 
						{

							addterm[cellI] = (v1.value()/(vp[i].value()-vp[i-1].value()))*Ngamma[i-1][cellI]*(N1[cellI]-NSati[i-1][cellI])*Nair[i-1][cellI];
							t1Plus[cellI] = t1Plus[cellI] + Ngamma[i-1][cellI]*Nair[i-1][cellI];
							t1Minus[cellI] = t1Minus[cellI] - Ngamma[i-1][cellI]*NSati[i-1][cellI]*Nair[i-1][cellI];
						}
					}

					//evaporation of monomers from node k+1 --> increase of Nk and N1						(decrease of N_(k-1))
					if (N1[cellI]<NSati[i+1][cellI]) 
					{
						addterm[cellI] = -(v1.value()/(vp[i+1].value()-vp[i].value()))*Ngamma[i+1][cellI]*(N1[cellI]-NSati[i+1][cellI])*Nair[i+1][cellI];
						t2Plus[cellI] = t2Plus[cellI] + Ngamma[i+1][cellI]*Nair[i+1][cellI];
						t2Minus[cellI] = t2Minus[cellI] - Ngamma[i+1][cellI]*NSati[i+1][cellI]*Nair[i+1][cellI];
					}

					//evaporation of monomers from node k --> decrease of Nk, increase of N1 				(increase of N_(k-1))
					if (N1[cellI]<NSati[i][cellI]) 
					{
						subterm[cellI] = -(v1.value()/(vp[i].value()-vp[i-1].value()))*Ngamma[i][cellI]*(N1[cellI]-NSati[i][cellI]);
						t3Plus[cellI] = t3Plus[cellI] + Ngamma[i][cellI]*Nair[i][cellI];
						t3Minus[cellI] = t3Minus[cellI] - Ngamma[i][cellI]*NSati[i][cellI]*Nair[i][cellI];
					}

					//condensation of monomers on node k --> decrease of Nk and N1
					if (N1[cellI]>NSati[i][cellI]) 
					{
						subterm[cellI] = (v1.value()/(vp[i+1].value()-vp[i].value()))*Ngamma[i][cellI]*(N1[cellI]-NSati[i][cellI]);
						t4Plus[cellI] = t4Plus[cellI] + Ngamma[i][cellI]*Nair[i][cellI];
						t4Minus[cellI] = t4Minus[cellI] - Ngamma[i][cellI]*NSati[i][cellI]*Nair[i][cellI];
					}

//		srcMonoNuc[cellI] = Jk[cellI]*kstar[cellI];
		srcMonoPlus[cellI] = t1Plus[cellI]+t2Plus[cellI]+t3Plus[cellI]+t4Plus[cellI];
		srcMonoMinus[cellI] = t1Minus[cellI]+t2Minus[cellI]+t3Minus[cellI]+t4Minus[cellI];
				}
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


	//Monomer Balance
	if (monomerBal)
	{
       	solve
		(
		fvm::div(phiAir, N1,"div(phi,n)")
//     	      	+ rhoAir*srcMonoNuc
     	      	+ fvm::Sp(rhoAir*srcMonoPlus, N1)
     	      	+ rhoAir*srcMonoMinus,
                mesh.solver("nPart")
		);


//		forAll(N1, cellI)
//    	{
//        	if (N1[cellI]<0)
//			{
//            	N1[cellI] = 0;
//			}
//		}
	}



}

