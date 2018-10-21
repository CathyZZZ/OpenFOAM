{
    Info << nl << "Solving air scalars" << endl;

forAll(dp, i) {
Info << "diameter for node " << i << " " << Nair[i][0] << " = " << dp[i].value() << " " <<endl;
}
	rhoTFlux = mesh.magSf()*fvc::snGrad(T);//-3*mu*rhoAir*rhoTFlux/4/T

// equations for dispersed phase represented by size nodes
	volScalarField& N1 = Nair[0];
	Info << nl << "N1 " << N1[0]  
	<< nl<< "N1 " << N1[472] 
	<< nl<< "N1 " << N1[1500] 
	<< nl<< "N1 " << N1[9999] << endl;
	Info << nl << "S " << S[0]  
	<< nl<< "S " << S[472] 
	<< nl<< "S " << S[1500] 
	<< nl<< "S " << S[9999] << endl;

//	sigma = (ACoeff-BCoeff*T)/1000;
//	pSat = pAir*Foam::exp(CCoeff-DCoeff/T);
//	S = N1*kb*T/pSat;
//	theta = s1*sigma/(T*kb);

	volScalarField term1 = S*0.0; //Foam::pow((pSat/(kb*T)),2);
	volScalarField term2 = S*0.0; //Foam::pow((2*sigma)/(pi*m1),0.5);
	volScalarField term3 = S*0.0; //Foam::log(S);
	volScalarField term4 = S*0.0;
	volScalarField term5 = S*0.0; //Foam::exp(term4);
	volScalarField term6 = S*0.0; //2*theta/(3*term3);
	volScalarField SN1 = S*0.0; //kb*T/pSat;
	volScalarField JkN1 = S*0.0; //SN1*v1*term1*term2*term5;



forAll(N1, cellI) {
	sigma[cellI] = (ACoeff.value()-BCoeff.value()*T[cellI])/1000;
	pSat[cellI] = pAir[cellI]*Foam::exp(CCoeff.value()-DCoeff.value()/T[cellI]);
	S[cellI] = N1[cellI]*kb.value()*T[cellI]/pSat[cellI];
	theta[cellI] = s1.value()*sigma[cellI]/(T[cellI]*kb.value());
	term1[cellI] = Foam::pow((pSat[cellI]/(kb.value()*T[cellI])),2);
	term2[cellI] = Foam::pow((2*sigma[cellI])/(pi*m1.value()),0.5);
	term3[cellI] = Foam::log(S[cellI]+0.01);
	term4[cellI] = theta[cellI] - (4*Foam::pow(theta[cellI],3))/(27*Foam::pow(term3[cellI],2));
	term5[cellI] = Foam::exp(term4[cellI]);
	term6[cellI] = 2*theta[cellI]/(3*term3[cellI]);

	SN1[cellI] = kb.value()*T[cellI]/pSat[cellI];
	JkN1[cellI] = SN1[cellI]*v1.value()*term1[cellI]*term2[cellI]*term5[cellI];

	Jk[cellI] = S[cellI]*v1.value()*term1[cellI]*term2[cellI]*term5[cellI];
	kstar[cellI] = Foam::pow(term6[cellI],3);
	dpstar[cellI] = 4*sigma[cellI]*v1.value()/(kb.value()*T[cellI]*term3[cellI]);
	vpstar[cellI] = pi/6*Foam::pow(dpstar[cellI],3);

}
//	Info<< nl << "N1 " << N1[0]
//	<< nl<< "N1 " << N1[472] 
//	<< nl<< "S " << S[0] 
//	<< nl << "S " << S[472] << nl << endl;
//	Info<< nl << "Jk " << Jk[0]
//	<< nl<< "Jk " << Jk[472] 
//	<< nl<< "theta " << theta[0] 
//	<< nl << "theta " << theta[472] << nl << endl;


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
	volScalarField srcMonoNuc = 0.0*Jk/N1;
	volScalarField srcNucleation = 0.0*Jk;


	// complicated calculation of equation [5] in Prakash et al. (2003) because pow(..,0.3333) does not give integer dimension exponents
	volScalarField term7 = 0.0*S; //sigma*m1*Foam::pow(d1,2.5);
    forAll(term7, cellI)
    {
       term7[cellI] = Foam::pow((3/(pi*4)),0.16666666666)*Foam::pow((6*kb.value()*T[cellI]/rhoDisp.value()),0.5);
    }

//saturation concentration over particles of size i
    for (label k = 0; k < Nair.size(); k++)
    {
		NSatk[k] = pSat/(kb*T)*Foam::exp(4*sigma*MWDisp/(R*T*rhoDisp*dp[k]));
	}

//collision frequency kernels used for coagulation and growth
    for (label k = 0; k < Nair.size(); k++)
    {
	    for (label j = 0; j < Nair.size(); j++)
    	{
			scalar gammaExpTerm = Foam::pow(vp[k].value(),0.3333333)+Foam::pow(vp[j].value(),0.3333333);
			
			Ngamma[k][j] = dimCorr*term7*Foam::pow(1/vp[j].value()+1/vp[k].value(),0.5)*Foam::pow(gammaExpTerm,2);
		}
	}



//for (label k = 0; k < Nair.size(); k++) {
//Info<< nl << k << ", vp[k] " << vp[k] << endl;
//}
//SOURCE TERM CALCULATION
    for (label k = 1; k < node-1; k++)
    {
		volScalarField& Nk = Nair[k];
//		volScalarField srcNucleation = 0.0*Jk;
		volScalarField srcCoagulation = 0.0*Jk;
		volScalarField addterm = 0.0*Jk;
		volScalarField subterm = 0.0*Jk2;
		volScalarField coagterm1 = 0.0*Jk;
		volScalarField coagterm2 = 0.0*Jk/N1;
		volScalarField coagterm3 = 0.0*Jk/N1;

		//nucleation source terms
		if (nucleation) {
			if (k>0) {
			forAll(Nzeta[k], cellI)
			{
				if (vpstar[cellI]<v1.value()) {
					Nzeta[1][cellI] = vpstar[cellI]/vp[1].value();
//					Info<< nl << "k=" << k << " vpstar<v1 " <<
//					nl << "k=" << k << "Cell=" << cellI << ", Nzeta[1][Cell] " << Nzeta[1][cellI] <<
//					nl << "k=" << k << "Cell=" << cellI << ", vpstar[Cell] " << vpstar[cellI] <<
//					nl << "k=" << k << "Cell=" << cellI << ", kstar[Cell] " << kstar[cellI] <<
//					nl << "k=" << k << "Cell=" << cellI << ", Jk[Cell] " << Jk[cellI] <<
//					nl << "k=" << k << "Cell=" << cellI << ", vp[k] " << vp[k].value() <<
//					nl << "k=" << k << "Cell=" << cellI << ", vp[k-1] " << vp[k-1].value() <<endl;
				}
				else {
					if ((vp[k-1].value() <= vpstar[cellI]) && (vpstar[cellI] < vp[k].value())) 
					{
					Nzeta[k][cellI] = vpstar[cellI]/vp[k].value();
					if (cellI == 1490) {
					Info<< nl << "k=" << k << " v(k-1) < vpstar<v1 < v(k)"  <<
					nl << "k=" << k << "Cell=" << cellI << ", Nzeta[1][Cell] " << Nzeta[1][cellI] <<
					nl << "k=" << k << "Cell=" << cellI << ", vpstar[Cell] " << vpstar[cellI] <<
					nl << "k=" << k << "Cell=" << cellI << ", kstar[Cell] " << kstar[cellI] <<
					nl << "k=" << k << "Cell=" << cellI << ", Jk[Cell] " << Jk[cellI] <<
					nl << "k=" << k << "Cell=" << cellI << ", vp[k] " << vp[k].value() <<
					nl << "k=" << k << "Cell=" << cellI << ", vp[k-1] " << vp[k-1].value() <<endl;
					}
					}
					else
					Nzeta[k][cellI] = 0;
				}
//				if (Nzeta[k][cellI]>1) 
//				{
//					Info<< nl << "k=" << k << "Cell=" << cellI << ", Nzeta[k][1] " << Nzeta[k][cellI] <<endl;
//					rc[k][cellI]=1; 
//				}
				srcNucleation[cellI] = Jk[cellI]*Nzeta[k][cellI];
				srcMonoNuc[cellI] = JkN1[cellI]*kstar[cellI];

//				if (srcMonoNuc[cellI] > 1)
//					srcMonoNuc[cellI]=0.01;
	
//			if (cellI==1) {
//				Info<< nl << k << ", Nzeta[k][1] " << Nzeta[k][cellI] <<endl;
//				Info<< nl << k << ", vpstar[1] " << vpstar[cellI] << nl <<endl;
//				Info<< nl << k << ", T[1] " << T[cellI] << nl <<endl;
//			}
//			if (cellI==2) {
//				Info<< nl << k << ", Nzeta[k][2] " << Nzeta[k][cellI] <<endl;
//				Info<< nl << k << ", vpstar[2] " << vpstar[cellI] << nl <<endl;
//			}
			
			}
			}
		}
		
	
		//growth source terms
		if (growth) {
		if (k < Nair.size()-1) 
		{
			    forAll(N1, cellI)
    			{
					//condensation of monomers on node k-1 --> increase of Nk, decrease of N1				(decrease of N_(k-1))
					if (k>1) 
					{
						if (N1[cellI]>NSatk[k-1][cellI]) 
						{

							addterm[cellI] = (v1.value()/(vp[k].value()-vp[k-1].value()))*Ngamma[0][k-1][cellI]*(N1[cellI]-NSatk[k-1][cellI])*Nair[k-1][cellI];
							t1Plus[cellI] = t1Plus[cellI] + Ngamma[0][k-1][cellI]*Nair[k-1][cellI];
							t1Minus[cellI] = t1Minus[cellI] - Ngamma[0][k-1][cellI]*NSatk[k-1][cellI]*Nair[k-1][cellI];
						}
					}

					//evaporation of monomers from node k+1 --> increase of Nk and N1						(decrease of N_(k-1))
					if (N1[cellI]<NSatk[k+1][cellI]) 
					{
						addterm[cellI] = -(v1.value()/(vp[k+1].value()-vp[k].value()))*Ngamma[0][k+1][cellI]*(N1[cellI]-NSatk[k+1][cellI])*Nair[k+1][cellI];
						t2Plus[cellI] = t2Plus[cellI] + Ngamma[0][k+1][cellI]*Nair[k+1][cellI];
						t2Minus[cellI] = t2Minus[cellI] - Ngamma[0][k+1][cellI]*NSatk[k+1][cellI]*Nair[k+1][cellI];
					}

					//evaporation of monomers from node k --> decrease of Nk, increase of N1 				(increase of N_(k-1))
					if (N1[cellI]<NSatk[k][cellI]) 
					{
						subterm[cellI] = -(v1.value()/(vp[k].value()-vp[k-1].value()))*Ngamma[0][k][cellI]*(N1[cellI]-NSatk[k][cellI]);
						t3Plus[cellI] = t3Plus[cellI] + Ngamma[0][k][cellI]*Nair[k][cellI];
						t3Minus[cellI] = t3Minus[cellI] - Ngamma[0][k][cellI]*NSatk[k][cellI]*Nair[k][cellI];
					}

					//condensation of monomers on node k --> decrease of Nk and N1
					if (N1[cellI]>NSatk[k][cellI]) 
					{
						subterm[cellI] = (v1.value()/(vp[k+1].value()-vp[k].value()))*Ngamma[0][k][cellI]*(N1[cellI]-NSatk[k][cellI]);
						t4Plus[cellI] = t4Plus[cellI] + Ngamma[0][k][cellI]*Nair[k][cellI];
						t4Minus[cellI] = t4Minus[cellI] - Ngamma[0][k][cellI]*NSatk[k][cellI]*Nair[k][cellI];
					}

				srcMonoPlus[cellI] = t1Plus[cellI]+t2Plus[cellI]+t3Plus[cellI]+t4Plus[cellI];
				srcMonoMinus[cellI] = t1Minus[cellI]+t2Minus[cellI]+t3Minus[cellI]+t4Minus[cellI];
				}
		}
		}


		//coagulation source terms
		if (coagulation) {

			forAll(Nk, cellI) 
			{
				for (label i = 1; i < node-1; i++) {
					coagterm2[cellI] += Ngamma[k][i][cellI]*Nair[i][cellI];
					for (label j = 1; j <= k; j++) {
//						Info << "Xijk " << Xijk[i][j][k].value() <<endl;
						coagterm1[cellI] += Xijk[i][j][k].value()*Ngamma[i][j][cellI]*Nair[i][cellI]*Nair[j][cellI];
					}
				}

			}
			forAll(srcCoagulation, cellI) 
			{
				srcCoagulation[cellI] = 0.5*coagterm1[cellI];
			}

		}

	    if (k > 0)
	    {
       	solve
		(
        	fvm::div(phiAir, Nk,"div(phi,n)")
//              - fvm::laplacian(rhoAir*diffPart, Nk, "laplacian(diff,n)")
//              + fvm::div(-kEPart*rhoAir*rhoEFlux*CPPPart*indx, Nk,"div(phi,n)")
			  	+ fvm::Sp(rhoAir*subterm, Nk)
			  	+ fvm::Sp(rhoAir*coagterm2, Nk)
//			  	- fvm::Sp(rhoAir*coagterm3, Nk)
				- rhoAir*srcNucleation
			  	- rhoAir*srcCoagulation
     	      	- rhoAir*addterm,
                mesh.solver("nPart")
		);

        Nk.max(Foam::doubleScalarSMALL);
		}

	}


	//Monomer Balance
	if (monomerBal)
	{
	volScalarField N1old = N1;

	fvScalarMatrix monoBalEqn(
	fvm::div(phiAir, N1,"div(phi,n)") 
	+ fvm::Sp(rhoAir*srcMonoNuc, N1)
	);
    // under-relax and solve
    monoBalEqn.relax(0.99);
	monoBalEqn.solve(mesh.solver("nPart"));
//    monoBalEqn.solve();


//       	solve
//		(
//		fvm::div(phiAir, N1,"div(phi,n)")
//     	      	+ fvm::Sp(rhoAir*srcMonoNuc, N1),
////     	      	+ fvm::Sp(rhoAir*srcMonoPlus, N1)
// //    	      	+ rhoAir*srcMonoMinus,
//                mesh.solver("nPart")
//		);


		forAll(N1, cellI)
    	{
        	if (srcMonoNuc[cellI]==0.01)
			{	
				Info<< nl << "Cell=" << cellI << ", srcMonoNuc " << srcMonoNuc[cellI]
				<< nl << "Cell=" << cellI << ", N1 " << N1old[cellI] 
				<< nl << "Cell=" << cellI << ", N1 " << N1[cellI] <<endl;
				Info<< nl << "Cell=" << cellI << ", N1 " << N1[cellI] <<endl;
			}
		}
	}



}

