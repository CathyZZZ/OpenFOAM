{

        Info << nl << "Solving air scalars" << endl;
forAll(dp, i) {
        Info << "diameter for node " << i << " " << Nair[i][0] << " = " << dp[i].value() << " " <<endl;
}

 //       rhoTFlux = mesh.magSf()*fvc::snGrad(T);//-3*mu*rhoAir*rhoTFlux/4/T

volScalarField& N1 = Nair[0];
Info << nl << "N1[400]="<< N1[400] << nl << "N1[800]=" << N1[800] << nl << "N1[5000]"<< N1[5000] << endl; 

Info << nl << "T[400]="<< T[400] << nl << "T[800]=" << T[800] << nl << "T[5000]"<< T[5000] << endl; 

volScalarField term1 = S*0.0;
volScalarField term2 = S*0.0; //Foam::pow((2*sigma)/(pi*m1),0.5);
volScalarField term3 = S*0.0; //Foam::log(S);
volScalarField term4 = S*0.0;
volScalarField term5 = S*0.0; //Foam::exp(term4);
volScalarField term6 = S*0.0; //2*theta/(3*term3);
volScalarField term7 = S*0.0;
volScalarField SN1 = S*0.0; //kb*T/pSat;
volScalarField JkN1 = S*0.0; //SN1*v1*term1*term2*term5;

forAll(N1,cellI) {
    sigma[cellI] = (ACoeff.value()-BCoeff.value()*T[cellI])/1000;
    pSat[cellI] = pAir[cellI]*Foam::exp(CCoeff.value()-DCoeff.value()/T[cellI]);
    S[cellI] = N1[cellI]* kb * T[cellI]/pSat[cellI];
    theta[cellI] = s1.value()*sigma[cellI]/(T[cellI] * kb );
    term1[cellI] = Foam::pow((pSat[cellI]/( kb * T[cellI])),2);
    term2[cellI] = Foam::pow((2*sigma[cellI])/(pi*m1.value()),0.5);
    term3[cellI] = Foam::log(S[cellI]+0.01);
    term4[cellI] = theta[cellI] - (4*Foam::pow(theta[cellI],3))/(27*Foam::pow(term3[cellI],2));
    term5[cellI] = Foam::exp(term4[cellI]);
    term6[cellI] = 2*theta[cellI]/(3*term3[cellI]);
    SN1[cellI] = kb * T[cellI]/pSat[cellI];//N1 at saturated condition
    JkN1[cellI] = SN1[cellI]*v1.value()*term1[cellI]*term2[cellI]*term5[cellI];
    Jk[cellI] = S[cellI]*v1.value()*term1[cellI]*term2[cellI]*term5[cellI];
    kstar[cellI] = Foam::pow(term6[cellI],3);
    dpstar[cellI] = 4*sigma[cellI]*v1.value()/(kb * T[cellI] * term3[cellI]);
    vpstar[cellI] = pi/6*Foam::pow(dpstar[cellI],3);
 //   Info << "RRRRRRRRRRRR"<< endl;
//    Info << sigma[cellI]<<", "<< term6[cellI]<<", "<< vpstar[cellI]<<endl;
}
/*
Info << nl << "pAir[400]="<< pAir[400] << nl << "pAir[800]=" << pAir[800] << nl << "pAir[5000]"<< pAir[5000] << endl;
Info << nl << "sigma[400]="<< sigma[400] << nl << "sigma[800]=" << sigma[800] << nl << "sigma[5000]"<< sigma[5000] << endl;
Info << nl << "pSat[400]="<< pSat[400] << nl << "pSat[800]=" << pSat[800] << nl << "pSat[1500]"<< pSat[1500] <<nl<< "pSat[5000]="<< pSat[5000]<< endl;
Info <<" s1"<< s1<<nl<<"kb"<<kb <<"m1"<<m1 << endl;
Info << nl << "S[5000]="<< S[5000] << nl << "theta[5000]=" << theta[5000] << nl << "term1[5000]"<< term1[5000] <<nl<< "term2[5000]="<< term2[5000]<< endl;

Info << nl << "term3[5000]="<< term3[5000] << nl << "term4[5000]=" << term4[5000] << nl << "term5[5000]"<< term6[5000] <<nl<< "SN1[5000]="<< SN1[5000]<< endl;

Info << nl << "JkN1[5000]="<< JkN1[5000] << nl << "Jk[5000]=" << Jk[5000] << nl << "kstar[5000]"<< kstar[5000] <<nl<< "dpstar[5000]="<< dpstar[5000]<< endl;

Info << nl << "Jk[400]="<< Jk[400] << nl << "Jk[800]=" << Jk[800] << nl << "Jk[5000]"<< Jk[5000] << endl;
Info << nl << "kstar[400]="<< kstar[400] << nl << "kstar[800]=" << kstar[800] << nl << "kstar[5000]"<< kstar[5000] << endl;

Info << nl << "dpstar[400]="<< dpstar[400] << nl << "dpstar[800]=" << dpstar[800] << nl << "dpstar[5000]"<< dpstar[5000] << endl;
Info << nl << "vpstar[400]="<< vpstar[400] << nl << "vpstar[800]=" << vpstar[800] << nl << "vpstar[5000]"<< vpstar[5000] << endl;
*/
//collision frequency beta calculation
forAll(term7,cellI) {
    term7[cellI] = Foam::pow((3/( pi * 4)),1 / 6) * Foam::pow(( 6* kb * T[cellI] /rhoDisp.value()), 0.5);
}
Info << nl << "term7[400]="<< term7[400] << nl << "term7[800]=" << term7[800] << nl << "term7[1500]"<< term7[1500] << endl;
//    Info << "PPPPPPPPP"<< endl;
//    Info << "MW"<< MWDisp << endl;
//saturation concentration over particles of size i
    //Info << T<< endal;
 
for (int k =0; k<node; k++) {
    forAll( NSatk[k], cellI) {       
         NSatk[k][cellI] =  pSat[cellI]/ (kb *T[cellI]) *Foam::exp(4000 * sigma[cellI] * MWDisp.value()/(R.value() * T[cellI] * rhoDisp.value() * dp[k].value()));        
}

  //       Info <<" k=" <<k << endl;

//Info <<"NSatl["<<k<<"]="<< NSatk[k][5000] << endl;
/*
double aa= pSat[5000] / (kb * T[5000]) ;
double bb = 4 * sigma[5000] * MWDisp.value() / (R.value() * T[5000] * rhoDisp.value() * dp[k].value());
double cc = Foam::exp(bb);
double dd = aa * cc;
*/
/*
Info << nl << "NSatk["<<k<<"][5000]="<< NSatk[k][5000] << endl;
//Info << aa << nl << bb << nl << cc << dd <<endl;
Info <<"NSatk_total"<< pSat[5000] / (kb * T[5000]) * Foam::exp(4 * sigma[5000] * MWDisp.value() / (R.value() * T[5000] * rhoDisp.value() * dp[k].value())) <<endl;
//Info << "1="<< dp[k].value() << endl;*/
/*Info << "pSat=" <<pSat[5000] <<nl<<"kb="<< kb <<nl<<"T="<< T[5000] <<nl<<"sigma="<< sigma[5000] <<"MWDisp="<< MWDisp.value() <<nl<<"R="<< R.value() <<nl<<"rhoDisp="<< rhoDisp.value() <<nl<<"dp[k]"<< dp[k].value() << endl;
*/
}


scalar beta = 0;
 //   Info << "XXXXXXXx"<< endl;
//collision frequency kernels used for coagulation and growth
for (label j = 0; j < node; j++) {
 //   Info<<"j="<< j<< endl;
    for (label i = 0; i < node; i++) {
        forAll(Beta[j][i], cellI) {

             beta = Foam::pow(vp[j].value(), 1.0/3)+Foam::pow(vp[i].value(),1.0/3);
//            Info<<"qqqq"<< endl;
            Beta[j][i][cellI] = dimCorr.value() *term7[cellI] *Foam::pow(1/vp[i].value()+1/vp[j].value(),0.5)*Foam::pow(beta,2);
        }
//        double tem2= Foam::pow(1/vp[i].value()+1/vp[j].value(),0.5);
//        Info <<" temporary_term["<<j<<"]["<<i<< "]="<< beta << endl;
//        Info << "dimCorr="<<dimCorr<<endl;
//        Info <<" temporary_term2[" << j << "][" <<i << "]="<< tem2 << endl;
//        Info << "Beta["<<j<<"]["<<i<<"]="<< Beta[j][i][5000] << endl;
                                                                                                                    
//        Info << "Beta[" << i << "][ " << j <<"]="<< Beta[i][j][5000] << endl;
       // Info <<"size_ofBeta"<< Beta.size() << endl;
       // Info <<"size_ofBeta[0]"<< Beta[0].size() << endl;
    }
}
//Info<<"1111111111"<< endl;

//Info << N1<< endl;
volScalarField srcMonoNuc = 0.0 * Jk ;
//Info << Jk<< endl;
//Info << N1<< endl;
//Info<<srcMonoNuc<< endl;
volScalarField srcNucleation = 0.0 * Jk;
//    Info<<"22222222"<< endl;

volScalarField tem_growth = Nair[1];
//    Info<<"22222222"<< endl;



//SOURCE TERM CALCULATION
for (label k = 1; k<node; k++) {
    volScalarField& Nk = Nair[k];
    volScalarField srcCoagulation = 0.0 *Jk;
    volScalarField srcgrowth = 0.0*Jk;
    volScalarField addterm = 0.0 *Jk;
    volScalarField subterm = 0.0 *Jk2;
    volScalarField coagterm1 = 0.0 *Jk;//addterm
    volScalarField coagterm2 = 0.0 *Jk;//sub term
    volScalarField coagterm3 = 0.0 *Jk;
//    Info<<"22222222"<< endl;
//nucleation source terms

    if (nucleation) {
        forAll (Nzeta[k], cellI) {
            if (vpstar[cellI] <= v1.value()) {
              Nzeta[1][cellI] = vpstar[cellI] / vp[1].value();
            }
            else if (( vp[k-1].value() <= vpstar[cellI] && (vpstar[cellI] <= vp[k].value()))) {
            Nzeta[k][cellI] = vpstar[cellI] / vp[k].value();
            }
            else {
                Nzeta[k][cellI] = 0;
            }
            srcNucleation[cellI] = Jk[cellI] * Nzeta[k][cellI];
            srcMonoNuc[cellI] = -Jk[cellI] * kstar[k];//different from Robert's
       }

 //       Info << nl << "srcNucleation[400]="<< srcNucleation[400] << nl << "srcNucleation[800]=" << srcNucleation[800] << nl << "srcNucleation[1500]"<< srcNucleation[1500] << endl;
    }

    if (coagulation) {
        forAll (Nk, cellI) {
            for ( label i=1; i < node; i++) {
                coagterm2[cellI] += Beta[i][k][cellI] * Nair[i][cellI];
                for ( label j=1; j < node; j++) {
                coagterm1[cellI] += Xijk[i][j][k].value() * Beta[i][j][cellI] * Nair[i][cellI] * Nair[j][cellI];
                }
            }
        }
        forAll(srcCoagulation, cellI) {
            srcCoagulation[cellI] = 0.5 * coagterm1[cellI] - Nair[k][cellI] * coagterm2[cellI];
        }
   //     Info << nl << "srcCoagulation[400]="<< srcCoagulation[400] << nl << "srcCoagulation[800]=" << srcCoagulation[800] << nl << "srcCoagulation[1500]"<< srcCoagulation[1500] << endl;
    }


    if (growth) {
        if (k < Nair.size()) {
            //double oem_growth[cellI] = Nair[1][cellI];

            forAll ( Nk, cellI) {
                if (N1[cellI] > NSatk[k-1][cellI]) {
                    if (k==1) {
                        addterm[cellI] = 0;
                    }
                    else
                     
                        addterm[cellI] = v1.value() /(vp[k].value() - vp[k-1].value()) * Beta[0][k-1][cellI]* (N1[cellI] -NSatk[k-1][cellI]) * tem_growth[cellI];
                        //addterm[cellI] = v1.value()/(vp[k].value() - vp[k-1].value()) * Beta[0][k-1][cellI] (N1[cellI] -NSatk[k-1][cellI]) * Nair[k-1][cellI];
                }
                if (N1[cellI]<NSatk[k+1][cellI]) {
                    addterm[cellI] = -(v1.value() /( vp[k+1].value() - vp[k].value())) * Beta[0][k+1][cellI] *(N1[cellI] - NSatk[k+1][cellI]) * Nair[k+1][cellI];
                }
                if (N1[cellI] < NSatk[k][cellI]) {
                    subterm[cellI] = -(v1.value() / (vp[k].value() - vp[k-1].value())) * Beta[0][k][cellI] * (N1[cellI] - NSatk[k][cellI]) * Nair[k][cellI];
                }

                if (N1[cellI] > NSatk[k][cellI]) {
                    subterm[cellI] = (v1.value() / (vp[k+1].value() - vp[k].value())) * Beta[0][k][cellI] * (N1[cellI] - NSatk[k][cellI]) * Nair[k][cellI];
                }
                if (k != 1) {
                    tem_growth[cellI] = Nair[k-1][cellI]; 
                }         
                srcgrowth[cellI] = addterm[cellI] + subterm[cellI];
            }
        }
    
     //   Info << nl << "srcgrowth[400]="<< srcgrowth[400] << nl << "srcgrowth[800]=" << srcgrowth[800] << nl << "srcgrowth[1500]"<< srcgrowth[1500] << endl;
    }

    solve
        (
         fvm::div(phiAir, Nk, "div(phi,n)")
         - rho * srcCoagulation
         - rho * srcNucleation
         - rho * srcgrowth,
         mesh.solver("nPart")
        );
    Nk.max(Foam::doubleScalarSMALL);
}
}
