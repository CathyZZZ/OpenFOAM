{

Info << nl << "Solving density" << endl;


	forAll(rho, cellI) {
		rho[cellI] = pAir[cellI]/(R.value() * T[cellI]);
	}


 Info << "Rho min mean max     = "
         << Foam::gMin(rho.internalField()) << "    "
         << Foam::gAverage(rho.internalField()) << "    "
         << Foam::gMax(rho.internalField())
         << nl <<endl;

Info << nl << "Solving specific heat" << endl;


        forAll(cp, cellI)
{
              double tem=T[cellI]/1000;
              double A,B,C,D,E;

              if (T[cellI] <=1000)
{               A=33.066178;
                B=-11.363417;
                C=11.432816;
                D=-2.772874;
                E=-0.158558;    
              
                cp[cellI] =0.5 * 1000*(A + B * tem + C * pow(tem, 2) + D * pow(tem, 3) + E / pow(tem, 2));
        }        
              else
{               A=18.563083;
                B=12.257357;
                C=-2.859786;
                D=0.268238;
                E=1.977990;

                cp[cellI] =0.5 * 1000 * (A + B * tem + C * pow(tem,2) + D * pow(tem,3) + E / pow(tem,2));
        }
}


Info << nl << "Solving dynamic viscosity" << endl;
//Sutherland's formula

        forAll(muAir, cellI)
{ 
              double a,b,tem;
              a= 0.555 * 528.93 + 72;
              b= 0.555 * T[cellI] * 1.8+ 72;
              tem=T[cellI] *1.8/528.93;
              muAir[cellI]= 8.76e-6 * a / b * pow(tem, 3/2);
}
}
/*
Info << nl << "Solving diffusivity" << endl;

        forAll(D, cellI)
{
             double d, M,tem1,tem2,tem3;
             d=2.89e-10;
             M=2;
             tem1=2/(3* pow(M_PI,3/2 ) * pow(d,2));
             tem2= R.value() * T[cellI] /( pAir[cellI] * Na.value() );
             tem3= pow(R.value() * T[cellI] /M, 1/2);
             D[cellI]=tem1 * tem2 * tem3;
}*/
