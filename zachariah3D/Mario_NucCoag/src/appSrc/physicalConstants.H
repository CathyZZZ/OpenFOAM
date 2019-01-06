/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright 2007 Hrvoje Jasak.  All rights reserved
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Namespace
    physicalConstants

Description
    Physical constants.

\*---------------------------------------------------------------------------*/

#ifndef physicalConstants_H
#define physicalConstants_H

#include "dimensionedTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace physicalConstant
{
    // Boltzmann constant [(m2.kg/(s2.K1))]
    // Google: 6.02214179e23
    const dimensionedScalar kb
    (
        "kb",
        dimensionSet(1, 2, -2, -1, 0, 0, 0),
        1.3806488e-23
    );


    // Vacuum dielectric constant (F/m)
    // NIST:  8.314472
    const dimensionedScalar epsilon0
    (
        "epsilon0",
        dimensionSet(-1, -3, 4, 0, 0, 2, 0),
        8.85418e-12
    );

    // Faraday's constant (charge of electron) [C mol^{-1}];  C = A s
    // NIST:  96485.3399
    const dimensionedScalar e
    (
        "e",
        dimensionSet(0, 0, 1, 0, 0, 1, 0),
	1.60217e-19
    );

  

    // Standard atmosphere pressure [Pa]
    // BIPM 10th Conferance Generale des Poids et Mesures (resolution 4)
    const dimensionedScalar pAtm
    (
        "pAtm",
        dimensionSet(1, -1, -2, 0, 0, 0, 0),
        1.01325e5
    );

    // Avogadro number [1/mol]
	// Wikipedia: 6.022140857e23
    const dimensionedScalar NA
    (
        "NA",
        dimensionSet(0, 0, 0, 0, -1, 0, 0),
        6.022140857e23
    );

    // universal gas constant
	// Wikipedia: 8.3144598 [kg*m^2/s^2/mol/K]
    const dimensionedScalar R
    (
        "R",
        dimensionSet(1, 2, -2, -1, -1, 0, 0),
        8.3144598
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
