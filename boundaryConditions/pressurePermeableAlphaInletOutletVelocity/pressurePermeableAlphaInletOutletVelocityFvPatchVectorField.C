/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "pressurePermeableAlphaInletOutletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressurePermeableAlphaInletOutletVelocityFvPatchVectorField::
pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    alphaName_("none"),
    alphaMin_(1.0)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1.0;
}


Foam::pressurePermeableAlphaInletOutletVelocityFvPatchVectorField::
pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
(
    const pressurePermeableAlphaInletOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    alphaName_(ptf.alphaName_),
    alphaMin_(ptf.alphaMin_)
{}


Foam::pressurePermeableAlphaInletOutletVelocityFvPatchVectorField::
pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    alphaName_(dict.lookupOrDefault<word>("alpha", "none")),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 1))
{
    patchType() = dict.lookupOrDefault<word>("patchType", word::null);
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1.0;
}


Foam::pressurePermeableAlphaInletOutletVelocityFvPatchVectorField::
pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
(
    const pressurePermeableAlphaInletOutletVelocityFvPatchVectorField& pivpvf
)
:
    mixedFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    alphaName_(pivpvf.alphaName_),
    alphaMin_(pivpvf.alphaMin_)
{}


Foam::pressurePermeableAlphaInletOutletVelocityFvPatchVectorField::
pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
(
    const pressurePermeableAlphaInletOutletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    alphaName_(pivpvf.alphaName_),
    alphaMin_(pivpvf.alphaMin_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::pressurePermeableAlphaInletOutletVelocityFvPatchVectorField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& phi = db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    const vectorField n(patch().nf());

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        refValue() = (phip/patch().magSf())*n;
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        refValue() = (phip/(rhop*patch().magSf()))*n;
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    valueFraction() = 1.0 - pos0(phip);

    if (alphaName_ != "none")
    {
        const scalarField& alphap =
            patch().lookupPatchField<volScalarField, scalar>(alphaName_);

        const scalarField alphaCut(pos(alphap - alphaMin_));
        valueFraction() = max(alphaCut, valueFraction());
        forAll (*this, faceI)
        {
            if (valueFraction()[faceI] == 1.0)
            {
                refValue()[faceI] = Zero;
            }
        }
    }

    mixedFvPatchVectorField::updateCoeffs();
}


void Foam::pressurePermeableAlphaInletOutletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    mixedFvPatchVectorField::write(os);
    writeEntryIfDifferent<word>(os,"phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os,"rho", "rho", rhoName_);
    writeEntryIfDifferent<word>(os,"alpha", "none", alphaName_);
    writeEntryIfDifferent<scalar>(os,"alphaMin", 1, alphaMin_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
::operator=
(
    const fvPatchField<vector>& pvf
)
{
    tmp<vectorField> n = patch().nf();

    fvPatchField<vector>::operator=
    (
        valueFraction()*(n()*(n() & pvf))
      + (1 - valueFraction())*pvf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        pressurePermeableAlphaInletOutletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
