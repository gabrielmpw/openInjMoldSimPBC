
#ifndef myConstrainPressure_H
#define myConstrainPressure_H

#include "boundaryConditions/prghPermeableAlphaTotalPressure/updateableSnGrad.H"
template<class TargetType, class Type>
inline const TargetType* myIsA(const Type& t)
{
    const Type* p = &t;
    return dynamic_cast<const TargetType*>(p);
}

template<class RhoType, class RAUType, class MRFType>
void myConstrainPressure
(
    volScalarField& p,
    const RhoType& rho,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rhorAU,
    const MRFType& MRF
)
{
    const fvMesh& mesh = p.mesh();

    volScalarField::Boundary& pBf = p.boundaryFieldRef();

    const volVectorField::Boundary& UBf = U.boundaryField();
    const surfaceScalarField::Boundary& phiHbyABf =
        phiHbyA.boundaryField();
    const typename RAUType::Boundary& rhorAUBf = rhorAU.boundaryField();
    const surfaceVectorField::Boundary& SfBf = mesh.Sf().boundaryField();
    const surfaceScalarField::Boundary& magSfBf =
        mesh.magSf().boundaryField();

    forAll(pBf, patchi)
    {
        typedef updateablePatchTypes::updateableSnGrad snGradType;
        const auto* snGradPtr = myIsA<snGradType>(pBf[patchi]);
        
        if (snGradPtr)
        {
            const_cast<snGradType&>(*snGradPtr).updateSnGrad
            (
                (
                    phiHbyABf[patchi]
                  - rho.boundaryField()[patchi]
                   *MRF.relative(SfBf[patchi] & UBf[patchi], patchi)
                )
               /(magSfBf[patchi]*rhorAUBf[patchi])
            );
        }
    }
}


#endif

