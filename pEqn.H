{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );

    surfaceScalarField phig
    (
        (
            interface.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;

    // Update the fixedFluxPressure BCs to ensure flux consistency
    
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - (mesh.Sf().boundaryField() & U.boundaryField())
        )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
    );
    myConstrainPressure(p_rgh, geometricOneField(), U, phiHbyA, rAUf, MRF);
    tmp<fvScalarMatrix> p_rghEqnComp1;
    tmp<fvScalarMatrix> p_rghEqnComp2;

    if (pimple.transonic())
    {
        surfaceScalarField phid1("phid1", fvc::interpolate(psi1)*phi);
        surfaceScalarField phid2("phid2", fvc::interpolate(psi2)*phi);

        p_rghEqnComp1 =
            fvc::ddt(rho1) + fvc::div(phi, rho1) - fvc::Sp(fvc::div(phi), rho1)
          + correction
            (
                psi1*fvm::ddt(p_rgh)
              + fvm::div(phid1, p_rgh) - fvm::Sp(fvc::div(phid1), p_rgh)
            );
        deleteDemandDrivenData(p_rghEqnComp1.ref().faceFluxCorrectionPtr());
        p_rghEqnComp1.ref().relax();

        p_rghEqnComp2 =
            fvc::ddt(rho2) + fvc::div(phi, rho2) - fvc::Sp(fvc::div(phi), rho2)
          + correction
            (
                psi2*fvm::ddt(p_rgh)
              + fvm::div(phid2, p_rgh) - fvm::Sp(fvc::div(phid2), p_rgh)
            );
        deleteDemandDrivenData(p_rghEqnComp2.ref().faceFluxCorrectionPtr());
        p_rghEqnComp2.ref().relax();
    }
    else
    {
        p_rghEqnComp1 =
            fvc::ddt(rho1) + psi1*correction(fvm::ddt(p_rgh))
          + fvc::div(phi, rho1) - fvc::Sp(fvc::div(phi), rho1);

        p_rghEqnComp2 =
            fvc::ddt(rho2) + psi2*correction(fvm::ddt(p_rgh))
          + fvc::div(phi, rho2) - fvc::Sp(fvc::div(phi), rho2);
    }

    // Cache p_rgh prior to solve for density update
    volScalarField p_rgh_0(p_rgh);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh)
        );

        solve
        (
            (
                (max(alpha1, scalar(0))/rho1)*p_rghEqnComp1()
              + (max(alpha2, scalar(0))/rho2)*p_rghEqnComp2()
            )
          + p_rghEqnIncomp
        );

        if (pimple.finalNonOrthogonalIter())
        {
            volScalarField pCalc = p_rgh + (alpha1*rho1 + alpha2*rho2)*gh;
            p = max(pCalc + pAux, pMin); // pActual = pCalc + pAux, p is what PVT gets
            pAux = 0.5*min(pCalc + pAux - pMin, 0.0*pMin) + 0.5*pAux;//Kristjan: transfer sub-pMin to pAux and use in UEqn
            p_rgh = p - (alpha1*rho1 + alpha2*rho2)*gh;

            dgdt =
            (
                pos(alpha2)*(p_rghEqnComp2 & p_rgh)/rho2
              - pos(alpha1)*(p_rghEqnComp1 & p_rgh)/rho1
            );

            phi = phiHbyA + p_rghEqnIncomp.flux();

            U = HbyA
              + rAU*fvc::reconstruct((phig + p_rghEqnIncomp.flux())/rAUf);
            U.correctBoundaryConditions();
        }
    }

    // p = max(p_rgh + (alpha1*rho1 + alpha2*rho2)*gh, pMin);

    // Update densities from change in p_rgh
    rho1 += psi1*(p_rgh - p_rgh_0);
    rho2 += psi2*(p_rgh - p_rgh_0);

    rho = alpha1*rho1 + alpha2*rho2;

    K = 0.5*magSqr(U);

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;
}
