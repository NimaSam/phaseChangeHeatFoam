//Info << "solve pEqn" << endl;
{

    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("Dp", fvc::interpolate(rAU));

   volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (fvc::interpolate(HbyA) & mesh.Sf())
       + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );
    adjustPhi(phiHbyA, U, p_rgh);
    phi = phiHbyA;

    surfaceScalarField phig
    (
        (
           fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          //  fvc::interpolate(interface.sigmaK())*fvc::snGrad(interface.smoothAlpha())
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf()
    );

    phiHbyA += phig;
/*
    // Update the fixedFluxPressure BCs to ensure flux consistency
    setSnGrad<fixedFluxPressureFvPatchScalarField>
      (
          p_rgh.boundaryField(),
          (
              phiHbyA.boundaryField()
            - (mesh.Sf().boundaryField() & U.boundaryField())
          )/(mesh.magSf().boundaryField()*rAUf.boundaryField())
      );
*/
    Pair<tmp<volScalarField> > vDotP = twoPhaseProperties->vDotP();
    const volScalarField& vDotcP = vDotP[0]();
    const volScalarField& vDotvP = vDotP[1]();

//#include "HelmholtzEqn.H"



    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
           fvc::div(phiHbyA) - fvm::laplacian(rAUf, p_rgh)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));

        if (phaseChange){

        		  solve
        		  (
        				  p_rghEqn == fvm::Sp(vDotcP-vDotvP,p_rgh) - (vDotcP-vDotvP)*(pSat-rho*gh) ,
        				  //p_rghEqn == vDotcP-vDotvP
        				 // p_rghEqn == fvm::Sp(vDotPs,p_rgh) - (vDotPs)*(pSat-rho*gh),
        				  mesh.solver(p_rgh.select(pimple.finalInnerIter()))
        		  );
        }
        else
        {
                solve
                (
                		p_rghEqn,
                		mesh.solver(p_rgh.select(pimple.finalInnerIter()))
                );
        }
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + p_rghEqn.flux();

            U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    Pair<tmp<volScalarField> > mDot = twoPhaseProperties->mDotP();
       const volScalarField& mDotcP = mDot[0]();
       const volScalarField& mDotvP = mDot[1]();

       mDotP = (mDotcP-mDotvP)*(p-pSat);
}


