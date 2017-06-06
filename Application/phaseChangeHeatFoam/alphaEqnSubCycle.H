//Info << "solve alphaEqn" << endl;
rhoPhi=0.0*rhoPhi;
{
	Pair<tmp<volScalarField> > mDot = twoPhaseProperties->mDotAlphal();
				const volScalarField& mDotcAlpha= mDot[0]();
				const volScalarField& mDotvAlpha= mDot[1]();

				mDotAlpha = mDotvAlpha*alpha1+mDotcAlpha*(1.0-alpha1);


    const dictionary& pimpleDict = pimple.dict();

    label nAlphaCorr(readLabel(pimpleDict.lookup("nAlphaCorr")));

    label nAlphaSubCycles(readLabel(pimpleDict.lookup("nAlphaSubCycles")));

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(interface.cAlpha()*phic, max(phic));

    volScalarField divU(fvc::div(phi));

    dimensionedScalar totalDeltaT = runTime.deltaT();

    if (nAlphaSubCycles > 1)
    {
        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
        }
    }
    else
    {
        #include "alphaEqn.H"
    }

}


rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;

//rho == interface.smoothAlpha()*rho1 + (scalar(1) - interface.smoothAlpha())*rho2;
