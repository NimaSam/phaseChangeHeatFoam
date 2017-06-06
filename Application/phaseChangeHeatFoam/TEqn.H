Info << "solve TEqn" << endl;
{

	Pair<tmp<volScalarField> > vDotT = twoPhaseProperties->vDotT();
		const volScalarField& vDotcT = vDotT[0]();
		const volScalarField& vDotvT = vDotT[1]();

	volScalarField k= K1*alpha1 + (1.0 - alpha1)*K2;
	volScalarField rhoC = rho1*C1*alpha1 + (1.0 - alpha1)*rho2*C2;


	K = 0.5*magSqr(U);

	fvScalarMatrix TEqn
	(
		  fvm::ddt(T)
		+ fvm::div(phi, T)
	    - fvm::Sp(fvc::div(phi),T)
		- fvm::laplacian( k/rhoC , T,"laplacian(alphaEff,T)")
	/*	+ (
			fvc::div(fvc::absolute(phi, U), p)
		  + fvc::ddt(rho, K) + fvc::div(rhoPhi, K)
		  )/rhoC
    */
	);

	TEqn.relax();

    if (phaseChange){
		solve
		 (
			TEqn
			==
			//fvm::Sp(vDotvT-vDotcT,T) - (vDotvT-vDotcT)*TSat
			fvm::Sp(vDotcT-vDotvT,T) - (vDotcT-vDotvT)*TSatLocal
		 );
	 }
	 else
	 {
		 TEqn.solve();
	 }


	Info << "TAve = "
	<< T.weightedAverage(mesh.V()).value()
	<< "  Min(T) = " << min(T).value()
	<< "  Max(T) = " << max(T).value()
	<< endl;

	Pair<tmp<volScalarField> > mDot = twoPhaseProperties->mDotT();
	const volScalarField& mDotcT = mDot[0]();
	const volScalarField& mDotvT = mDot[1]();

	mDotT = (mDotcT-mDotvT)*(T-TSatLocal);
}

