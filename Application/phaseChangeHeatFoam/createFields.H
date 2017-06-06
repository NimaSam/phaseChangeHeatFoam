    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField p
       (
       	IOobject
       	(
       		"p",
       		runTime.timeName(),
       		mesh,
       		IOobject::NO_READ,
       		IOobject::AUTO_WRITE
       	),
       	p_rgh
       );

    #include "createPhi.H"

    Info<< "Creating phaseChangeTwoPhaseMixture\n" << endl;
    autoPtr<phaseChangeTwoPhaseMixture> twoPhaseProperties =
        phaseChangeTwoPhaseMixture::New(U, phi);

    	const dimensionedScalar& rho1 = twoPhaseProperties->rho1();
        const dimensionedScalar& rho2 = twoPhaseProperties->rho2();

        const dimensionedScalar& K1 = twoPhaseProperties->K1();
        const dimensionedScalar& K2 = twoPhaseProperties->K2();

        const dimensionedScalar& C1 = twoPhaseProperties->C1();
        const dimensionedScalar& C2 = twoPhaseProperties->C2();

        const volScalarField& TSatLocal = twoPhaseProperties->TSatLocal();
        const dimensionedScalar& pSat = twoPhaseProperties->pSat();


      // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();

    // Construct interface from alpha1 distribution
    //interfaceProperties interface(alpha1, U, twoPhaseProperties());
    smoothInterfaceProperties interface(alpha1, U, twoPhaseProperties());

    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties())
    );


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"
    p = p + rho*gh;

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

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


    fv::IOoptionList fvOptions(mesh); //add

    surfaceScalarField rhoPhi
    (
    		IOobject
    		(
    				"rhoPhi",
    				runTime.timeName(),
    				mesh
    		),
    		mesh,
    		dimensionedScalar("0", dimMass/dimTime, 0)
    );

#include "readTransportProperties.H"

    Info<< "Creating field kinetic energy K\n" << endl;
        volScalarField K("K", 0.5*magSqr(U));

        volScalarField mDotT
           (
           		IOobject
                   (
                       "mDotT",
                       runTime.timeName(),
                       mesh,
                       IOobject::NO_READ,
                       IOobject::AUTO_WRITE
                   ),
                mesh,
                dimensionedScalar ("mDotT",dimDensity/dimTime,0)
           );

        volScalarField mDotP
                  (
                  		IOobject
                          (
                              "mDotP",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ,
                              IOobject::AUTO_WRITE
                          ),
                       mesh,
                       dimensionedScalar ("mDotP",dimDensity/dimTime,0)
                  );
        volScalarField mDotAlpha
                  (
                  		IOobject
                          (
                              "mDotAlpha",
                              runTime.timeName(),
                              mesh,
                              IOobject::NO_READ,
                              IOobject::AUTO_WRITE
                          ),
                       mesh,
                       dimensionedScalar ("mDotAlpha",dimDensity/dimTime,0)
                  );
