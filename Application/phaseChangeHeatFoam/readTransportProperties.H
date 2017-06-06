 IOdictionary transportPropertiesDict
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word phaseChangeTwoPhaseMixtureTypeName
    (
        transportPropertiesDict.lookup("phaseChangeTwoPhaseMixture")
    );

    bool phaseChange = true;
    if (transportPropertiesDict.found("phaseChange"))
    {
    	phaseChange = Switch(transportPropertiesDict.lookup("phaseChange"));
    }
/*
    dimensionedScalar D
    (
    		transportPropertiesDict.lookup("D")
    );
*/
