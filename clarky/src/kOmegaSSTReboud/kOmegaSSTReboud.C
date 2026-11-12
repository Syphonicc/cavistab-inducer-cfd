#include "kOmegaSSTReboud.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
void kOmegaSSTReboud<BasicTurbulenceModel>::applyReboud()
{
    const auto* alphaPtr =
        this->mesh_.template findObject<volScalarField>(alphaName_);

    if (!alphaPtr)
    {
        // Phase fraction not yet registered (e.g. during construction)
        return;
    }

    const volScalarField alphaL(max(min(*alphaPtr, scalar(1)), scalar(0)));
    const volScalarField rhoM(alphaL*rhoL_ + (1 - alphaL)*rhoV_);
    const volScalarField fRho(rhoV_ + pow(alphaL, nReboud_)*(rhoL_ - rhoV_));

    this->nut_.primitiveFieldRef() *= (fRho/rhoM)().primitiveField();
    this->nut_.correctBoundaryConditions();
}


template<class BasicTurbulenceModel>
void kOmegaSSTReboud<BasicTurbulenceModel>::correctNut(const volScalarField& S2)
{
    kOmegaSST<BasicTurbulenceModel>::correctNut(S2);
    applyReboud();
    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void kOmegaSSTReboud<BasicTurbulenceModel>::correctNut()
{
    correctNut(2*magSqr(symm(fvc::grad(this->U_))));
}


template<class BasicTurbulenceModel>
kOmegaSSTReboud<BasicTurbulenceModel>::kOmegaSSTReboud
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSST<BasicTurbulenceModel>
    (
        alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName, type
    ),
    alphaName_
    (
        this->coeffDict_.template getOrDefault<word>("alphaName", "alpha.water")
    ),
    rhoL_
    (
        dimensionedScalar::getOrAddToDict
        (
            "rhoL", this->coeffDict_, dimDensity, 999.19
        )
    ),
    rhoV_
    (
        dimensionedScalar::getOrAddToDict
        (
            "rhoV", this->coeffDict_, dimDensity, 0.02308
        )
    ),
    nReboud_(this->coeffDict_.template getOrDefault<scalar>("nReboud", 10))
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


template<class BasicTurbulenceModel>
bool kOmegaSSTReboud<BasicTurbulenceModel>::read()
{
    if (kOmegaSST<BasicTurbulenceModel>::read())
    {
        this->coeffDict_.readIfPresent("alphaName", alphaName_);
        rhoL_.readIfPresent(this->coeffDict_);
        rhoV_.readIfPresent(this->coeffDict_);
        this->coeffDict_.readIfPresent("nReboud", nReboud_);
        return true;
    }
    return false;
}

} // End namespace RASModels
} // End namespace Foam
