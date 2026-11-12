/*---------------------------------------------------------------------------*\
Class
    Foam::RASModels::kOmegaSSTReboud

Description
    k-omega SST with the Reboud et al. (1998) density correction of the eddy
    viscosity for homogeneous-mixture cavitating flow:

        mu_t  = f(rho) Cmu k/omega-type SST expression,
        f(rho) = rho_v + alpha_l^n (rho_l - rho_v),

    applied to the kinematic eddy viscosity as nut *= f(rho)/rho_m,
    with rho_m = alpha_l rho_l + (1 - alpha_l) rho_v. Recovers standard SST in
    pure liquid and pure vapour; strongly reduces nut in the mixture region.

    Coefficients (kOmegaSSTReboudCoeffs / RAS dict):
        alphaName   alpha.water;   // liquid volume fraction field
        rhoL        999.19;
        rhoV        0.02308;
        nReboud     10;
\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTReboud_H
#define kOmegaSSTReboud_H

#include "kOmegaSST.H"

namespace Foam
{
namespace RASModels
{

template<class BasicTurbulenceModel>
class kOmegaSSTReboud
:
    public kOmegaSST<BasicTurbulenceModel>
{
    // Private Data

        word alphaName_;
        dimensionedScalar rhoL_;
        dimensionedScalar rhoV_;
        scalar nReboud_;

        kOmegaSSTReboud(const kOmegaSSTReboud&) = delete;
        void operator=(const kOmegaSSTReboud&) = delete;


protected:

        //- Apply the density correction to nut
        void applyReboud();

        virtual void correctNut(const volScalarField& S2);
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    TypeName("kOmegaSSTReboud");

    kOmegaSSTReboud
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName
    );

    virtual ~kOmegaSSTReboud() = default;

    virtual bool read();
};

} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "kOmegaSSTReboud.C"
#endif

#endif
