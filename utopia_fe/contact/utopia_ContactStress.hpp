#ifndef UTOPIA_CONTACT_STRESS_HPP
#define UTOPIA_CONTACT_STRESS_HPP

#include "utopia_LameeParameters.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_QuadratureUtils.hpp"
#include "utopia_LinearElasticity.hpp"

#include "utopia_TransferUtils.hpp"

namespace utopia {

    template<class FunctionSpaceT, class Matrix, class Vector>
    class ContactStress final {
    public:
        typedef utopia::Traits<FunctionSpaceT> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;
        typedef typename TraitsT::Vector ElementVector;
        using ElasticMaterialT  = utopia::ElasticMaterial<Matrix, Vector>;
        using LinearElasticityT = utopia::LinearElasticity<FunctionSpaceT, Matrix, Vector>;

        ContactStress(FunctionSpaceT &V, const LameeParameters &params)
        : V_(V), elast_(utopia::make_unique<LinearElasticityT>(P1_, params))
        {
            init();
        }

        template<class ElasticityT>
        ContactStress(FunctionSpaceT &V, const LameeParameters &params)
        : V_(V), elast_(utopia::make_unique<ElasticityT>(P1_, params))
        {
            init();
        }

        void init();
        bool assemble(const UVector &x, UVector &result);

    private:
        FunctionSpaceT &V_;
        FunctionSpaceT P1_;
        AssemblyContext<LIBMESH_TAG> ctx_;

        USparseMatrix P1toV_, VtoP1_;

        std::unique_ptr<ElasticMaterialT> elast_;
        UVector x_p1_, stress_p1_;

        UVector inverse_mass_vector_;
    };

}

#endif //UTOPIA_CONTACT_STRESS_HPP
