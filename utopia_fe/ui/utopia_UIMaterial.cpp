#include "utopia_UIMaterial.hpp"

#include "utopia_ElasticMaterial.hpp"
#include "utopia_StabilizedMaterial.hpp"
#include "utopia_LameeParameters.hpp"
#include "utopia_LinearElasticity.hpp"
#include "utopia_NewLinearElasticity.hpp"
#include "utopia_NeoHookean.hpp"
#include "utopia_NewNeoHookean.hpp"
#include "utopia_SaintVenantKirchoff.hpp"

#include "utopia_ui.hpp"

namespace utopia {

    template<class FunctionSpace, class Matrix, class Vector>
    void UIMaterial<FunctionSpace, Matrix, Vector>::read(Input &is)
    {
        std::string material = "LinearElasticity";
        std::string stabilization = "none";
        Scalar stabilization_mag = 0.0001;
        Scalar rescaling = 1.0;

        is.get("material", material);
        is.get("stabilization", stabilization);
        is.get("stabilization-mag", stabilization_mag);
        is.get("parameters", params);
        is.get("rescaling", rescaling);

        params.describe(std::cout);

        if(material == "NeoHookean") {
            std::cout << "Using: NeoHookean" << std::endl;
            material_ = std::make_shared<NeoHookean<FunctionSpace, Matrix, Vector>>(V_, params);
        } else if(material == "NewNeoHookean") {
            std::cout << "Using: NewNeoHookean" << std::endl;
            material_ = std::make_shared<NewNeoHookean<FunctionSpace, Matrix, Vector>>(V_, params);
        } else if(material == "NewLinearElasticity") {
            std::cout << "Using: NewLinearElasticity" << std::endl;
            material_ = std::make_shared<NewLinearElasticity<FunctionSpace, Matrix, Vector>>(V_, params);
        } else if(material == "SaintVenantKirchoff") {
            std::cout << "Using: SaintVenantKirchoff" << std::endl;
            material_ = std::make_shared<SaintVenantKirchoff<FunctionSpace, Matrix, Vector>>(V_, params);
        } else /*if(material == "LinearElasticity")*/ {
        	std::cout << "Using: LinearElasticity" << std::endl;
            material_ = std::make_shared<LinearElasticity<FunctionSpace, Matrix, Vector>>(V_, params);
        }

        if(stabilization != "none") {
            std::cout << "using stabilization: " << stabilization << " mag: " << stabilization_mag << std::endl;
            // StabilizedMaterial<FunctionSpace, Matrix, Vector> sm(V_, stabilization_mag, material_, stabilization);
            material_ = std::make_shared<StabilizedMaterial<FunctionSpace, Matrix, Vector>>(V_, stabilization_mag, material_, stabilization);
        }

        material_->rescaling(rescaling);
    }

    template class UIMaterial<ProductFunctionSpace<LibMeshFunctionSpace>, USparseMatrix, UVector>;
}

