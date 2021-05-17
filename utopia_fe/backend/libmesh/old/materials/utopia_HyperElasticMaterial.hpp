#ifndef UTOPIA_HYPER_ELASTIC_MATERIAL_HPP
#define UTOPIA_HYPER_ELASTIC_MATERIAL_HPP

#include "utopia_ElasticMaterial.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class HyperElasticMaterial : public ElasticMaterial<Matrix, Vector> {
    public:
        virtual ~HyperElasticMaterial() {}
    };
}  // namespace utopia

#endif  // UTOPIA_HYPER_ELASTIC_MATERIAL_HPP
