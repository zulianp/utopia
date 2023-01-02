#ifndef UTOPIA_ELASTIC_MATERIAL_HPP
#define UTOPIA_ELASTIC_MATERIAL_HPP

#include "utopia_Model.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_Types.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class ElasticMaterial : public Model<Matrix, Vector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        virtual ~ElasticMaterial() {}
        // virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) = 0;
        // virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) = 0;

        virtual bool stress(const Vector &x, Vector &result) {
            assert(false && "implement me");
            return false;
        }

        virtual bool normal_stress(const UVector &x, UVector &out, const int subspace = 0) {
            // out = local_values(local_size(x).get(0), -1);
            out.values(layout(x), -1);
            return false;
        }

        virtual bool von_mises_stress(const UVector &x, UVector &out, const int subspace = 0) {
            // out = local_values(local_size(x).get(0), -1);
            out.values(layout(x), -1);
            return false;
        }

        // virtual void clear() {}
        // virtual bool is_linear() const { return false; }

        virtual Scalar rescaling() const { return 1.0; }

        virtual void rescaling(const Scalar &) {}
    };

}  // namespace utopia

#endif  // UTOPIA_ELASTIC_MATERIAL_HPP
