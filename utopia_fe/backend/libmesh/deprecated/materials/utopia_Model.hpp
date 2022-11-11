#ifndef UTOPIA_MODEL_HPP
#define UTOPIA_MODEL_HPP

#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_Types.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class Model : public Configurable {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        virtual ~Model() {}
        virtual bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) = 0;
        virtual bool is_linear() const { return false; }

        virtual bool assemble_hessian(const Vector &, Matrix &) = 0;
        // {
        //     assert(false);
        //     return false;
        // }

        virtual bool assemble_gradient(const Vector &, Vector &) = 0;
        // {
        //     assert(false);
        //     return false;
        // }

        virtual bool assemble_hessian(Matrix &) {
            assert(false);
            return false;
        }

        virtual void clear() {}
        void read(Input &in) override { UTOPIA_UNUSED(in); }

        virtual void initialize() {}
    };

}  // namespace utopia

#endif  // UTOPIA_MODEL_HPP
