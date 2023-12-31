#ifndef UTOPIA_UFLOW_HPP
#define UTOPIA_UFLOW_HPP

#include "utopia_Model.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"

namespace utopia {
    template <class FunctionSpace, class Matrix, class Vector>
    class UFlow final : public Model<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        typedef utopia::Traits<FunctionSpace> TraitsT;
        typedef typename TraitsT::Matrix ElementMatrix;

        UFlow(FunctionSpace &space);

        void read(Input &in) override;
        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override;
        bool assemble_hessian(Matrix &hessian) override;
        bool assemble_hessian(const Vector &x, Matrix &hessian) override;
        bool assemble_gradient(const Vector &x, Vector &gradient) override;

        inline bool is_linear() const override {
            assert(is_valid());
            return flow_model_->is_linear();
        }

        inline void clear() override {
            assert(is_valid());
            flow_model_->clear();
        }

        inline void rescale(const Scalar rescale) {
            assert(!is_valid() && "rescale has to be called before read(Input &in)");
            rescale_ = rescale;
        }

        inline bool is_valid() const { return static_cast<bool>(flow_model_); }

        inline Scalar rescaling() const { return rescale_; }

    private:
        FunctionSpace &space_;
        std::unique_ptr<Model<Matrix, Vector>> flow_model_;
        Scalar rescale_{1.0};
    };
}  // namespace utopia

#endif  // UTOPIA_UFLOW_HPP
