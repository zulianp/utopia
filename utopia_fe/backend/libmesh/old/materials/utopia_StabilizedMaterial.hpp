#ifndef UTOPIA_STABILIZED_MATERIAL_HPP
#define UTOPIA_STABILIZED_MATERIAL_HPP

#include "utopia.hpp"
#include "utopia_ElasticMaterial.hpp"
#include "utopia_fe_EDSL.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"

namespace utopia {

    template <class FunctionSpace, class Matrix, class Vector>
    class StabilizedMaterial final : public ElasticMaterial<Matrix, Vector> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        enum StabilizationType { L2 = 0, H1 = 1, L2_LUMPED = 2 };

        StabilizedMaterial(FunctionSpace &V,
                           const double stabilization_mag,
                           const std::shared_ptr<ElasticMaterial<Matrix, Vector>> &material,
                           const StabilizationType type = H1)
            : V_(V), stabilization_mag_(stabilization_mag), material_(material), type_(type) {}

        StabilizedMaterial(FunctionSpace &V,
                           const double stabilization_mag,
                           const std::shared_ptr<ElasticMaterial<Matrix, Vector>> &material,
                           const std::string type)
            : V_(V), stabilization_mag_(stabilization_mag), material_(material) {
            if (type == "L2" || type == "l2") {
                type_ = L2;
            } else if (type == "L2_LUMPED" || type == "l2_lumped") {
                type_ = L2_LUMPED;
            } else {
                type_ = H1;
            }
        }

        bool assemble_hessian(Matrix &hessian) override {
            if (!material_->assemble_hessian(hessian_)) {
                std::cerr << "decorated material failed to assemble" << std::endl;
                return false;
            }

            if (empty(stab_)) {
                auto u = trial(V_);
                auto v = test(V_);

                switch (type_) {
                    case L2: {
                        utopia::assemble(inner(u, v) * dX, stab_);
                        break;
                    }
                    case L2_LUMPED: {
                        Matrix temp;
                        utopia::assemble(inner(u, v) * dX, temp);
                        Vector d = sum(temp, 1);
                        stab_ = diag(d);
                        break;
                    }
                    default: {
                        utopia::assemble(inner(grad(u), grad(v)) * dX, stab_);
                        break;
                    }
                }

                stab_ *= stabilization_mag_;
            }

            hessian = hessian_ + stab_;
            return true;
        }

        bool assemble_hessian(const Vector &x, Matrix &hessian) override {
            if (!material_->assemble_hessian(x, hessian_)) {
                std::cerr << "decorated material failed to assemble" << std::endl;
                return false;
            }

            if (empty(stab_)) {
                auto u = trial(V_);
                auto v = test(V_);

                switch (type_) {
                    case L2: {
                        utopia::assemble(inner(u, v) * dX, stab_);
                        break;
                    }
                    case L2_LUMPED: {
                        Matrix temp;
                        utopia::assemble(inner(u, v) * dX, temp);
                        Vector d = sum(temp, 1);
                        stab_ = diag(d);
                        break;
                    }
                    default: {
                        utopia::assemble(inner(grad(u), grad(v)) * dX, stab_);
                        break;
                    }
                }

                stab_ *= stabilization_mag_;
            }

            hessian = hessian_ + stab_;
            return true;
        }

        bool assemble_gradient(const Vector &x, Vector &gradient) override {
            return material_->assemble_gradient(x, gradient);
        }

        bool assemble_hessian_and_gradient(const Vector &x, Matrix &hessian, Vector &gradient) override {
            if (!material_->assemble_hessian_and_gradient(x, hessian_, gradient)) {
                std::cerr << "decorated material failed to assemble" << std::endl;
                return false;
            }

            if (empty(stab_)) {
                auto u = trial(V_);
                auto v = test(V_);

                switch (type_) {
                    case L2: {
                        utopia::assemble(inner(u, v) * dX, stab_);
                        break;
                    }
                    case L2_LUMPED: {
                        Matrix temp;
                        utopia::assemble(inner(u, v) * dX, temp);
                        Vector d = sum(temp, 1);
                        stab_ = diag(d);
                        break;
                    }
                    default: {
                        utopia::assemble(inner(grad(u), grad(v)) * dX, stab_);
                        break;
                    }
                }

                stab_ *= stabilization_mag_;
            }

            hessian = hessian_ + stab_;
            return true;
        }

        bool stress(const Vector &x, Vector &result) override { return material_->stress(x, result); }

        inline Scalar rescaling() const override { return material_->rescaling(); }

        inline void rescaling(const Scalar &value) override { material_->rescaling(value); }

    private:
        FunctionSpace &V_;
        double stabilization_mag_;
        std::shared_ptr<ElasticMaterial<Matrix, Vector>> material_;

        Matrix hessian_;
        Matrix stab_;
        StabilizationType type_;
    };
}  // namespace utopia

#endif  // UTOPIA_NEOOHOOKEAN_HPP
