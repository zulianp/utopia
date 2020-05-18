#ifndef UTOPIA_WEAK_DIRICHLET_BOUNDARY_CONDITIONS_HPP
#define UTOPIA_WEAK_DIRICHLET_BOUNDARY_CONDITIONS_HPP

#include "utopia_ForcedMaterial.hpp"
#include "utopia_ui.hpp"

namespace utopia {

    template <class FunctionSpace, class Matrix, class Vector>
    class WeakDirichletBoundaryConditions final : public Configurable {
    public:
        WeakDirichletBoundaryConditions(FunctionSpace &V) : V_(V) {}

        ~WeakDirichletBoundaryConditions() {}

        void read(Input &is) override {
            is.get_all([this](Input &is) {
                int side = -1;
                double penalty = 1e8;
                std::string function_type = "expr";

                is.get("side", side);
                is.get("penalty", penalty);

                assert(side != -1);

                if (side == -1) {
                    std::cerr << "[Error] side not specified" << std::endl;
                    return;
                }

#ifdef WITH_TINY_EXPR
                std::string value;
                is.get("value", value);
                auto f = symbolic(value);
#else
                double value = 0.;
                is.get("value", value);
                auto f = coeff(value);
#endif  // WITH_TINY_EXPR

                std::cout << "[Status] weak bc side = " << side << ", value = " << value << ", penalty = " << penalty
                          << std::endl;

                auto u = trial(V_);
                auto v = test(V_);

                auto b_form = surface_integral(penalty * inner(u, v), side);
                auto l_form = surface_integral(penalty * inner(f, v), side);

                Matrix B_temp;
                Vector Bf_temp;
                utopia::assemble(b_form == l_form, B_temp, Bf_temp);

                if (empty(B_)) {
                    B_ = std::move(B_temp);
                } else {
                    B_ += B_temp;
                }

                if (empty(Bf_)) {
                    Bf_ = std::move(Bf_temp);
                } else {
                    Bf_ += Bf_temp;
                }
            });
        }

        inline void apply(Matrix &A, Vector &b) const {
            if (empty(B_)) return;

            A += B_;
            b += Bf_;
        }

    private:
        FunctionSpace &V_;

        Matrix B_;
        Vector Bf_;
    };
}  // namespace utopia

#endif  // UTOPIA_WEAK_DIRICHLET_BOUNDARY_CONDITIONS_HPP
