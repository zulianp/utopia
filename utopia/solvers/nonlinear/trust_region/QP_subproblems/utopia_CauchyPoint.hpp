#ifndef UTOPIA_TR_SUBPROBLEM_CAUCHY_POINT_HPP
#define UTOPIA_TR_SUBPROBLEM_CAUCHY_POINT_HPP
#include "utopia_TRSubproblem.hpp"

namespace utopia {

    /**
     * @brief
               This class evaluates the Cauchy Point.
            \details
                \f$ p_k^C = - \tau_k \frac{\Delta_k}{||g_k||} g_k \f$, where
                \f$ \tau_k  =         \begin{cases}
                1 &      g_k^T B_k g_k \leq 0   \\
                \min(\frac{ ||g_k||^3 }{\Delta_k g_k^T B_k g_k}, 1 )& otherwise
            \end{cases}  \f$
     */
    template <class Matrix, class Vector>
    class CauchyPoint final : public TRSubproblem<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

    public:
        CauchyPoint() : TRSubproblem<Matrix, Vector>() {}

        bool apply(const Vector &b, Vector &x) override { return aux_solve(*this->get_operator(), b, x); }

        inline CauchyPoint *clone() const override { return new CauchyPoint(); }

        void init_memory(const Layout &layout) override { Bg_.zeros(layout); }

    private:
        bool aux_solve(const Matrix &B, const Vector &g, Vector &p_k) {
            Scalar g_norm, g_B_g;
            Bg_ = B * g;

            dots(g, g, g_norm, g, Bg_, g_B_g);
            g_norm = std::sqrt(g_norm);

            Scalar tau = std::min(1.0, pow(g_norm, 3.0) / (this->current_radius() * g_B_g));
            p_k = tau * (this->current_radius() / g_norm) * (g);
            return true;
        }

        Vector Bg_;
    };
}  // namespace utopia

#endif  // UTOPIA_TR_SUBPROBLEM_CAUCHY_POINT_HPP
