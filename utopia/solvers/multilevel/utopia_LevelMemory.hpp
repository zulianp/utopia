#ifndef UTOPIA_LEVEL_MEMORY_HPP
#define UTOPIA_LEVEL_MEMORY_HPP

#include <limits>
#include <vector>
#include "utopia_Traits.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class FASLevelMemory {
    public:
        void init(const int n_levels) {
            x.resize(n_levels);
            x_0.resize(n_levels);

            g.resize(n_levels);
            g_diff.resize(n_levels);

            c.resize(n_levels);
        }

        std::vector<Vector> x, x_0, g, g_diff, c;
    };

    template <class Matrix, class Vector>
    class JFNKLevelMemory {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

    public:
        void init_memory(const std::vector<Layout> &layouts) {
            const SizeType n_levels = layouts.size();
            x.resize(n_levels);
            g.resize(n_levels);

            res.resize(n_levels);
            rhs.resize(n_levels);
            c.resize(n_levels);

            for (SizeType l = 0; l < n_levels; l++) {
                x[l].zeros(layouts[l]);
                g[l].zeros(layouts[l]);

                res[l].zeros(layouts[l]);
                rhs[l].zeros(layouts[l]);
                c[l].zeros(layouts[l]);
            }
        }

        std::vector<Vector> x, rhs, g, res, c;
    };

    template <class Matrix, class Vector>
    class RMTRLevelMemory {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

    public:
        void init_memory(const std::vector<Layout> &layouts) {
            const SizeType n_levels = layouts.size();

            x.resize(n_levels);
            x_0.resize(n_levels);

            s.resize(n_levels);
            s_working.resize(n_levels);
            help.resize(n_levels);

            delta.resize(n_levels);
            energy.resize(n_levels);
            gnorm.resize(n_levels);

            for (SizeType l = 0; l < n_levels; l++) {
                x[l].zeros(layouts[l]);
                x_0[l].zeros(layouts[l]);
                s[l].zeros(layouts[l]);
                s_working[l].zeros(layouts[l]);
                help[l].zeros(layouts[l]);
            }
        }

        std::vector<Scalar> delta, energy, gnorm;
        std::vector<Vector> x, x_0, s, s_working, help;
    };

    template <class Vector>
    class ConstraintsLevelMemory {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

    public:
        void init_memory(const std::vector<Layout> &layouts) {
            const Scalar inf = std::numeric_limits<Scalar>::infinity();
            const SizeType n_levels = layouts.size();
            active_lower.resize(n_levels);
            active_upper.resize(n_levels);

            for (SizeType l = 0; l < n_levels; l++) {
                active_lower[l].values(layouts[l], -inf);
                active_upper[l].values(layouts[l], inf);
            }
        }

        std::vector<Vector> active_lower, active_upper;
    };

}  // namespace utopia

#endif  // UTOPIA_LEVEL_MEMORY_HPP
