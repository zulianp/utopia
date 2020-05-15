#ifndef UTOPIA_POLYGON_HPP
#define UTOPIA_POLYGON_HPP

namespace utopia {
    template <int Dimensions, typename T = double>
    class Polygon {
    public:
        Polygon(const int n_vertices) : coordinates_(n_vertices * Dimensions) {}

        inline bool is_convex() const { return false; }

        inline T *point(const int index) { return &coordinates_[index * Dimensions]; }

        inline const T *point(const int index) const { return &coordinates_[index * Dimensions]; }

    private:
        std::vector<T> coordinates_;

        inline std::vector<T> &coordinates() { return coordinates_; }

        inline const std::vector<T> &coordinates() const { return coordinates_; }
    };

    // returns true if the result is usefull
    template <int Dimensions, typename T, class Fun>
    bool split_curve(Fun fun,
                     const T *from,
                     const T *to,
                     const int polynomial_order,
                     std::vector<T> &result,
                     const T tol) {
        T f_from[Dimensions];
        T f_to[Dimensions];
        T f_q[Dimensions];
        T f_q2[Dimensions];
        const T *f_last;

        fun(from, f_from);
        fun(to, f_to);

        const int n_new_points = polynomial_order - 1;
        result.resize(n_new_points * Dimensions);

        T dist = 0;
        T v[Dimensions];
        for (int d = 0; d < Dimensions; ++d) {
            v[d] = to[d] - from[d];
            const T dd = f_to[d] - f_from[d];
            dist += dd * dd;
        }

        dist = std::sqrt(dist);

        const T h = 1. / polynomial_order;

        f_last = f_from;
        T cum_dist = 0.0;

        for (int i = 0; i < n_new_points; ++i) {
            T *q = &result[i * Dimensions];

            for (int d = 0; d < Dimensions; ++d) {
                q[d] = from[d] + ((i + 1) * h) * v[d];
            }

            fun(q, f_q2);

            T ddd = 0.0;
            for (int d = 0; d < Dimensions; ++d) {
                const T dd = f_q2[d] - f_last[d];
                ddd += dd * dd;

                //! overwriting f_last[d] (which is not need anymore but )
                f_q[d] = f_q2[d];
            }

            cum_dist += std::sqrt(ddd);
            f_last = f_q;
        }

        T ddd = 0.0;
        for (int d = 0; d < Dimensions; ++d) {
            const T dd = f_to[d] - f_last[d];
            ddd += dd * dd;
        }

        assert(ddd > 0);

        cum_dist += std::sqrt(ddd);
        const T diff = std::abs(cum_dist - dist);
        return diff > tol;
    }

    template <int Dimensions, typename T, class Fun>
    void discretize_curve(Fun fun,
                          const T *from,
                          const T *to,
                          const int polynomial_order,
                          std::vector<T> &param_points,
                          std::vector<T> &result,
                          const T tol = 1e-4,
                          const int recursion = 0) {
        // yes I know max harcoded 20 recursions not that nice
        if (recursion == 20) {
            return;
        }

        if (recursion == 0) {
            param_points.clear();
            param_points.insert(param_points.end(), from, from + Dimensions);

            T f_from[Dimensions];
            fun(from, f_from);
            result.clear();
            result.insert(result.end(), f_from, f_from + Dimensions);
        }

        std::vector<T> split_res;
        if (split_curve<Dimensions>(fun, from, to, polynomial_order, split_res, tol)) {
            split_res.insert(split_res.begin(), from, from + Dimensions);
            split_res.insert(split_res.end(), to, to + Dimensions);

            for (int i = 0; i < split_res.size() - Dimensions; i += Dimensions) {
                if (i > 0) {
                    param_points.insert(param_points.end(), &split_res[i], &split_res[i] + Dimensions);

                    T f_p[Dimensions];
                    fun(&split_res[i], f_p);
                    result.insert(result.end(), f_p, f_p + Dimensions);
                }

                discretize_curve<Dimensions>(fun,
                                             &split_res[i],
                                             &split_res[i + Dimensions],
                                             polynomial_order,
                                             param_points,
                                             result,
                                             tol,
                                             recursion + 1);
            }
        }

        if (recursion == 0) {
            param_points.insert(param_points.end(), to, to + Dimensions);

            T f_to[Dimensions];
            fun(to, f_to);
            result.insert(result.end(), f_to, f_to + Dimensions);
        }
    }

}  // namespace utopia

#endif  // UTOPIA_POLYGON_HPP
