#ifndef UTOPIA_ASSEMBLE_LAPLACIAN_1D
#define UTOPIA_ASSEMBLE_LAPLACIAN_1D

#include "utopia_Base.hpp"
#include "utopia_Range.hpp"
#include "utopia_Writable.hpp"

namespace utopia {

    template <class Matrix>
    void assemble_laplacian_1D(Matrix &m, const bool bc = false) {
        // n x n matrix with maximum 3 entries x row
        Write<Matrix> w(m);
        Range r = row_range(m);
        auto n = size(m).get(0);

        for (SizeType i = r.begin(); i != r.end(); ++i) {
            if (bc) {
                if (i == 0) {
                    m.set(0, 0, 1.);
                    continue;
                }

                if (i == n - 1) {
                    m.set(n - 1, n - 1, 1.);
                    continue;
                }
            }

            if (i > 0) {
                m.set(i, i - 1, -1.0);
            }

            if (i < n - 1) {
                m.set(i, i + 1, -1.0);
            }

            if (i == 0 || i == n - 1) {
                m.set(i, i, 1.);
            } else {
                m.set(i, i, 2.0);
            }
        }
    }

    template <class Matrix>
    void assemble_symmetric_laplacian_1D(Matrix &m, const bool bc = false) {
        // n x n matrix with maximum 3 entries x row
        Write<Matrix> w(m);
        Range r = row_range(m);
        auto n = size(m).get(0);

        for (SizeType i = r.begin(); i != r.end(); ++i) {
            if (i > 0) {
                m.set(i, i - 1, -1.0);
            }

            if (i < n - 1) {
                m.set(i, i + 1, -1.0);
            }

            if (i == 0 || i == n - 1) {
                m.set(i, i, 1.);
            } else {
                m.set(i, i, 2.0);
            }
        }

        if (bc) {
            if (r.inside(0)) {
                m.set(0, 0, 1.);
                m.set(0, 1, 0.);
            }

            if (r.inside(1)) {
                m.set(1, 0, 0.);
            }

            if (r.inside(n - 1)) {
                m.set(n - 1, n - 1, 1.);
                m.set(n - 1, n - 2, 0.);
            }

            if (r.inside(n - 2)) {
                m.set(n - 2, n - 1, 0.);
            }
        }
    }

    template <class Matrix>
    void assemble_laplacian_1D_with_scaling(Matrix &m, const double scale_factor, const bool bc = false) {
        // n x n matrix with maximum 3 entries x row
        Write<Matrix> w(m);
        Range r = row_range(m);
        auto n = size(m).get(0);

        for (SizeType i = r.begin(); i != r.end(); ++i) {
            if (bc && (i == 0 || i == n - 1)) {
                m.set(i, i, 1.);
                continue;
            } else {
                m.set(i, i, 2.0 * scale_factor);
            }

            if (i > 0) {
                m.set(i, i - 1, -1.0 * scale_factor);
            }

            if (i < n - 1) {
                m.set(i, i + 1, -1.0 * scale_factor);
            }
        }
    }

    template <class Matrix, class Vector>
    void assemble_poisson_problem_1D(const typename Traits<Vector>::Scalar &f,
                                     Matrix &A,
                                     Vector &b,
                                     const bool use_mesh_size = true) {
        auto n = A.rows();
        auto h = 1. / (n - 1);
        assemble_laplacian_1D_with_scaling(A, use_mesh_size ? (1. / h) : 1.0, true);

        b.set(f * (use_mesh_size ? h : 1.0));

        Write<Vector> w(b);

        auto r = b.range();

        if (r.inside(0)) {
            b.set(0, 0.0);
        }

        if (r.inside(n - 1)) {
            b.set(n - 1, 0.0);
        }
    }

}  // namespace utopia

#endif  // UTOPIA_ASSEMBLE_LAPLACIAN_1D
