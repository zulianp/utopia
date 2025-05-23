#ifndef UTOPIA_EXTENDED_ROSENBROCK_21
#define UTOPIA_EXTENDED_ROSENBROCK_21

#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_Layout.hpp"
#include "utopia_TestFunctions.hpp"

namespace utopia {
    template <class Matrix, class Vector>
    class ExtendedRosenbrock21 final : public UnconstrainedTestFunction<Matrix, Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        ExtendedRosenbrock21(const SizeType &n_loc) : n_loc_(n_loc) { init(Comm::get_default(), n_loc); }

        ExtendedRosenbrock21(const Comm &comm = Comm::get_default(), const SizeType &n_loc = 2) : n_loc_(n_loc) {
            init(comm, n_loc);
        }

        void init(const Comm &comm, const SizeType &n_loc) {
            x_init_.values(layout(comm, n_loc, Traits::determine()), 1.0);

            {
                auto x_view = view_device(x_init_);

                parallel_for(
                    range_device(x_init_),
                    UTOPIA_LAMBDA(const SizeType &i) { x_view.set(i, (i % 2 == 1) ? -1.2 : 1.0); });
            }

            x_exact_.values(layout(x_init_), 1.0);
        }

        std::string name() const override { return "Extended Rosenbrock"; }

        SizeType dim() const override { return n_loc_; }

        bool parallel() const override { return true; }

        bool update(const Vector &point) override {
            init_perm(point);
            x_perm_ = perm_ * point;

            {
                Read<Vector> read_x(x_perm_);
                Range r_x = range(x_perm_);
                xm1 = x_perm_.get(r_x.begin());
                xp1 = x_perm_.get(r_x.begin() + 1);
            }

            // Read<Vector> read_p(point);
            // Range r = r_p(point);
            //
            // // computing value
            // value_ = 0;
            // for (SizeType i = r.begin(); i < r.end(); i++) {
            //     Scalar xi = point.get(i);
            //     Scalar xnext = point.get(i+1);
            //     sum += 100.0 * pow(xnext - xi * xi, 2.0) + pow(xi - 1, 2.0);
            // }

            return true;
        }

        bool value(const Vector &point, Scalar &result) const override {
            const Read<Vector> read(point);
            const Range r = range(point);

            Scalar sum = 0;
            SizeType endm1 = r.end() - 1;
            for (SizeType i = r.begin(); i < endm1; i++) {
                Scalar xi = point.get(i);
                Scalar xnext = point.get(i + 1);
                sum += 100.0 * pow(xnext - xi * xi, 2.0) + pow(xi - 1, 2.0);
            }
            if (r.end() != point.size())
                sum += 100.0 * pow(xp1 - point.get(endm1) * point.get(endm1), 2.0) + pow(point.get(endm1) - 1, 2.0);
            result = sum;
            return true;
        }

        bool gradient(const Vector &point, Vector &result) const override {
            Read<Vector> read(point);

            SizeType d = point.size();
            if (empty(result)) {
                result.zeros(layout(point));
            } else {
                result.set(0.0);
            }

            // for parallel access to result
            Range r_result = range(result);
            if (r_result.empty()) return false;
            SizeType rr_begin = r_result.begin();
            SizeType rr_end = r_result.end();

            // for parallel access to point
            Range r_point = range(point);
            SizeType rp_begin = r_point.begin();
            SizeType rp_end = r_point.end();

            // is this true?
            // assert(rr_begin == rp_begin && rr_end == rp_end);

            // must be put after range otherwise it will crash (segfault)
            Write<Vector> write(result);

            if (rr_begin == 0) {
                Scalar x0 = point.get(0);
                Scalar x1 = (1 == rp_end ? xp1 : point.get(1));
                result.set(0, -400.0 * (x1 - x0 * x0) * x0 + 2.0 * (x0 - 1));
                rr_begin++;
            }

            if (rr_end == d) {
                Scalar dm2 = (d - 1 == rp_begin ? xm1 : point.get(d - 2));
                result.set(d - 1, 200.0 * (point.get(d - 1) - dm2 * dm2));
                rr_end--;
            }

            for (SizeType i = rr_begin; i < rr_end; i++) {
                Scalar xi = point.get(i);
                Scalar xprev = (i == rp_begin ? xm1 : point.get(i - 1));
                Scalar xnext = (i + 1 == rp_end ? xp1 : point.get(i + 1));

                result.set(i, 200.0 * (xi - xprev * xprev) - 400.0 * (xnext - xi * xi) * xi + 2.0 * (xi - 1));
            }

            return true;
        }

        bool hessian(const Vector &point, Matrix &result) const override {
            Read<Vector> read(point);

            SizeType d = point.size();

            auto vl = layout(point);

            if (empty(result)) {
                result.dense(layout(vl.comm(), vl.local_size(), vl.local_size(), vl.size(), vl.size()));
            }

            result.set(0.0);

            // for parallel access to Matrix
            Range r_result = row_range(result);
            if (r_result.empty()) return false;
            SizeType r_begin = r_result.begin();
            SizeType r_end = r_result.end();

            // for parallel access to point
            Range r_point = range(point);
            SizeType rp_begin = r_point.begin();
            SizeType rp_end = r_point.end();

            Write<Matrix> write(result);

            // first row
            if (r_begin == 0 && r_begin != r_end) {
                Scalar x0 = point.get(0);
                Scalar x1 = (1 == rp_end ? xp1 : point.get(1));
                result.set(0, 0, -400.0 * ((x1 - x0 * x0) - 2.0 * x0 * x0) + 2.0);
                result.set(0, 1, -400.0 * x0);
                r_begin++;
            }

            // last row
            if (r_end == d && r_begin != r_end) {
                result.set(d - 1, d - 2, -400.0 * (d - 1 == rp_begin ? xm1 : point.get(d - 2)));
                result.set(d - 1, d - 1, 200.0);
                r_end--;
            }

            // all other rows
            for (SizeType i = r_begin; i < r_end; i++) {
                Scalar xi = point.get(i);
                Scalar xprev = (i == rp_begin ? xm1 : point.get(i - 1));
                Scalar xnext = (i + 1 == rp_end ? xp1 : point.get(i + 1));

                result.set(i, i - 1, -400.0 * xprev);
                result.set(i, i, 200.0 - 400.0 * ((xnext - xi * xi) - 2.0 * xi * xi) + 2.0);
                result.set(i, i + 1, -400.0 * xi);
            }

            return true;
        }

        void clear() {}

        Vector initial_guess() const override { return x_init_; }

        const Vector &exact_sol() const override { return x_exact_; }

    private:
        SizeType n_loc_;
        Matrix perm_;
        Vector x_perm_;
        Scalar xm1, xp1;

        Scalar value_;

        void init_perm(const Vector &x) {
            auto vl = layout(x);

            const SizeType rows_local = 2; 
            const SizeType cols_local =  vl.local_size(); 
            
            const SizeType rows_global = 2 * vl.comm().size(); 
            const SizeType cols_global = vl.size(); 

            auto mat_layout = layout(vl.comm(), rows_local, cols_local, rows_global, cols_global);

            auto r = range(x);
            // long n = local_size(x).get(0);
            long N = size(x).get(0);

            if (is_sparse<Matrix>::value) {
                perm_.sparse(mat_layout, 2, 2);
            } else {
                perm_.dense(mat_layout);
                perm_.set(0.0);
            }

            auto r_perm = row_range(perm_);
            const long i1 = r_perm.begin();
            const long i2 = r_perm.begin() + 1;

            {
                Write<Matrix> w_p(perm_);

                if (r.begin() > 0) {
                    perm_.set(i1, r.begin() - 1, 1);
                }

                if (r.end() != N) {
                    perm_.set(i2, r.end(), 1);
                }
            }
            // disp(perm_);
        }

    private:
        Vector x_init_;
        Vector x_exact_;
    };

}  // namespace utopia
#endif  // UTOPIA_EXTENDED_ROSENBROCK_21