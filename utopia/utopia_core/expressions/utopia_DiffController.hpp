#ifndef UTOPIA_UTOPIA_DIFFCONTROLLER_HPP
#define UTOPIA_UTOPIA_DIFFCONTROLLER_HPP

#include "utopia_FiniteDifference.hpp"
#include "utopia_IOStream.hpp"
#include "utopia_Input.hpp"
#include "utopia_Rename.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Wrapper.hpp"

namespace utopia {

    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class DiffController : public Configurable {
    public:
        using Scalar = typename Traits<Vector>::Scalar;

        DiffController(const Scalar spacing = 1e-5) : spacing_(spacing) {}

        inline void hessian_from_grad(const bool val) { hessian_from_grad_ = val; }

        void read(Input &in) override {
            in.get("spacing", spacing_);
            in.get("hessian_from_grad", hessian_from_grad_);
        }

        template <class Fun>
        bool check(const Fun &fun, const Vector &x, const Vector &g, const Matrix &H) const {
            return check_grad(fun, x, g) && check_hessian(fun, x, H);
        }

        template <class Fun>
        bool check_grad(const Fun &fun, const Vector &x, const Vector &g) const {
            FiniteDifference<typename Vector::Scalar> fd(spacing_);
            Vector gfd;
            fd.grad(fun, x, gfd);
            bool ok = approxeq(gfd, g, 1e-7);

            Scalar diff = norm1(gfd - g);

            if (!ok) {
                utopia::err() << "------ Failure -------\n";
                utopia::err() << "------ Gradient ------\n";
                std::cerr << "error: " << diff << std::endl;

                utopia::err() << "Expected:\n";
                rename("g_fd", gfd);
                disp(gfd);
                write("G_fd.m", gfd);

                utopia::err() << "Actual:\n";
                 disp(g);
                rename("g", const_cast<Vector &>(g));
                write("G.m", g);

                utopia::err() << "----------------------\n";
                assert(false);
            } else {
                utopia::err() << "------ Gradient Check Passed -------\n";
                std::cerr << "error: " << diff << std::endl;
            }

            return ok;
        }

        template <class Fun>
        bool check_hessian(const Fun &fun, const Vector &x, const Matrix &H) const {
            FiniteDifference<Scalar> fd(spacing_);
            Matrix Hfd = H;
            Hfd *= 0.0;

            if (!fd.hessian(fun, x, Hfd, hessian_from_grad_)) {
                return true;
            }

            Matrix diff_mat = H - Hfd;
            Scalar diff = norm1(diff_mat);

            bool ok = diff < 1e-7;

            if (!ok) {
                std::cerr << "------- Failure -------\n";
                std::cerr << "------- Hessian -------\n";

                std::cerr << "error: " << diff << std::endl;

                std::cerr << "Diff:\n";
                disp(diff_mat);

                rename("D", diff_mat);
                write("diff_mat.m", diff_mat);

                rename("H_fd", Hfd);
                write("Hfd.m", Hfd);

                rename("h", const_cast<Matrix &>(H));
                write("H_utopia.m", H);

                 std::cerr << "-------------------------\n";
                 std::cerr << "Expected:\n";
                 std::cerr << "-------------------------\n";
                 disp(Hfd);

                 std::cerr << "-------------------------\n";
                 std::cerr << "Actual:\n";
                 std::cerr << "-------------------------\n";
                 disp(H);

                std::cerr << "----------------------\n";
                assert(false);
            } else {
                utopia::err() << "------ Hessian Check Passed -------\n";
                std::cerr << "error: " << diff << std::endl;
            }

            return ok;
        }

        void spacing(const Scalar &s) { spacing_ = s; }

    private:
        Scalar spacing_;
        bool hessian_from_grad_{true};
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_DIFFCONTROLLER_HPP
