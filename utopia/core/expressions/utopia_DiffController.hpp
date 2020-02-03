#ifndef UTOPIA_UTOPIA_DIFFCONTROLLER_HPP
#define UTOPIA_UTOPIA_DIFFCONTROLLER_HPP

#include "utopia_FiniteDifference.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Input.hpp"

namespace utopia {

    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class DiffController : public Configurable {
    public:
        using Scalar = typename Traits<Vector>::Scalar;

        DiffController(const Scalar spacing = 1e-5)
        : spacing_(spacing)
        {}

        void read(Input &in) override
        {
            in.get("spacing", spacing_);
        }

        template<class Fun>
        bool check(Fun &fun, const Vector &x, const Vector &g, const Matrix &H) {
            return check_grad(fun, x, g) && check_hessian(fun, x, H);
        }

        template<class Fun>
        bool check_grad(Fun &fun, const Vector &x, const Vector &g) {
            FiniteDifference<typename Vector::Scalar> fd(spacing_);
            Vector gfd;
            fd.grad(fun, x, gfd);
            bool ok = approxeq(gfd, g, 1e-2);

            if (!ok) {
                std::cout << "------ Failure -------\n";
                std::cout << "------ Gradient ------\n";
                std::cout << "Expected:\n";
                rename("g_fd", gfd);
                disp(gfd);

                std::cout << "Actual:\n";
                disp(g);

                std::cout << "----------------------\n";
                assert(false);
            }


            return ok;
        }

        template<class Fun>
        bool check_hessian(Fun &fun, const Vector &x, const Matrix &H) {

            FiniteDifference<Scalar> fd(spacing_);
            Matrix Hfd = H;
            Hfd *= 0.0;

            if(!fd.hessian(fun, x, Hfd)) {
                return true;
            }

            Matrix diff_mat = H - Hfd;
            Scalar diff = norm1(diff_mat);

            bool ok = diff < 1e-2;

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
                write("H.m",   H);

                // std::cerr << "Expected:\n";
                // disp(Hfd);

                // std::cerr << "Actual:\n";
                // disp(H);

                std::cerr << "----------------------\n";
                assert(false);
            }

            return ok;
        }

        void spacing(const Scalar &s)
        {
            spacing_ = s;
        }

    private:
        Scalar spacing_;
    };
}

#endif //UTOPIA_UTOPIA_DIFFCONTROLLER_HPP
