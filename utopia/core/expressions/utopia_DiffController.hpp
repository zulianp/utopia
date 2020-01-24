//
// Created by Patrick Zulian on 29/05/15.
//

#ifndef UTOPIA_UTOPIA_DIFFCONTROLLER_HPP
#define UTOPIA_UTOPIA_DIFFCONTROLLER_HPP

#include "utopia_FiniteDifference.hpp"
#include "utopia_Wrapper.hpp"

namespace utopia {
    class DiffController {
    public:
        template<class Fun, class Vector, class Matrix>
        bool check(Fun &fun, const Vector &x, const Vector &g, const Matrix &H) {
            return check_grad(fun, x, g) && check_hessian(fun, x, H);
        }

        template<class Fun, class Vector>
        bool check_grad(Fun &fun, const Vector &x, const Vector &g) {
            FiniteDifference<typename Vector::Scalar> fd;
            Vector gfd;
            fd.grad(fun, x, gfd);
            bool ok = approxeq(gfd, g, 1e-2);

            if (!ok) {
                std::cout << "------ Failure -------\n";
                std::cout << "------ Gradient ------\n";
                std::cout << "Expected:\n";
                disp(gfd);

                std::cout << "Actual:\n";
                disp(g);

                std::cout << "----------------------\n";
                assert(false);
            }


            return ok;
        }

        template<class Fun, class Vector, class Matrix>
        bool check_hessian(Fun &fun, const Vector &x, const Matrix &H) {

            using Scalar = typename Traits<Vector>::Scalar;

            FiniteDifference<Scalar> fd;
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

                write("diff_mat.m", diff_mat);
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
    };
}

#endif //UTOPIA_UTOPIA_DIFFCONTROLLER_HPP
