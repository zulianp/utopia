//
// Created by Patrick Zulian on 29/05/15.
//

#ifndef UTOPIA_UTOPIA_DIFFCONTROLLER_HPP
#define UTOPIA_UTOPIA_DIFFCONTROLLER_HPP

#include "utopia_FiniteDifference.hpp"

namespace utopia {
    class DiffController {
    public:
        template<class Fun, class Vector, class Matrix>
        bool check(Fun &fun, const Vector &x, const Vector &g, const Matrix &H) {
            return checkGrad(fun, x, g) && checkHessian(fun, x, H);
        }

        template<class Fun, class Vector>
        bool checkGrad(Fun &fun, const Vector &x, const Vector &g) {
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
            }


            return ok;
        }

        template<class Fun, class Vector, class Matrix>
        bool checkHessian(Fun &fun, const Vector &x, const Matrix &H) {

            FiniteDifference<typename Vector::Scalar> fd;
            Matrix Hfd;
            fd.hessian(fun, x, Hfd);

            bool ok = approxeq(Hfd, H, 1e-2);

            if (!ok) {
                std::cout << "------- Failure -------\n";
                std::cout << "------- Hessian -------\n";
                std::cout << "Expected:\n";
                disp(Hfd);

                std::cout << "Actual:\n";
                disp(H);

                std::cout << "----------------------\n";
            }

            return ok;
        }
    };
}

#endif //UTOPIA_UTOPIA_DIFFCONTROLLER_HPP
