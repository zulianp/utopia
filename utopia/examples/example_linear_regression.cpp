#ifndef UTOPIA_LINEAR_REGRESSION
#define UTOPIA_LINEAR_REGRESSION

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "utopia.hpp"
#include "utopia_Base.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Version.hpp"
using namespace std;

// Contains a file whose first columns contain a population/10000 and the profit
// of a foodtrack in a place with that population/10000
const char* file = "./../examples/profit.csv";

namespace utopia {

    template <class Matrix, class Vector>
    class LinearRegression final : public FunctionBase<Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        // Through this constructor, we read the file csv and we save the data in
        // a matrix X and a vector y.
        LinearRegression() : num_of_samples_(96) {
            fstream data_file;
            string line, number;
            int counter = 0;

            Comm comm = Comm::get_default();

            SizeType n_local = num_of_samples_;
            SizeType n_global = n_local * comm.size();  // TODO:: fix

            auto l = layout(comm, n_local, n_global);
            y_.values(l, 0);
            auto y_view = local_view_device(y_);

            // local entries
            SizeType n_local_rows = num_of_samples_;
            SizeType n_local_cols = 2;

            X_.dense(layout(comm, n_local_rows, n_local_cols, n_local_rows * comm.size(), n_local_cols * comm.size()),
                     1.0);

            {
                Write<Matrix> w(X_);
                Range r = row_range(X_);
                SizeType j = r.begin();

                // Check if the csv file opens.
                int i = 0;
                data_file.open(file, ios::in);
                if (!data_file) {
                    cout << "Unable to open file";
                    exit(1);  // terminate with error
                }

                // Saving all the data
                while (!data_file.eof()) {
                    getline(data_file, line);
                    stringstream s(line);

                    while (getline(s, number, ',')) {
                        if (counter == 0 && j != r.end()) {
                            X_.set(j, 1, stof(number));
                            counter = 1;
                            j++;
                        } else if (counter == 1) {
                            y_view.set(i, stof(number));
                            counter = 0;
                            i++;
                        }
                    }
                }
                data_file.close();
            }

            Xt_ = transpose(X_);

            // Initialise the vector fo the squared errors and
            // the hypothesis function
            errors_ = make_unique<Vector>(l, 0.0);
            h_ = make_unique<Vector>(l, 0.0);
        }

        bool value(const Vector& theta, Scalar& cost) const override {
            if (theta.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            // Our initial hypothesis
            *h_ = X_ * theta;

            // Evaluating the hypothesis by summing all the differences before the actual value and expected
            // value
            {
                auto squared_errors_view = local_view_device(*errors_);
                auto h_view = const_local_view_device(*h_);
                auto y_view = const_local_view_device(y_);

                parallel_for(
                    local_range_device(*errors_), UTOPIA_LAMBDA(const SizeType& i) {
                        const auto yi = y_view.get(i);
                        const auto hi = h_view.get(i);
                        squared_errors_view.set(i, (hi - yi) * (hi - yi));
                    });
            }

            // Computing the cost, which should go down at each new iteration of the
            // gradient we are using.
            cost = 1. / (2. * num_of_samples_) * sum(*errors_);

            return true;
        }

        bool gradient(const Vector& theta, Vector& g) const override {
            if (theta.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            // Computing the gradient.
            *h_ = X_ * theta;
            *errors_ = *h_ - y_;

            g = 1. / num_of_samples_ * Xt_ * (*errors_);

            return true;
        }

    private:
        SizeType num_of_samples_;
        Matrix X_;
        Matrix Xt_;
        Vector y_;

        std::unique_ptr<Vector> h_;
        std::unique_ptr<Vector> errors_;
    };

}  // namespace utopia
#endif

template <class Matrix, class Vector>
void test() {
    using namespace utopia;

    // construct our problem
    LinearRegression<Matrix, Vector> newLinearRegression;

    // constructing initial guess
    Vector x;
    x.values(serial_layout(2), 0);

    // // // nonlinear solve
    // GradientDescent<Vector> solver;

    // // set parameters
    // solver.verbose(true);
    // solver.max_it(1500);
    // // solver.dumping_parameter(0.01);
    // auto strategy_to_choose_alpha = std::make_shared<utopia::Backtracking<Vector>>();
    // solver.set_line_search_strategy(strategy_to_choose_alpha);

    // // solve
    // solver.solve(newLinearRegression, x);

    auto hess_approx = std::make_shared<LBFGS<Vector>>(5);
    auto lsolver = std::make_shared<EmptyPrecondMatrixFreeLinearSolver<Vector>>();

    auto precond = hess_approx->build_Hinv_precond();
    lsolver->set_preconditioner(precond);

    QuasiNewton<Vector> solver(hess_approx, lsolver);
    // solve
    solver.verbose(true);
    solver.solve(newLinearRegression, x);

    // display solution
    disp(x, "my_solution");
}

int main(int argc, char** argv) {
    using namespace utopia;

// Checking which backend is active.
#ifdef WITH_PETSC
    using MatrixT = PetscMatrix;
    using VectorT = PetscVector;
#else
#ifdef WITH_TRILINOS
    using MatrixT = TpetraMatrixd;
    using VectorT = TpetraVectord;
#else
    using MatrixT = BlasMatrixd;
    using VectorT = BlasVectord;
#endif
#endif

    Utopia::Init(argc, argv);

    test<MatrixT, VectorT>();

    return Utopia::Finalize();
}