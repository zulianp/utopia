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

/**
 * With linear regression, we are searching for parameters, such that our model
 * gives the most reasonable prediction. To do this, we need to minimise the cost
 * function J(theta), which basically
 * tells us the difference between the actual value y and the value predicted by
 * our hypothesis function.
 * For this example, we rely on a dataset consisting of s samples. Each sample is
 * defined by an input xi and an output yi. The input xi expresses a
 * population and the output yi expresses the profit earned by a foodtruck of
 * a certain company. Supposing that we would like to have a foodtruck in a city
 * with a certain population, we would be able to predict the profit for the food
 * truck, thanks to our model parametrized by theta. For examples with a population
 * of 43'000, we will be able to write a vector [1, 4.3], - with 1 being the bias,
 * -  and by multiplying this vector by theta we will obtain an estimation
 * of the profit. The cost function is defined as:
 */

/**
 * Contains a file whose first columns contain a population/10000 and the profit
 * of a foodtrack in a place with that population/10000.
 */
const char* file = "./../examples/profit.csv";

namespace utopia {

    template <class Matrix, class Vector>
    class LinearRegression final : public FunctionBase<Vector> {
    public:
        using Traits = utopia::Traits<Vector>;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using Comm = typename Traits::Communicator;

        /**
         * Through this constructor, we read the file csv and we save the data in
         * a matrix X and a vector y.
         */
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

            /**
             * local entries
             */
            SizeType n_local_rows = num_of_samples_;
            SizeType n_local_cols = 2;

            X_.dense(layout(comm, n_local_rows, n_local_cols, n_local_rows * comm.size(), n_local_cols * comm.size()),
                     1.0);

            {
                Write<Matrix> w(X_);
                Range r = row_range(X_);
                SizeType j = r.begin();

                int i = 0;
                data_file.open(file, ios::in);
                if (!data_file) {
                    cout << "Unable to open file";
                    exit(1);
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

            /**
             * Initialise the vector fo the squared errors and
             * the hypothesis function.
             */
            errors_ = make_unique<Vector>(l, 0.0);
            h_ = make_unique<Vector>(l, 0.0);
        }

        // \f\begin{equation*}
        // 	J(\Theta)
        // =\frac{1}{2s}\mathlarger{\mathlarger{‎‎\sum}}_{i=1}^{s}(h_{\Theta}({\Theta}^{(i)})-y^{(i)})^2
        // \end{equation*} \f
        bool value(const Vector& theta, Scalar& cost) const override {
            if (theta.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            // Our initial hypothesis
            *h_ = X_ * theta;

            ////////////////////////////////////////////////////////////////////////////////
            // Evaluating the hypothesis by summing all the differences before the actual
            // value and expected value.
            ////////////////////////////////////////////////////////////////////////////////
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

            ////////////////////////////////////////////////////////////////////////////////
            // Computing the cost, which should go down at each new iteration of the
            // gradient we are using.
            ////////////////////////////////////////////////////////////////////////////////
            cost = 1. / (2. * num_of_samples_) * sum(*errors_);

            return true;
        }

        ////////////////////////////////////////////////////////////////////////////////
        // The theta values are the paramater of our model and we try to find optimal
        // some theta and we keep updating them until we find a minumum.
        // For example, to reach
        // this minimum, we can use a gradient based algorithm which update theta at each
        // iteration using gradient information.values which minimise the cost function
        // J(theta) which is our goal. So we start with
        ////////////////////////////////////////////////////////////////////////////////
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

    ////////////////////////////////////////////////////////////////////////////////
    // construct our problem
    ////////////////////////////////////////////////////////////////////////////////
    LinearRegression<Matrix, Vector> newLinearRegression;

    ////////////////////////////////////////////////////////////////////////////////
    // constructing initial guess
    Vector x;
    ////////////////////////////////////////////////////////////////////////////////
    x.values(serial_layout(2), 0);

    ////////////////////////////////////////////////////////////////////////////////
    // Gradient Descent solver
    ////////////////////////////////////////////////////////////////////////////////
    Chrono GD;
    GD.start();
    GradientDescent<Vector> solverGD;

    ////////////////////////////////////////////////////////////////////////////////
    // set parameters
    // solverGD.verbose(true); // uncomment to see full output
    ////////////////////////////////////////////////////////////////////////////////
    solverGD.max_it(1500);
    ////////////////////////////////////////////////////////////////////////////////
    // solverGD.dumping_parameter(0.01);
    ////////////////////////////////////////////////////////////////////////////////
    auto strategy_to_choose_alpha = std::make_shared<utopia::Backtracking<Vector>>();
    solverGD.set_line_search_strategy(strategy_to_choose_alpha);

    ////////////////////////////////////////////////////////////////////////////////
    // solve
    ////////////////////////////////////////////////////////////////////////////////
    solverGD.solve(newLinearRegression, x);
    GD.stop();

    ////////////////////////////////////////////////////////////////////////////////
    // display the solution
    ////////////////////////////////////////////////////////////////////////////////
    disp(x, "my_solution with Gradient");

    ////////////////////////////////////////////////////////////////////////////////
    // Quasi Newton solver
    ////////////////////////////////////////////////////////////////////////////////
    Chrono QN;
    QN.start();
    auto hess_approx = std::make_shared<LBFGS<Vector>>(5);
    auto lsolver = std::make_shared<EmptyPrecondMatrixFreeLinearSolver<Vector>>();

    auto precond = hess_approx->build_Hinv_precond();
    lsolver->set_preconditioner(precond);

    QuasiNewton<Vector> solverQN(hess_approx, lsolver);
    // solve
    // solverQN.verbose(true);
    solverQN.solve(newLinearRegression, x);
    QN.stop();

    // display solution
    disp(x, "my_solution with QuasiNewton");

    cout << "Gradient Descent took " << GD.get_seconds() << " seconds to converge." << endl;
    cout << "Quasi Newton took " << QN.get_seconds() << " seconds to converge and was " << flush;
    cout << GD.get_seconds() - QN.get_seconds() << " seconds faster than graident descent." << endl;
}

// This example does not work in parallel.
// Please run the example with the following command:
// make -j 4  utopia_examples && ./examples/example_linear_regression
int main(int argc, char** argv) {
    using namespace utopia;

#ifdef UTOPIA_WITH_PETSC
    using MatrixT = PetscMatrix;
    using VectorT = PetscVector;
#else
#ifdef UTOPIA_WITH_TRILINOS
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