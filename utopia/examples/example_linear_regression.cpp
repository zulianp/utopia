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
 * gives the most reasonable prediction. To do this, we need to minimise the
 * cost function \f$J({\Theta}) : \mathbb{R}^2 \rightarrow  \mathbb{R}\f$, which
 * basically tells us the difference between the actual value y and the value
 * predicted by our hypothesis function.
 */

/**
 * Contains a file whose first column states a population/10000 and the second
 * columns states the profit of a foodtrack in a place with that
 * population/10000.
 */
const char *file = "./../examples/data_for_examples/profit.csv";

/**
 * For this example, we rely on a dataset consisting of s samples. Each sample
 * is defined by an input \f$x_i$ and an output $y_{i}\f$. The input $x_{i}$
 * expresses a population and the output \f$y_{i}\f$ expresses the profit earned
 * by a foodtruck of a certain company. Supposing that we would like to have a
 * foodtruck in a city with a certain population, we would be able to predict
 * the profit for the food truck, thanks to our model parametrized by
 * \f${\Theta} \in \mathbb{R}^{2}\f$. For examples with a population of 43'000,
 * we will be able to write a vector [1, 4.3], - with 1 being the bias, -  and
 * by multiplying this vector by \f${\Theta}\f$ we will obtain an estimation of
 * the profit.
 */

namespace utopia {

    /**
     * @brief A Class to apply linear regression to a one variable sample of data.
     * In our implementation, each population example corresponds to a row in
     * our vector \f$x \in \mathbb{R}^{s}\f$, to which we add an additional columns
     * and we set it to all ones to add the bias. The vector x then becomes a matrix
     * \f$X \in \mathbb{R}^{s \times 2}$\f. Similarly, \f$y\in \mathbb{R}^{s}\f$ is
     * a vector with each row correspondig to an example of the profit. \class
     * \tparam Matrix the matrix type
     * \tparam Vector the vector type
     * @ Filippo Cesana
     * @ August 2020
     */
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

        /**
         * Our cost function is defined as it follows:
         * \f{equation*}{
         *    J(\Theta)
         * =\frac{1}{2s}\mathlarger{\mathlarger{‎‎\sum}}_{i=1}^{s}(h_{\Theta}({\Theta}^{(i)})-y^{(i)})^2
         *  \f}
         * Also, the mean si halved since it makes more convienent to compute the
         * gradient descent: the derivative term will delete the \f$\frac{1}{2}\f$
         * term.
         */

        /** Cost function
         * @param param1 Vector &theta, through which we minimise the cost function
         * @param param2 Sclar & cost: it should go down at each new iteration of the
         * gradient we are using.
         * @return a boolean
         */
        bool value(const Vector &theta, Scalar &cost) const override {
            if (theta.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            /**
             * Initial hypothesis.
             */
            *h_ = X_ * theta;

            /**
             * Evaluating the hypothesis by summing all the differences before the
             * actual value and expected value.
             */
            {
                auto squared_errors_view = local_view_device(*errors_);
                auto h_view = const_local_view_device(*h_);
                auto y_view = const_local_view_device(y_);

                parallel_for(
                    local_range_device(*errors_), UTOPIA_LAMBDA(const SizeType &i) {
                        const auto yi = y_view.get(i);
                        const auto hi = h_view.get(i);
                        squared_errors_view.set(i, (hi - yi) * (hi - yi));
                    });
            }

            /**
             * Computing the cost
             */
            cost = 1. / (2. * num_of_samples_) * sum(*errors_);

            return true;
        }

        /**
         * \f$ h_{\Theta}(x, {\Theta})\f$ is the hypothesis function  \f$h_{\Theta} :
         * \mathbb{R}^{2} \rightarrow \mathbb{R}\f$ and it is given by the following
         * linear model: \f{equation*}{ h_{\Theta}({\Theta}) = {\Theta}X \f}
         */

        /**
         * The \f${\Theta}\f$ values are the paramaters of our model and we try to
         * find optimal values which minimise the cost function \f$J({\Theta})\f$,
         * which is our goal. So we start with some \f$\Theta\f$ and we keep updating
         * them until we find a minimum.
         */

        /**
         * For example, to reach this minimum, we can use a gradient based algorithm
         * which update \f${\Theta}\f$ at each iteration using gradient information.
         */

        /**
         * For \f$j = \{0,1\}\f$, the gradient has the following form:
         */

        /**
         * \f{equation*}{
         * \frac{\partial}{\partial \Theta_{0}}J(\Theta) = \frac{1}{s}
         * \mathlarger{\mathlarger{‎‎\sum}}_{i=1}^{s}(h_{\Theta}({\Theta}^{(i)})-y^{(i)}),
         * \f}
         */

        /**
         * \f{equation*}{
         * \frac{\partial}{\partial \Theta_{1}}J(\Theta) = \frac{1}{s}
         * \mathlarger{\mathlarger{‎‎\sum}}_{i=1}^{s}(h_{\Theta}({\Theta}^{(i)})-y^{(i)})x^{(i)}
         * \f}
         */

        /** Gradient function
         * @param param1 Vector &theta, which we update at each iteration
         * @param param2 second Vector &g the predicted values
         * @return a boolean
         */
        bool gradient(const Vector &theta, Vector &g) const override {
            if (theta.comm().size() > 1) {
                utopia_error("Function is not supported in parallel... \n");
                return false;
            }

            /**
             * Computing the gradient
             */
            *h_ = X_ * theta;
            *errors_ = *h_ - y_;

            g = 1. / num_of_samples_ * Xt_ * (*errors_);

            return true;
        }

    private:
        SizeType num_of_samples_;        /**< Number of samples in the files cv */
        Matrix X_;                       /**< The matrix in which we save the x variables and the bias of one */
        Matrix Xt_;                      /**< The same matrix of above but transposed */
        Vector y_;                       /**< The output y with which we train the model  */
        std::unique_ptr<Vector> h_;      /**< The vector where we saved our hypothesis */
        std::unique_ptr<Vector> errors_; /**< The errors between each hypothesis and the actual output */
    };

}  // namespace utopia
#endif

template <class Matrix, class Vector>
void test() {
    using namespace utopia;

    /** Construct the problem */
    LinearRegression<Matrix, Vector> newLinearRegression;

    /** Constructing the initial guess */
    Vector x;
    x.values(serial_layout(2), 0);

    /** Gradient Descent Solver */
    Chrono GD;
    GD.start();
    GradientDescent<Vector> solverGD;

    /**
     * set parameters
     * solverGD.verbose(true); // uncomment to see full output
     */
    solverGD.max_it(1500);
    /** solverGD.dumping_parameter(0.01); */
    auto strategy_to_choose_alpha = std::make_shared<utopia::Backtracking<Vector>>();
    solverGD.set_line_search_strategy(strategy_to_choose_alpha);

    /** Solve */
    solverGD.solve(newLinearRegression, x);
    GD.stop();

    /** Display the solution*/
    disp(x, "The solution with Gradient");

    /** Quasi Newton Solver */
    Chrono QN;
    QN.start();
    auto hess_approx = std::make_shared<LBFGS<Vector>>(5);
    auto lsolver = std::make_shared<EmptyPrecondMatrixFreeLinearSolver<Vector>>();

    auto precond = hess_approx->build_Hinv_precond();
    lsolver->set_preconditioner(precond);

    QuasiNewton<Vector> solverQN(hess_approx, lsolver);
    /** Constructing the initial guess */
    /** solverQN.verbose(true); */
    solverQN.solve(newLinearRegression, x);
    QN.stop();

    /** Display solution */
    disp(x, "The solution with with QuasiNewton");

    cout << "Gradient Descent took " << GD.get_seconds() << " seconds to converge." << endl;
    cout << "Quasi Newton took " << QN.get_seconds() << " seconds to converge and was " << flush;
    cout << GD.get_seconds() - QN.get_seconds() << " seconds faster than graident descent." << endl;
}

/** This example does not work in parallel.
 * Please compile and run the example with the following command:
 * make utopia_examples && ./examples/example_linear_regression
 */

int main(int argc, char **argv) {
    using namespace utopia;

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