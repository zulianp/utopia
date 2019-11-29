#ifndef UTOPIA_BACKTRACKING_HPP
#define UTOPIA_BACKTRACKING_HPP
#include "utopia_LS_Strategy.hpp"
#include "utopia_PrintInfo.hpp"


namespace utopia
{

    /**
     * @brief      This class implements very simple backtracking, without any interpolation with Wolfe SDC.
     *             It is mainly made for learning and testing purposes. \n
     *             For more serious implementation check Backtracing class.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Vector>
    class SimpleBacktracking final : public LSStrategy<Vector>
    {
        typedef UTOPIA_SCALAR(Vector) Scalar;

    public:
        SimpleBacktracking():  LSStrategy<Vector>()

        {

        }

        bool get_alpha(LeastSquaresFunctionBase<Vector> &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) override
        {
            return get_alpha_aux(fun, g, x, d, alpha);
        }


        bool get_alpha(FunctionBase<Vector> &fun, const Vector &g, const Vector& x, const Vector &d, Scalar &alpha) override
        {
            return get_alpha_aux(fun, g, x, d, alpha);
        }

        template<class FunctionT>
        bool get_alpha_aux(FunctionT &fun, const Vector &g, const Vector& x, const Vector &p_k, Scalar &alpha_k)
        {
            Scalar E_k, E_k1, g_p;

            fun.value(x, E_k);
            alpha_k = 1.0;
            g_p =  dot(g, p_k);

            E_k1 = E_k;

            x_k = x + alpha_k * p_k;
            fun.value(x_k, E_k1);

            SizeType it = 0;

            if(this->verbose())
                PrintInfo::print_init("SIMPLE_BACKTRACKING_LS_INNER_ITERATIONS", {" it. ", "|| E_k1 ||"});

            // Wolfe conditions
            while( E_k1 >(E_k + this->c1() * alpha_k * g_p) && it < this->max_it()  && alpha_k > this->alpha_min())
            {
                x_k = x + alpha_k * p_k;
                fun.value(x_k, E_k1);
                it++;
                alpha_k *= this->rho();
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {E_k1});

            }

           // std::cout<<"it:  "<< it << "  \n";
            return true;
        }


    private:
        Vector x_k; 

    };

}


#endif //UTOPIA_BACKTRACKING_HPP