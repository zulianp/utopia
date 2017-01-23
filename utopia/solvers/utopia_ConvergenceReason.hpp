/*
* @Author: Alena Kopanicakova
* @Date:   2016-07-31
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-10-11
*/

#ifndef UTOPIA_SOLVER_CONVERGENCE_REASON
#define UTOPIA_SOLVER_CONVERGENCE_REASON


#include <iomanip>
#include <limits> 
#include <chrono>   
    
namespace utopia 
{
    /**
     * @brief      The class helping to print-out information about solver: initialization messages, prinitning iteration status, time-stats and exit/convergance messages. 
     *             It also helps pass solution and informations about solve back into FEM packages. 
     */
    class  ConvergenceReason
    {
        typedef double Scalar;

    public:
        
        // TODO:  put all reasons here ... 
        static void exitMessage(const long &num_it, const int &convergence_reason)
        {
            std::cout << std::endl;
            if(convergence_reason == DIVERGED_MAX_IT )
            {
                std::cerr << "\033[1;31m [Error] Maximum number of iteration reached (" << num_it << "). \033[0m\n"; 
            }
            else if(convergence_reason == CONVERGED_SNORM_RELATIVE)
            {
                std::cout << "\033[1;32m  LinearSolver terminated at iteration " << num_it << ", (step_length< tol). \033[0m\n";
            }
            else if(convergence_reason == CONVERGED_TR_DELTA)
            {
                std::cout << "\033[1;32m  LinearSolver terminated at iteration " << num_it << ", (radius < 1e-12). \033[0m\n";
            }
            else if(convergence_reason == CONVERGED_SNORM_RELATIVE)
            {
                std::cout << "\033[1;32m  LinearSolver terminated at iteration " << num_it << ", no more energy reduction. \033[0m\n";
            }
            else if(convergence_reason < 0)
            {
                std::cerr << "\033[1;31m [Error] LinearSolver stopped at iteration " << num_it << " . \033[0m\n";
            }
            else
            {
                std::cout << "\033[1;32m  LinearSolver converged in " << num_it << " iterations.\033[0m\n";
            }
        }
    
        // success 
        static const int CONVERGED_FNORM_ABS          = 2;   /* ||g|| < atol */
        static const int CONVERGED_FNORM_RELATIVE     = 3;   /* ||g|| < rtol*||g_0|| */
        static const int CONVERGED_SNORM_RELATIVE     = 4;   /* Newton computed step size small; || delta x || < stol || x ||*/
        static const int CONVERGED_ITS                = 5;   /* maximum iterations reached */
        static const int CONVERGED_TR_DELTA           = 7;

        // fail
        static const int DIVERGED_FUNCTION_DOMAIN     = -1;  /* the new x location passed the function is not in the domain of F */
        static const int DIVERGED_FUNCTION_COUNT      = -2;
        static const int DIVERGED_LINEAR_SOLVE        = -3;  /* the linear solve failed */
        static const int DIVERGED_FNORM_NAN           = -4;
        static const int DIVERGED_MAX_IT              = -5;
        static const int DIVERGED_LINE_SEARCH         = -6;  /* the line search failed */
        static const int DIVERGED_INNER               = -7;  /* inner solve failed */
        static const int DIVERGED_LOCAL_MIN           = -8;  /* || J^T b || is small, implies converged to local minimum of F() */
        static const int CONVERGED_ITERATING          =  0;


    private: 
        ConvergenceReason()
        {

        }

    };

}
#endif //UTOPIA_SOLVER_CONVERGENCE_REASON
