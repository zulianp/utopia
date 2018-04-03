/*
* @Author: Alena Kopanicakova
* @Date:   2016-07-31
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-01-31
*/

#ifndef UTOPIA_SOLVER_CONVERGENCE_REASON
#define UTOPIA_SOLVER_CONVERGENCE_REASON


#include <iomanip>
#include <limits> 
#include <chrono> 
#include <string>  
    
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
            if(convergence_reason == CONVERGED_ITERATING)
            {
                std::cout << "\033[1;32m  LinearSolver converged in " << num_it << " iterations.\033[0m\n";
            }
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
                std::cerr << "\033[1;31m [Error] LinearSolver diverged at iteration " << num_it << " reason = " << diverged_reason_string(convergence_reason) << ". \033[0m\n";
            }
            else
            {
                std::cerr << "\033[1;31m LinearSolver converged at iteration " << num_it << " for reason = " << convergence_reason << " \033[0m\n";
            }
        }

        static void exitMessage_nonlinear(const long &num_it, const int &convergence_reason)
        {
            std::cout << std::endl;
            if(convergence_reason == CONVERGED_ITERATING)
            {
                std::cout << "\033[1;32m  NonlinearSolver converged in " << num_it << " iterations.\033[0m\n";
            }
            else if(convergence_reason == DIVERGED_MAX_IT )
            {
                std::cerr << "\033[1;31m [Error] Nonlinear solver: Maximum number of iteration reached (" << num_it << "). \033[0m\n"; 
            }
            else if(convergence_reason == CONVERGED_SNORM_RELATIVE)
            {
                std::cout << "\033[1;32m  NonlinearSolver terminated at iteration " << num_it << ", (step_length< tol). \033[0m\n";
            }
            else if(convergence_reason == CONVERGED_TR_DELTA)
            {
                std::cout << "\033[1;32m  NonlinearSolver terminated at iteration " << num_it << ", (radius < 1e-12). \033[0m\n";
            }
            else if(convergence_reason == CONVERGED_SNORM_RELATIVE)
            {
                std::cout << "\033[1;32m  NonlinearSolver terminated at iteration " << num_it << ", no more energy reduction. \033[0m\n";
            }
            else if(convergence_reason == CONVERGED_FNORM_RELATIVE)
            {
                std::cout << "\033[1;32m  NonlinearSolver terminated at iteration " << num_it << ", (|| F || < atol). \033[0m\n";
            }            
            else if(convergence_reason < 0)
            {
                std::cerr << "\033[1;31m [Error] NonlinearSolver stopped at iteration " << num_it << " . \033[0m\n";
            }
            else
            {
                 std::cerr << "\033[1;31m NonlinearSolver converged at iteration " << num_it << " for reason = " << convergence_reason << " \033[0m\n";
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

        static std::string diverged_reason_string(const int reason)
        {
            if(reason >= 12) {
                return "UNDEFINED";
            }


            static const std::string string_reason[12] =
            {
                "KSP_CONVERGED_ITERATING",
                "UNDEFINED",
                "KSP_DIVERGED_NULL",
                "KSP_DIVERGED_ITS",
                "KSP_DIVERGED_DTOL",
                "KSP_DIVERGED_BREAKDOWN",
                "KSP_DIVERGED_BREAKDOWN_BICG",
                "KSP_DIVERGED_NONSYMMETRIC",
                "KSP_DIVERGED_INDEFINITE_PC",
                "KSP_DIVERGED_NANORINF",
                "KSP_DIVERGED_INDEFINITE_MAT",
                "KSP_DIVERGED_PCSETUP_FAILED"
            };   

            return string_reason[std::abs(reason)];
        }

    private: 
        ConvergenceReason()
        {

        }

    };

}
#endif //UTOPIA_SOLVER_CONVERGENCE_REASON
