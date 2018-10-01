#ifndef UTOPIA_TRILINOS_LINEAR_SOLVER_FACTORY_HPP
#define UTOPIA_TRILINOS_LINEAR_SOLVER_FACTORY_HPP


#include "utopia_Traits.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverFactory.hpp"

#include "utopia_DirectSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_ConjugateGradient.hpp"

#include <map>
#include <string>
#include <memory>

namespace utopia
{

template<typename Matrix, typename Vector>
class LinearSolverFactory<Matrix, Vector, TRILINOS>
    {
    public:
        typedef std::shared_ptr< LinearSolver<Matrix, Vector> >  LinearSolverPtr;
        std::map<std::string, LinearSolverPtr> solvers_;
        ///See new implementation in utopia_petsc_LinearSolverFactory.hpp

        inline static LinearSolverPtr new_linear_solver(const SolverTag &tag)
            {
            auto it = instance().solvers_.find(tag);
            if(it == instance().solvers_.end())
                {
                //std::cout<<"LinearSolver not avaialble, solve with CG.  \n";   // TODO fix tests and put back
                return std::make_shared<ConjugateGradient<Matrix, Vector> >();
                }
            else
                {
                return it->second;
                }
            }


    private:

        inline static const LinearSolverFactory &instance()
            {
            static LinearSolverFactory instance_;
            return instance_;

            }

            LinearSolverFactory()
            {
                init();
            }


        void init()
            {
//TODO
            }

        /*        virtual bool apply(const Vectord &rhs, Vectord &sol)
                    {
        //TODO
                    }*/

    };
}
#endif //UTOPIA_TRILINOS_LINEAR_SOLVER_FACTORY_HPP

