/*
* @Author: alenakopanicakova
* @Date:   2016-10-03
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-10-14
*/
#ifndef UTOPIA_DIRECT_SOLVER_HPP
#define UTOPIA_DIRECT_SOLVER_HPP 


#include <iostream>
#include <algorithm>
#include "utopia_LinearSolver.hpp"
#include "utopia_Preconditioner.hpp"

namespace utopia 
{   
    template<typename Matrix, typename Vector>
    class DirectSolver : public LinearSolver<Matrix, Vector> {
    public:
        virtual bool apply(const Vector &rhs, Vector &sol) override
        {
            std::cerr << "[Warning] Subclass of DirectSolver did not override apply(.,.). Using solve..." << std::endl;
            std::cerr << "[Task]    You need to fix your direct solver and override apply(.,.) instead of solve(.,.,.)." << std::endl;
            assert(false && "override this method");
            return this->solve(*this->get_operator(), rhs, sol);
        }
    };

}

#endif //UTOPIA_DIRECT_SOLVER_HPP
