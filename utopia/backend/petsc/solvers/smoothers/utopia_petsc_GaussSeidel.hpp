/*
* @Author: alenakopanicakova
* @Date:   2016-04-04
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-04-21
*/

#ifndef UTOPIA_PETSC_GS_HPP
#define UTOPIA_PETSC_GS_HPP

#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"


#include <petscpc.h>
#include <petscksp.h>

// extern "C" 
// {
#include "petscmat.h"
#include "petscvec.h"
// }


namespace utopia {
      template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
  class GaussSeidel {};




        /**
         * @brief      Wrapper for PETSC implementation of SOR. 
         *             Be aware, that this function doesn't run in parallel. 
         *  
         * @tparam     Matrix  
         * @tparam     Vector  
         */
        template<class Matrix, class Vector>
  class GaussSeidel<Matrix, Vector, PETSC> : public IterativeSolver<Matrix, Vector>, 
  public Smoother<Matrix, Vector>
  {
    typedef UTOPIA_SCALAR(Vector)                   Scalar;
    typedef UTOPIA_SIZE_TYPE(Vector)                SizeType;
    typedef utopia::IterativeSolver<Matrix, Vector> Solver;
    typedef utopia::Smoother<Matrix, Vector>        Smoother;

public:
    GaussSeidel(const Parameters params = Parameters()) 
    { 
        set_parameters(params); 
    }

        /**
         * @brief      Smoothing of GS from Petsc. Currently we are using symmetric block GS (builds block jacobi and on blocks calls GS). 
         *
         * @param[in]  A     The stiffness matrix. 
         * @param[in]  rhs   The right hand side. 
         * @param      x     The solution. 
         */
    bool smooth(const Matrix &A, const Vector &rhs, Vector &x) override
    {        
        MatSOR( raw_type(A), 
            raw_type(rhs), 
            1, 
                    //  SOR_FORWARD_SWEEP,
            SOR_LOCAL_SYMMETRIC_SWEEP,  
            0, 
            this->sweeps(), 
            this->sweeps(), 
            raw_type(x)); 

        return true; 
    }

        /**
         * @brief      Solving system with Gauss-Seidel method. 
         *
         * @param[in]  A     The stifness matrix. 
         * @param[in]  rhs   The right hand side.
         * @param      x     The solution. 
         */
    bool solve(const Matrix &A, const Vector &rhs, Vector &x) override
    {
        SizeType it = 0; 
        Scalar r_norm = 9999; 
        this->init_solver("Petsc Gauss-Seidel", {"it. ", "||r||" }); 

        while(it < this->max_it() && r_norm > this->rtol()) {
            MatSOR(
                raw_type(A), 
                raw_type(rhs), 
                1, 
                // SOR_FORWARD_SWEEP,
                SOR_LOCAL_SYMMETRIC_SWEEP,    // parallel implementation - builds block jacobi and on blocks it calls GS 
                0, 
                this->sweeps(), 
                this->sweeps(), 
                raw_type(x)); 

            r_norm = norm2(A*x - rhs);

            it += this->sweeps();

            if(this->verbose())
                PrintInfo::print_iter_status(it, {r_norm}); 
        }

        return true; 
    }


    void set_parameters(const Parameters params) override
    {
        Smoother::set_parameters(params); 
        Solver::set_parameters(params); 
    }



        /**
         * @brief      Nonlinear smoothing. First it calls GS smoother and then it project constrains. 
         * !!!! THIS DOES NOT GIVE CORRECT RESULTS IN A SIMULATION THE PROJECTION CANNOT HAPPEN OUTSIDE A GS STEP
         *
         * @param[in]  A            The stiffness matrix. 
         * @param[in]  rhs          The right hand side. 
         * @param      x            The solution. 
         * @param[in]  ub           The upper bound. 
         * @param[in]  lb           The lower bound. 
         * @param      active_set   The vector containing indices of active set. 
         */
        // bool nonlinear_smooth(const Matrix &A, const Vector &rhs, const Vector& ub, const Vector& lb, Vector &x, std::vector<SizeType>& active_set) override
        // {

        //     smooth(A, rhs, x);
        //     project_constraints(ub, lb, x, active_set); 
        //     return true; 

        // }


private:
        /**
         * @brief      Function projects constraints, such that \f$ lb < x < ub \f$
         *
         * @param[in]  ub           The upper bound. 
         * @param[in]  lb           The lower bound. 
         * @param      x            The solution. 
         * @param      active_set   The vector containing indices of active set. 
         * !!!! THIS DOES NOT GIVE CORRECT RESULTS IN A SIMULATION THE PROJECTION CANNOT HAPPEN OUTSIDE A GS STEP
         */
    bool project_constraints(const Vector& ub, const Vector& lb, Vector &x, std::vector<SizeType>& active_set)
    {
        Vector x_0 = x; 
        {
            Read<Vector> r_ub(ub), r_lb(lb);
            each_transform(x_0, x, [&ub, &lb, &active_set](const SizeType i, const Scalar entry) -> double  
            { 
                Scalar ui = ub.get(i), li = lb.get(i); 
                if(entry > ui && entry < li)
                {
                    active_set.push_back(i);
                    return (std::max(li, std::min(ui, entry)));
                }
                else
                    return entry;   
            }    );
        }
        return true;
    }

};

}

#endif //UTOPIA_PETSC_GS_HPP

