#ifndef UTOPIA_PETSC_NEWTON_HPP
#define UTOPIA_PETSC_NEWTON_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_Newton.hpp"

#include <iomanip>
#include <limits>


namespace utopia
{

    template<class Matrix, class Vector>
    class Newton<Matrix, Vector, PETSC_EXPERIMENTAL> final: public SNESSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                   Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                SizeType;

        typedef utopia::SNESSolver<Matrix, Vector>      SNESSolver;
        typedef utopia::LSStrategy<Vector>              LSStrategy; 
        typedef utopia::LinearSolver<Matrix, Vector>    LinearSolver;


        public:
        Newton( const std::shared_ptr <LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(), 
                const Scalar & alpha = 1.0, 
                const SizeType & order = 3, 
                const std::vector<std::string> ls_types    = {"basic", "bt", "l2", "cp", "nleqerr"}) :   
                    
                    SNESSolver(linear_solver), alpha_(alpha), line_search_types_(ls_types)
        { 
            line_search_type_ = line_search_types_.at(0); 
            this->set_snes_type("newtonls"); 

            line_search_order_ = (order < 4 && order > 0) ?  order : 3; 
        }

        void line_search_type(const std::string & ls_type)
        {
            line_search_type_ = in_array(ls_type, line_search_types_) ? ls_type : line_search_types_.at(0);; 
        }

        void line_search_order(const SizeType & order)
        {
            line_search_order_ = (order < 4 && order > 0) ?  order : 3; 
        }

        void damping_parameter(const Scalar & alpha)
        {
            alpha_ = alpha; 
        }


        void read(Input &in) override
        {
            SNESSolver::read(in); 

            in.get("damping_parameter", alpha_);
            in.get("line_search_order", line_search_order_);
            in.get("line_search_type", line_search_type_);
        }


        void print_usage(std::ostream &os) const override
        {
            SNESSolver::print_usage(os); 

            this->print_param_usage(os, "damping_parameter", "double", "Constant damping parameter.", "1"); 
            this->print_param_usage(os, "line_search_order", "int", "Interpolation order used inside of linea-search-strategy.", "3"); 
            this->print_param_usage(os, "line_search_type", "string", "Type of line-search to be used.", "basic"); 
        }        



    protected: 
        void set_snes_options(SNES & snes,  const Scalar & atol     = SNESSolver::atol(), 
                                            const Scalar & rtol     = SNESSolver::rtol(), 
                                            const Scalar & stol     = SNESSolver::stol(), 
                                            const SizeType & max_it = SNESSolver::max_it()) override 
        {
            SNESSolver::set_snes_options(snes, atol, rtol, stol, max_it); 
            
            SNESLineSearch linesearch; 
            SNESGetLineSearch(snes, &linesearch);
            SNESLineSearchSetType(linesearch,   line_search_type_.c_str() ); 
            SNESLineSearchSetOrder(linesearch, line_search_order_); 
            SNESLineSearchSetDamping(linesearch, alpha_); 
        }


    private:
        Scalar alpha_; 
        const std::vector<std::string> line_search_types_; 
        std::string line_search_type_; 
        SizeType line_search_order_; 

    };



}
#endif //UTOPIA_PETSC_NEWTON_HPP
