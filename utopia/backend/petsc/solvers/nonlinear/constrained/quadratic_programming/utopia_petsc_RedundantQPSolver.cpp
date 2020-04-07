#include "utopia_petsc_RedundantQPSolver.hpp"
#include <cassert>

namespace utopia {

    void RedundantQPSolver<PetscMatrix, PetscVector>::read(Input &in)
    {
        Super::read(in);
        if(qp_solver_) {
            qp_solver_->read(in);
        }
    }

    bool RedundantQPSolver<PetscMatrix, PetscVector>::valid() const
    {
        return static_cast<bool>(qp_solver_);
    }

    RedundantQPSolver<PetscMatrix, PetscVector>::RedundantQPSolver(const RedundantQPSolver &other)
    {
        if(other.qp_solver_) {
            qp_solver_ = std::shared_ptr<OperatorBasedQPSolver>(other.qp_solver_->clone());
        }
    }

    void RedundantQPSolver<PetscMatrix, PetscVector>::update(const Operator<PetscVector> &A)
    {
        assert(valid());

        auto &&comm = A.comm();
        if(comm.size() == 1) {
            // qp_solver_->update(A);

        } else {
            // Super::update(A);
            auto * mat = dynamic_cast<const PetscMatrix *>(&A);

            if(mat) {
                if(red_.empty()) {
                    red_.init(row_layout(*mat));
                    red_.create_sub_matrix(*mat, redundant_matrix_);
                    red_.create_sub_vector(redundant_sol_);
                    red_.create_sub_vector(redundant_rhs_);

                    redundant_lower_bound_ = std::make_shared<PetscVector>();
                    redundant_upper_bound_ = std::make_shared<PetscVector>();
                    red_.create_sub_vector(*redundant_lower_bound_);
                    red_.create_sub_vector(*redundant_upper_bound_);
                } else {
                    red_.super_to_sub(*mat, redundant_matrix_);
                }
            } else {
                assert(false && "IMPLEMENT ME");
            }
        }
    }

    bool RedundantQPSolver<PetscMatrix, PetscVector>::solve(const Operator<PetscVector> &A, const PetscVector &rhs, PetscVector &sol)
    {
        if(!qp_solver_) {
            return false;
        }

        update(A);

        auto &&comm = A.comm();
        if(comm.size() == 1) {
            qp_solver_->lower_bound() = this->lower_bound();
            qp_solver_->upper_bound() = this->upper_bound();
            return qp_solver_->solve(A, rhs, sol);

        } else {
            red_.super_to_sub(rhs, redundant_rhs_);
            if(!empty(sol)) {
                red_.super_to_sub(sol, redundant_sol_);
            } else {
                redundant_sol_.set(0.0);
            }

            if(this->lower_bound()) {
                red_.super_to_sub(*this->lower_bound(), *redundant_lower_bound_);
                qp_solver_->lower_bound() = redundant_lower_bound_;
            }

            if(this->upper_bound()) {
                red_.super_to_sub(*this->upper_bound(), *redundant_upper_bound_);
                qp_solver_->upper_bound() = redundant_upper_bound_;
            }


            bool ok = qp_solver_->solve(redundant_matrix_, redundant_rhs_, redundant_sol_);
            red_.sub_to_super(redundant_sol_, sol);
            return ok;
        }

    }

}

