#ifndef UTOPIA_NONLINEAR_ML_BASE_HPP
#define UTOPIA_NONLINEAR_ML_BASE_HPP
#include "utopia_Level.hpp"
#include "utopia_Transfer.hpp"
#include "utopia_MultiLevelBase.hpp"
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_SolutionStatus.hpp"
#include "utopia_MatrixTransfer.hpp"
#include "utopia_MultiLevelEvaluations.hpp"

#include <algorithm>
#include <vector>

namespace utopia {
#define CHECK_NUM_PRECISION_mode


    /**
     * @brief      Base class for all nonlinear multilevel solvers. \n
     *             Takes care of inializing multilevel hierarchy - calls into assembly routines on each level. \n
     *             Different levels are created by interpolation and restriction.\n
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class NonlinearMultiLevelBase : public MultiLevelBase<Matrix, Vector>, public NonLinearSolver<Vector>
    {

    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::Transfer<Matrix, Vector> Transfer;
        typedef utopia::MatrixTransfer<Matrix, Vector> MatrixTransfer;

        typedef utopia::ExtendedFunction<Matrix, Vector> Fun;
        typedef std::shared_ptr<Fun> FunPtr;

        using MultiLevelBase<Matrix, Vector>::set_transfer_operators;

        NonlinearMultiLevelBase(const SizeType & n_levels)
        {
            this->n_levels(n_levels);
        }

        virtual ~NonlinearMultiLevelBase(){}


        virtual void read(Input &in) override
        {
          MultiLevelBase<Matrix, Vector>::read(in);
          NonLinearSolver<Vector>::read(in);
        }

        virtual void print_usage(std::ostream &os) const override
        {
          MultiLevelBase<Matrix, Vector>::print_usage(os);
          NonLinearSolver<Vector>::print_usage(os);
        }


        /**
         * @brief      The solve function for nonlinear multilevel solvers.
         *
         * @param[in]  rhs   Function to be minimized.
         * @param      x_h   The initial guess.
         *
         */
        virtual bool solve( Vector &x_h) = 0;


        /**
         * @brief      Fnction inits functions associated with assemble on each level.
         *
         * @param[in]  level_functions  The level functions
         *
         */
        virtual bool set_functions(const std::vector<FunPtr> &level_functions)
        {
            level_functions_.clear();

            if(this->n_levels() != static_cast<SizeType>(level_functions.size())){
                utopia_error("utopia::NonlinearMultilevelBase:: Number of levels and level_functions do not match. \n");
            }

            level_functions_.insert(level_functions_.begin(), level_functions.begin(), level_functions.end());
            return true;
        }

        /* @brief
         Function initializes transfer  operators.
         Operators need to be ordered FROM COARSE TO FINE.
         *
         * @param[in]  interpolation_operators                The interpolation operators.
         * @param[in]  projection_operators                   The projection operators.
         *
         */
        virtual bool set_transfer_operators(const std::vector<std::shared_ptr<Matrix>> &interpolation_operators,
                                            const std::vector<std::shared_ptr<Matrix>> &projection_operators)
        {
            if(interpolation_operators.size()!=projection_operators.size()){
                utopia_error("utopia::NonlinearMultilevelBase::set_transfer_operators:: Number of interpolation_operators and projection_operators do not match. \n");
            }

            if(this->n_levels() != static_cast<SizeType>(interpolation_operators.size()) + 1){
                utopia_error("utopia::NonlinearMultilevelBase:: Number of levels and transfers do not match. \n");
            }

            this->transfers_.clear();
            for(auto I = interpolation_operators.begin(), P = projection_operators.begin(); I != interpolation_operators.end() && P != projection_operators.end(); ++I, ++P )
                this->transfers_.push_back(std::make_shared<MatrixTransfer>(*I, *P));

            return true;
        }

        /**
         * @brief      Sets the transfer operators.
         *
         * @param[in]  interpolation_operators  The interpolation operators
         * @param[in]  restriction_operators    The restriction operators
         * @param[in]  projection_operators     The projection operators
         *
         */
        virtual bool set_transfer_operators(const std::vector<std::shared_ptr<Matrix>> &interpolation_operators,
                                            const std::vector<std::shared_ptr<Matrix>> &restriction_operators,
                                            const std::vector<std::shared_ptr<Matrix>> &projection_operators)
        {

            if(interpolation_operators.size()!=restriction_operators.size() || interpolation_operators.size()!=projection_operators.size()){
                utopia_error("utopia::NonlinearMultilevelBase::set_transfer_operators:: Number of interpolation_operators and projection_operators do not match. \n");
            }

            if(this->n_levels() != static_cast<SizeType>(interpolation_operators.size()) + 1){
                utopia_error("utopia::NonlinearMultilevelBase:: Number of levels and transfers do not match. \n");
            }

            this->transfers_.clear();
            for(auto I = interpolation_operators.begin(), R = restriction_operators.begin(), P = projection_operators.begin(); I != interpolation_operators.end() && R != restriction_operators.end() &&  P != projection_operators.end(); ++I, ++R, ++P )
                this->transfers_.push_back(std::make_shared<MatrixTransfer>(*I, *R, *P));

            return true;
        }

    protected:

        /**
         * @brief      Function looks up for ids, where we should apply Dirichlet BC and set value to required one
         *
         * @param      fun   The fun
         * @param      x
         *
         */
        virtual bool make_iterate_feasible(Fun & fun, Vector & x)
        {
            const auto &bc_values   = fun.get_eq_constrains_values();
            const auto &bc_ids      = fun.get_eq_constrains_flg(); 

            {
                auto d_bc_ids       = const_device_view(bc_ids);
                auto d_bc_values    = const_device_view(bc_values);

                parallel_transform(x, UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar {
                    Scalar id =  d_bc_ids.get(i);
                    return (id==1.0) ? d_bc_values.get(i) : xi;
                });
            }

            return true;
        }


        /**
         * @brief      Function zeors correction, where we have Dirichlet BC aplied.
         *
         * @param      fun   The fun
         * @param      c     The correction
         */
        virtual bool zero_correction_related_to_equality_constrain(const Fun & fun, Vector & c) 
        {
            fun.zero_contribution_to_equality_constrains(c); 
            return true;
        }


        /**
         * @brief      Function zeors correction, where we have Dirichlet BC aplied.
         *
         * @param      fun   The fun
         * @param      M     matrix
         *
         */
        virtual bool zero_correction_related_to_equality_constrain_mat(const Fun & fun, Matrix & M)
        {
            const std::vector<SizeType> & index = fun.get_indices_related_to_BC(); 
            set_zero_rows(M, index, 1.);

            return true;
        }



        /**
         * @return     Name of solver - to have nice printouts
         */
        virtual std::string name() = 0;

        /**
         * @brief      Init internal memory used for implementation of given multilevel solver
         *
         * @param[in]  fine_local_size  The local size of fine level problem
         */
        virtual void init_memory(const SizeType & fine_local_size) = 0;



        inline Fun &function(const SizeType level)
        {
            assert(level < static_cast<SizeType>(level_functions_.size()));
            assert(level_functions_[level]);

            return *level_functions_[level];
        }

        inline const Fun &function(const SizeType level) const
        {
            assert(level < level_functions_.size());
            assert(level_functions_[level]);

            return *level_functions_[level];
        }


        /**
         * @brief      Writes CSV file with iteration info
         *
         * @param[in]  it_global  The iterator global
         */
        virtual void print_statistics(const SizeType & it_global) override
        {
            std::string path = "log_output_path";
            auto non_data_path = Utopia::instance().get(path);

            if(!non_data_path.empty())
            {
                CSVWriter writer;
                if (mpi_world_rank() == 0)
                {
                    if(!writer.file_exists(non_data_path))
                    {
                        writer.open_file(non_data_path);
                        writer.write_table_row<std::string>({"v_cycles", "time"});
                    }
                    else
                        writer.open_file(non_data_path);

                    writer.write_table_row<Scalar>({Scalar(it_global), this->get_time()});
                    writer.close_file();
                }
            }
        }



    protected:
        std::vector<FunPtr>                      level_functions_;

    };

}

#endif //UTOPIA_NONLINEAR_ML_BASE_HPP

