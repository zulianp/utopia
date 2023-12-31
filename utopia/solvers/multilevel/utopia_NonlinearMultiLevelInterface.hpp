#ifndef UTOPIA_NONLINEAR_ML_INTERFACE_HPP
#define UTOPIA_NONLINEAR_ML_INTERFACE_HPP
#include "utopia_Core.hpp"
#include "utopia_Function.hpp"
#include "utopia_IPRTransfer.hpp"
#include "utopia_Level.hpp"
#include "utopia_MultiLevelBase.hpp"
#include "utopia_MultiLevelEvaluations.hpp"
#include "utopia_SolutionStatus.hpp"
#include "utopia_Transfer.hpp"

#include <algorithm>
#include <vector>

namespace utopia {
#define CHECK_NUM_PRECISION_mode

    /**
     * @brief      Base class for all nonlinear multilevel solvers. \n
     *             Takes care of inializing multilevel hierarchy - calls into
     * assembly routines on each level. \n Different levels are created by
     * interpolation and restriction.\n
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector>
    class NonlinearMultiLevelInterface : public MultiLevelBase<Matrix, Vector> {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;

        typedef utopia::Transfer<Matrix, Vector> Transfer;
        typedef utopia::IPRTransfer<Matrix, Vector> IPRTransfer;

        typedef utopia::ExtendedFunction<Matrix, Vector> Fun;
        using FunPtr = std::shared_ptr<Fun>;

        using MultiLevelBase<Matrix, Vector>::set_transfer_operators;

        NonlinearMultiLevelInterface(const SizeType &n_levels) { this->n_levels(n_levels); }

        ~NonlinearMultiLevelInterface() override = default;

        void read(Input &in) override { MultiLevelBase<Matrix, Vector>::read(in); }

        void print_usage(std::ostream &os) const override { MultiLevelBase<Matrix, Vector>::print_usage(os); }

        /**
         * @brief      Fnction inits functions associated with assemble on each level.
         *
         * @param[in]  level_functions  The level functions
         *
         */
        virtual bool set_functions(const std::vector<FunPtr> &level_functions) {
            level_functions_.clear();
            local_level_layouts_.clear();

            if (this->n_levels() != static_cast<SizeType>(level_functions.size())) {
                utopia_error(
                    "utopia::NonlinearMultilevelBase:: Number of levels and "
                    "level_functions do not match. \n");
            }

            level_functions_.insert(level_functions_.begin(), level_functions.begin(), level_functions.end());

            for (auto l = 0; l < this->n_levels(); l++) {
                local_level_layouts_.push_back(level_functions_[l]->layout());
            }

            return true;
        }

        /* @brief
         Function initializes transfer  operators.
         Operators need to be ordered FROM COARSE TO FINE.
         *
         * @param[in]  interpolation_operators                The interpolation
         operators.
         * @param[in]  projection_operators                   The projection
         operators.
         *
         */
        virtual bool set_transfer_operators(const std::vector<std::shared_ptr<Matrix>> &interpolation_operators,
                                            const std::vector<std::shared_ptr<Matrix>> &projection_operators) {
            if (interpolation_operators.size() != projection_operators.size()) {
                utopia_error(
                    "utopia::NonlinearMultilevelBase::set_transfer_operators:: Number of "
                    "interpolation_operators and "
                    "projection_operators do not match. \n");
            }

            if (this->n_levels() != static_cast<SizeType>(interpolation_operators.size()) + 1) {
                utopia_error(
                    "utopia::NonlinearMultilevelBase:: Number of levels and transfers do "
                    "not match. \n");
            }

            this->transfers_.clear();
            for (auto I = interpolation_operators.begin(), P = projection_operators.begin();
                 I != interpolation_operators.end() && P != projection_operators.end();
                 ++I, ++P)
                this->transfers_.push_back(std::make_shared<IPRTransfer>(*I, *P));

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
                                            const std::vector<std::shared_ptr<Matrix>> &projection_operators) {
            if (interpolation_operators.size() != restriction_operators.size() ||
                interpolation_operators.size() != projection_operators.size()) {
                utopia_error(
                    "utopia::NonlinearMultilevelBase::set_transfer_operators:: Number of "
                    "interpolation_operators and "
                    "projection_operators do not match. \n");
            }

            if (this->n_levels() != static_cast<SizeType>(interpolation_operators.size()) + 1) {
                utopia_error(
                    "utopia::NonlinearMultilevelBase:: Number of levels and transfers do "
                    "not match. \n");
            }

            this->transfers_.clear();
            for (auto I = interpolation_operators.begin(),
                      R = restriction_operators.begin(),
                      P = projection_operators.begin();
                 I != interpolation_operators.end() && R != restriction_operators.end() &&
                 P != projection_operators.end();
                 ++I, ++R, ++P)
                this->transfers_.push_back(std::make_shared<IPRTransfer>(*I, *R, *P));

            return true;
        }

        /**
         * @brief      Function looks up for ids, where we should apply Dirichlet BC
         * and set value to required one
         *
         * @param      fun   The fun
         * @param      x
         *
         */
        virtual bool make_iterate_feasible(Fun &fun, Vector &x) {
            const auto &bc_values = fun.get_eq_constrains_values();
            const auto &bc_ids = fun.get_eq_constrains_flg();

            {
                auto d_bc_ids = const_local_view_device(bc_ids);
                auto d_bc_values = const_local_view_device(bc_values);
                auto x_view = local_view_device(x);

                parallel_for(
                    local_range_device(x), UTOPIA_LAMBDA(const SizeType &i) {
                        Scalar id = d_bc_ids.get(i);
                        auto xi = x_view.get(i);

                        x_view.set(i, (id == 1.0) ? d_bc_values.get(i) : xi);
                    });
            }

            return true;
        }

    protected:
        /**
         * @brief      Function zeors correction, where we have Dirichlet BC aplied.
         *
         * @param      fun   The fun
         * @param      c     The correction
         */
        virtual bool zero_correction_related_to_equality_constrain(const Fun &fun, Vector &c) {
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
        virtual bool zero_correction_related_to_equality_constrain_mat(const Fun &fun, Matrix &M) {
            const std::vector<SizeType> &index = fun.get_indices_related_to_BC();
            set_zero_rows(M, index, 1.);

            return true;
        }

        /**
         * @return     Name of solver - to have nice printouts
         */
        virtual std::string name() = 0;

        /**
         * @brief      Init internal memory used for implementation of given
         * multilevel solver
         *
         * @param[in]  fine_local_size  The local size of fine level problem
         */
        virtual void init_memory() = 0;

        inline Fun &function(const SizeType level) {
            assert(level < static_cast<SizeType>(level_functions_.size()));
            assert(level_functions_[level]);

            return *level_functions_[level];
        }

        inline const Fun &function(const SizeType level) const {
            assert(level < level_functions_.size());
            assert(level_functions_[level]);

            return *level_functions_[level];
        }

        const Layout &local_layouts(const SizeType &level) const { return local_level_layouts_[level]; }

        const std::vector<Layout> &local_level_layouts() const { return local_level_layouts_; }

        const std::vector<FunPtr> &level_functions() { return level_functions_; }

    protected:
        std::vector<FunPtr> level_functions_;
        std::vector<Layout> local_level_layouts_;
    };

}  // namespace utopia

#endif  // UTOPIA_NONLINEAR_ML_INTERFACE_HPP
