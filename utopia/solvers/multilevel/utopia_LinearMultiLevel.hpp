#ifndef UTOPIA_LINEAR_MULTI_LEVELHPP
#define UTOPIA_LINEAR_MULTI_LEVELHPP

#include "utopia_MatrixTransfer.hpp"
#include "utopia_MultiLevelBase.hpp"
#include "utopia_MultiLevelMask.hpp"
#include "utopia_Path.hpp"
#include "utopia_Recorder.hpp"

#include <iostream>

namespace utopia {

    /**
     * @brief      Base class for all linear multilevel solvers. \n
     *             Takes care of inializing multilevel hierarchy. \n
     *             Different levels are created by interpolation and restriction operators.\n
     *             Additionally, it provides stifness matrices on each level, created by using Galerkin assembly. \n
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector>
    class LinearMultiLevel : public MultiLevelBase<Matrix, Vector>, public IterativeSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        typedef utopia::Level<Matrix, Vector> Level;
        typedef utopia::Transfer<Matrix, Vector> Transfer;
        typedef utopia::MatrixTransfer<Matrix, Vector> MatrixTransfer;

    public:
        using MultiLevelBase<Matrix, Vector>::set_transfer_operators;

        LinearMultiLevel() : MultiLevelBase<Matrix, Vector>() {}

        void read(Input &in) override {
            IterativeSolver<Matrix, Vector>::read(in);
            MultiLevelBase<Matrix, Vector>::read(in);

            bool flg_masks;
            in.get("must_generate_masks", flg_masks);
            mask_.active(flg_masks);
        }

        LinearMultiLevel *clone() const override = 0;

        void print_usage(std::ostream &os) const override {
            IterativeSolver<Matrix, Vector>::print_usage(os);
            MultiLevelBase<Matrix, Vector>::print_usage(os);
            this->print_param_usage(
                os, "must_generate_masks", "bool", "Flag deciding if masks should be generated.", "-");
        }

        ~LinearMultiLevel() override = default;

        /**
         * @brief
         * Function initializes restriction transfer operators.
         * Operators need to be ordered FROM COARSE TO FINE.
         *
         * @param[in]  operators                The inteprolation operators.
         *
         */
        virtual bool set_transfer_operators(const std::vector<std::shared_ptr<Matrix>> &interpolation_operators) {
            if (this->n_levels() <= 0) {
                this->n_levels(interpolation_operators.size() + 1);
            } else if (this->n_levels() != static_cast<SizeType>(interpolation_operators.size()) + 1) {
                utopia_error("utopia::MultilevelBase:: number of levels and transfer operators do not match ... \n");
            }

            this->transfers_.clear();
            for (auto I = interpolation_operators.begin(); I != interpolation_operators.end(); ++I)
                this->transfers_.push_back(std::make_shared<MatrixTransfer>(*I));

            return true;
        }

        bool must_generate_masks() { return mask_.active(); }

        void must_generate_masks(const bool must_generate_masks) { mask_.active(must_generate_masks); }

        void add_level(Level &&level) { levels_.push_back(std::move(level)); }

        void add_level(const Level &level) { levels_.push_back(level); }

        static void fix_semidefinite_operator(Matrix &A) {
            Vector d(row_layout(A), 0.0);

            // TODO(Patrick) FIXME use atomic_store once available

            {
                auto d_view = view_device(d);

                A.read(UTOPIA_LAMBDA(const SizeType &i, const SizeType &j, const Scalar &val) {
                    if (i == j && device::abs(val) < 1e-12) {
                        d_view.atomic_add(i, 1.0);
                    }
                });
            }

            d.transform_values(UTOPIA_LAMBDA(const Scalar &v) {
                if (v >= 1.0) {
                    return 1.0;
                } else {
                    return 0.0;
                }
            });

            A += Matrix(diag(d));
        }

        /**
         * @brief
         *        The function creates corser level operators provided by assembling on differnet levels of MG hierarchy
         *		The first element of the vector is the coarset matrix and the last is the fienst
         * @param[in]  stifness matrix for finest level
         *
         */
        inline bool set_linear_operators(const std::vector<std::shared_ptr<const Matrix>> &A) {
            if (this->n_levels() <= 0) {
                this->n_levels(A.size());
            } else if (this->n_levels() != A.size()) {
                utopia_error("utopia::MultilevelBase:: number of levels and linear operators do not match ... \n");
            }

            levels_.clear();
            levels_.insert(levels_.begin(), A.begin(), A.end());
            return true;
        }

        inline const Level &level(const SizeType l) const { return levels_[l]; }

        void describe(std::ostream &os = std::cout) const override {
            if (levels_.empty()) {
                return;
            }

            SizeType i = 0;
            for (const auto &l : levels_) {
                const auto &A = l.A();
                auto s = size(A);

                os << "level: " << ++i << ", n_dofs: " << s.get(0) << std::endl;
            }
        }

        virtual void update_transfer(const SizeType level, const std::shared_ptr<Transfer> &t) {
            this->transfers_[level] = t;
        }

        virtual bool write(const Path &path) const {
            int l = 0;
            for (auto &level : levels_) {
                Path file_name = path / ("mat_l" + std::to_string(l++) + ".m");
                if (!level.write(file_name)) {
                    return false;
                }
            }

            return true;
        }

    protected:
        std::vector<Level> levels_; /*!< vector of level operators     */
        MultiLevelMask<Matrix, Vector> mask_;

        /**
         * @brief
         *        The function creates corser level operators by using Galerkin assembly.
         *
         *        $\f J_{i-1} = R * J_{i} * I  $\f
         *
         *        Resulting operators are assumed to go from fines = 0 to coarse = numlevels_
         *
         * @param[in]  stifness matrix for finest level
         *
         */
        virtual bool galerkin_assembly(const std::shared_ptr<const Matrix> &A) {
            levels_.clear();
            SizeType t_s = this->transfers_.size();
            if (t_s <= 0) std::cerr << "Provide interpolation operators first!  \n";

            levels_.push_back(Level(A));

            auto L = this->n_levels();

            mask_.generate_masks(*A, this->transfers_);

            for (SizeType i = 1; i < L; i++) {
                // J_{i-1} = R * J_{i} * I
                std::shared_ptr<Matrix> J_h = std::make_shared<Matrix>();
                this->transfer(t_s - i).restrict(levels_[i - 1].A(), *J_h);

                if (this->fix_semidefinite_operators()) {
                    fix_semidefinite_operator(*J_h);
                }

                levels_.push_back(Level(J_h));
            }

            std::reverse(std::begin(levels_), std::end(levels_));
            return true;
        }

        virtual void apply_mask(const SizeType l, Vector &v) const { mask_.apply(l, v); }
    };

}  // namespace utopia

#endif  // UTOPIA_LINEAR_MULTI_LEVELHPP
