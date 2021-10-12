#include "utopia_petsc_BDDOperator.hpp"

// Concrete includes
#include "utopia_DeviceView.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_petsc_Factorization.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class BDDOperator<Matrix, Vector>::Impl {
    public:
        std::shared_ptr<Matrix> A_GG_, A_GI_, A_II_, A_IG_;
        std::shared_ptr<Vector> secant_G_;
        std::shared_ptr<Factorization<Matrix, Vector>> A_II_inv_;
        std::shared_ptr<Vector> xL_, rhsL_, A_IG_x_, sol_I_;

        static void parallel_to_serial(const Vector &x_from, Vector &x_to) {
            {
                // FIXME (copying stuff just because of abstractions!)
                auto x_from_view = local_view_device(x_from);
                auto x_to_view = local_view_device(x_to);
                parallel_for(
                    local_range_device(x_to),
                    UTOPIA_LAMBDA(const SizeType i) { x_to_view.set(i, x_from_view.get(i)); });
            }
        }

        static void serial_to_parallel(const Vector &x_from, Vector &x_to) { parallel_to_serial(x_from, x_to); }

        void init(const std::shared_ptr<Matrix> &A_GG,
                  const std::shared_ptr<Matrix> &A_GI,
                  const std::shared_ptr<Matrix> &A_II,
                  const std::shared_ptr<Matrix> &A_IG) {
            A_GG_ = A_GG;
            A_GI_ = A_GI;
            A_II_ = A_II;
            A_IG_ = A_IG;

            // A_II_inv_ = std::make_shared<Factorization<Matrix, Vector>>("superlu", "lu");
            A_II_inv_ = std::make_shared<Factorization<Matrix, Vector>>("mumps", "cholesky");
            A_II_inv_->update(A_II_);
        }

        void init_rhs(const Vector &rhs_G, const Vector &rhs_I) {
            assert(A_II_inv_);

            // Compute rhs
            secant_G_ = std::make_shared<Vector>(layout(rhs_G));
            xL_ = std::make_shared<Vector>(row_layout(*A_GI_));
            rhsL_ = std::make_shared<Vector>(row_layout(*A_GI_));
            A_IG_x_ = std::make_shared<Vector>(row_layout(*A_II_));
            sol_I_ = std::make_shared<Vector>(row_layout(*A_II_));

            Vector inv_A_II_rhs_I(layout(rhs_I), 0);
            A_II_inv_->apply(rhs_I, inv_A_II_rhs_I);

            (*secant_G_) = rhs_G;

            Vector temp = (*A_GI_) * inv_A_II_rhs_I;

            assert(secant_G_->local_size() == temp.local_size());

            {
                auto sG_view = local_view_device(*secant_G_);
                auto temp_view = local_view_device(temp);

                parallel_for(
                    local_range_device(temp),
                    UTOPIA_LAMBDA(const SizeType i) { sG_view.set(i, sG_view.get(i) - temp_view.get(i)); });
            }

            // disp(rhs_G);
            // disp(*secant_G_);

            // for (SizeType r = 0; r < comm().size(); ++r) {
            //     comm().barrier();

            //     if (r == comm().rank()) {
            //         std::cout << "=========================================\n";
            //         std::cout << "Rank: " << r << "\n";
            //         std::cout << "=========================================\n";
            //         std::cout << "A_GI * A_II^-1 * rhs_I:\n";
            //         disp(temp);

            //         std::cout << "\n" << std::flush;
            //     }
            // }
        }

        bool apply(const Vector &x_G, Vector &rhs_G) const {
            parallel_to_serial(x_G, *xL_);

            *A_IG_x_ = (*A_IG_) * (*xL_);
            sol_I_->set(0);

            A_II_inv_->apply(*A_IG_x_, *sol_I_);

            (*rhsL_) = (*A_GI_) * (*sol_I_);

            if (empty(rhs_G)) {
                rhs_G.zeros(layout(x_G));
            }

            (*rhsL_) *= -1;

            serial_to_parallel(*rhsL_, rhs_G);
            rhs_G += (*A_GG_) * x_G;
            return true;
        }

        void finalize(const Vector &x_G, const Vector &rhs_I, Vector &x_I) {
            parallel_to_serial(x_G, *xL_);

            *A_IG_x_ = (*A_IG_) * (*xL_);
            *rhsL_ = rhs_I;
            *rhsL_ -= *A_IG_x_;

            if (empty(x_I)) {
                x_I.zeros(layout(*rhsL_));
            } else {
                x_I.set(0);
            }

            A_II_inv_->apply(*rhsL_, x_I);
        }
    };

    template <class Matrix, class Vector>
    BDDOperator<Matrix, Vector>::BDDOperator() : impl_(utopia::make_unique<Impl>()) {}

    template <class Matrix, class Vector>
    BDDOperator<Matrix, Vector>::~BDDOperator() = default;

    template <class Matrix, class Vector>
    bool BDDOperator<Matrix, Vector>::update(const std::shared_ptr<Matrix> &matrix) {
        return false;
    }

    template <class Matrix, class Vector>
    bool BDDOperator<Matrix, Vector>::apply(const Vector &x_G, Vector &rhs_G) const {
        return impl_->apply(x_G, rhs_G);
    }

    template <class Matrix, class Vector>
    Size BDDOperator<Matrix, Vector>::size() const {
        return impl_->A_GG_->size();
    }

    template <class Matrix, class Vector>
    Size BDDOperator<Matrix, Vector>::local_size() const {
        return impl_->A_GG_->local_size();
    }

    template <class Matrix, class Vector>
    typename Traits<Vector>::Communicator &BDDOperator<Matrix, Vector>::comm() {
        return impl_->A_GG_->comm();
    }

    template <class Matrix, class Vector>
    const typename Traits<Vector>::Communicator &BDDOperator<Matrix, Vector>::comm() const {
        return impl_->A_GG_->comm();
    }

    template class BDDOperator<PetscMatrix, PetscVector>;

}  // namespace utopia