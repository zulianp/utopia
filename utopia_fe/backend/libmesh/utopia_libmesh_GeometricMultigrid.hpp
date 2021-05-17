#ifndef UTOPIA_LIBMESH_GEOMETRIC_MULTIGRID_HPP
#define UTOPIA_LIBMESH_GEOMETRIC_MULTIGRID_HPP

#include "utopia_Input.hpp"
#include "utopia_Multigrid.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_libmesh_FunctionSpace_new.hpp"

namespace utopia {

    namespace libmesh {
        class GeometricMultigrid final : public IterativeSolver<Traits<utopia::libmesh::FunctionSpace>::Matrix,
                                                                Traits<utopia::libmesh::FunctionSpace>::Vector> {
        public:
            using Matrix = Traits<utopia::libmesh::FunctionSpace>::Matrix;
            using Vector = Traits<utopia::libmesh::FunctionSpace>::Vector;
            using FunctionSpace = utopia::libmesh::FunctionSpace;

            bool init(const std::shared_ptr<FunctionSpace> &space, const SizeType n_levels);
            void read(Input &in) override { algebraic_mg_->read(in); }

            void update(const std::shared_ptr<const Matrix> &op) override { algebraic_mg_->update(op); }

            bool apply(const Vector &rhs, Vector &sol) override { return algebraic_mg_->apply(rhs, sol); }

            GeometricMultigrid *clone() const override {
                auto gmg = utopia::make_unique<GeometricMultigrid>(*this);
                return gmg.release();
            }

            inline void verbose(const bool &val) override { algebraic_mg_->verbose(val); }

            inline void max_it(const SizeType &it) override { algebraic_mg_->max_it(it); }

            inline SizeType max_it() const override { return algebraic_mg_->max_it(); }

            inline void atol(const double &tol) override { algebraic_mg_->atol(tol); }

        private:
            std::shared_ptr<LinearMultiLevel<Matrix, Vector>> algebraic_mg_;
        };

    }  // namespace libmesh
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_GEOMETRIC_MULTIGRID_HPP
