#ifndef UTOPIA_LIBMESH_ASSEMBLER_HPP
#define UTOPIA_LIBMESH_ASSEMBLER_HPP

#include "utopia_petsc.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Vector.hpp"

#include "utopia_Assembler.hpp"

#include "utopia_kernels.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/reference_elem.h"

namespace utopia {

    bool libmesh_each(const LMFunctionSpace &space, std::function<bool(const libMesh::Elem &)> f);

    template <class Kernel>
    class LMLocalAssembler {};

    template <typename Scalar>
    class LMLocalAssembler<Poisson<Scalar, Scalar>> {
    public:
        using Model = Poisson<Scalar, Scalar>;
        using ElementMatrix = Traits<LMFunctionSpace>::ElementMatrix;
        using ElementVector = Traits<LMFunctionSpace>::ElementVector;

        void init_bilinear(const LMFunctionSpace &space) {
            fe_ = space.make_fe(model_.var);

            auto order = libMesh::Order(std::max(2 * (space.order(model_.var) - 1), 0));
            q_gauss_ = utopia::make_unique<libMesh::QGauss>(space.dim(), order);
            fe_->attach_quadrature_rule(q_gauss_.get());
            fe_->get_dphi();
            fe_->get_JxW();
        }

        void resize_vector(const LMFunctionSpace &space, ElementVector &el_vec) {}
        void resize_matrix(const LMFunctionSpace &space, ElementMatrix &el_mat) {}

        void reinit(const libMesh::Elem &el) { fe_->reinit(&el); }

        void init_linear(const LMFunctionSpace &space) {
            fe_ = space.make_fe(model_.var);

            auto order = libMesh::Order(std::max(2 * (space.order(model_.var)), 0));
            q_gauss_ = utopia::make_unique<libMesh::QGauss>(space.dim(), order);
            fe_->attach_quadrature_rule(q_gauss_.get());
            fe_->get_JxW();
            fe_->get_phi();
        }

        bool assemble_bilinear(ElementMatrix &el_mat) { return false; }

        bool assemble_linear(ElementVector &el_vec) { return false; }

        LMLocalAssembler(const Model &model) : model_(model) {}

    private:
        Model model_;
        std::unique_ptr<libMesh::FEBase> fe_;
        std::unique_ptr<libMesh::QGauss> q_gauss_;
    };

    template <>
    class Assembler<LMFunctionSpace> {
    public:
        using Traits = utopia::Traits<LMFunctionSpace>;
        using Scalar = Traits::Scalar;
        using SizeType = Traits::SizeType;
        using Vector = Traits::Vector;
        using Matrix = Traits::Matrix;
        using IndexArray = Traits::IndexArray;

        using ElementMatrix = Traits::ElementMatrix;
        using ElementVector = Traits::ElementVector;

        template <class Model>
        bool assemble(Model &model, Matrix &mat, Vector &rhs) {
            return assemble(model, mat) && assemble(model, rhs);
        }

        template <class Model>
        bool assemble(Model &model, Matrix &mat) {
            if (mat.empty()) {
                space_.create_matrix(mat);
            } else {
                mat *= 0.0;
            }

            LMLocalAssembler<Model> assembler(model);
            assembler.init_bilinear(space_);

            ElementMatrix el_mat;
            assembler.resize_matrix(space_, el_mat);

            std::vector<libMesh::dof_id_type> dofs;
            auto &dof_map = space_.dof_map();

            Write<Matrix> w_m(mat, utopia::GLOBAL_ADD);
            return libmesh_each(space_, [&](const libMesh::Elem &e) -> bool {
                assembler.reinit(e);
                el_mat.set(0.0);
                if (assembler.assemble_bilinear(el_mat)) {
                    dof_map.dof_indices(&e, dofs);
                    el_mat.read([&](const SizeType &i, const SizeType &j, const Scalar &val) {
                        mat.c_add(dofs[i], dofs[j], val);
                        return true;
                    });
                } else {
                    return false;
                }
            });
        }

        template <class Model>
        bool assemble(Model &model, Vector &vec) {
            if (vec.empty()) {
                space_.create_vector(vec);
            } else {
                vec.set(0.0);
            }

            LMLocalAssembler<Model> assembler(model);
            assembler.init_linear(space_);

            ElementVector el_vec;
            assembler.resize_vector(space_, el_vec);

            std::vector<libMesh::dof_id_type> dofs;
            auto &dof_map = space_.dof_map();

            Write<Vector> w_m(vec, utopia::GLOBAL_ADD);
            return libmesh_each(space_, [&](const libMesh::Elem &e) -> bool {
                assembler.reinit(e);
                el_vec.set(0.0);
                if (assembler.assemble_linear(el_vec)) {
                    dof_map.dof_indices(&e, dofs);

                    const SizeType n = el_vec.size();

                    assert(n == SizeType(dofs.size()));
                    for (SizeType i = 0; i < n; ++i) {
                        vec.c_add(dofs[i], el_vec.get(i));
                    }

                } else {
                    return false;
                }
            });
        }

        Assembler(const LMFunctionSpace &space) : space_(space) {}

    private:
        const LMFunctionSpace &space_;
    };
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_ASSEMBLER_HPP
