#ifndef UTOPIA_PETSC_DM_PLEX_FUNCTION_SPACE_HPP
#define UTOPIA_PETSC_DM_PLEX_FUNCTION_SPACE_HPP

#include "utopia_FunctionSpace.hpp"
#include "utopia_MultiVariateElement.hpp"
#include "utopia_petsc_DMPlex.hpp"

// TOBEREMOVED
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_MakeElem.hpp"

namespace utopia {

    template <class Point_, class IntArray_, int NComponents_, class UniVarElem_>
    class FunctionSpace<PetscDMPlex<Point_, IntArray_>, NComponents_, UniVarElem_> : public Configurable {
    public:
        // concrete types
        using Vector = utopia::PetscVector;
        using Matrix = utopia::PetscMatrix;
        using Comm = utopia::PetscCommunicator;

        static const int Dim = UniVarElem_::Dim;

        // from template arg list
        using Point = Point_;
        using Mesh = utopia::PetscDMPlex<Point_, IntArray_>;
        static const int NComponents = NComponents_;
        using Shape = UniVarElem_;
        using Elem = MultiVariateElem<UniVarElem_, NComponents_>;
        // static const int Dim = Mesh::StaticDim;

        using NodeIndex = typename Mesh::NodeIndex;

        using Scalar = typename Mesh::Scalar;
        using SizeType = typename Mesh::SizeType;

        using ViewDevice = FunctionSpace;
        using Device = typename Mesh::Device;

        using DirichletBC = utopia::DirichletBoundaryCondition<FunctionSpace>;
        using MemType = utopia::Varying<>;

        template <int NSubVars>
        using Subspace = FunctionSpace<Mesh, NSubVars, UniVarElem_>;

        //////////////////////////////////////////

        inline static constexpr int n_components() { return NComponents; }

        //////////////////////////////////////////

        void create_matrix(PetscMatrix &mat) const { mesh_->create_matrix(mat); }

        void create_vector(PetscVector &vec) const { mesh_->create_vector(vec); }

        void create_local_vector(PetscVector &vec) const { mesh_->create_local_vector(vec); }

        template <class ElementMatrix, class MatView>
        void add_matrix(const Elem &e, const ElementMatrix &el_mat, MatView &mat) const {
            ArrayView<SizeType, Elem::NFunctions> dofs;

            // FIXME
            mesh_->fields_local(e.idx(), 0, dofs);

            // disp(dofs);

            if (NComponents > 1) {
                for (SizeType c = NComponents_ - 1; c >= 0; --c) {
                    for (SizeType i = 0; i < Elem::NNodes; ++i) {
                        dofs[c * Elem::NNodes + i] = dofs[i] + c;
                    }
                }
            }

            // https://www.mcs.anl.gov/petsc/petsc-current/src/tao/tutorials/ex3.c.html

            //   for (p = 0, k = 0; p < n; ++p) {
            // 267:     PetscSectionGetDof(section, points[p], &dof);
            // 268:     PetscSectionGetOffset(section, points[p], &off);
            // 269:     for (d = 0; d < dof; ++d) user->bc_indices[k++] = off+d;
            // 270:   }

            // disp(dofs);
            mat.atomic_add_matrix(dofs, dofs, &el_mat(0, 0));
        }

        inline SizeType component(const SizeType &idx) const { return mesh_->component(idx); }

        template <class... Args>
        void emplace_dirichlet_condition(Args &&... args) {
            dirichlet_bcs_.push_back(utopia::make_unique<DirichletBC>(*this, std::forward<Args>(args)...));
        }

        void apply_constraints(PetscMatrix &mat, PetscVector &vec) const {
            for (const auto &bc : dirichlet_bcs_) {
                bc->apply(mat, vec);
            }
        }

        void apply_constraints(PetscMatrix &mat) const {
            for (const auto &bc : dirichlet_bcs_) {
                bc->apply(mat);
            }
        }

        void apply_constraints(PetscVector &vec) const {
            for (const auto &bc : dirichlet_bcs_) {
                bc->apply(vec);
            }
        }

        void copy_at_constrained_dofs(const PetscVector &in, PetscVector &vec) const {
            for (const auto &bc : dirichlet_bcs_) {
                bc->copy(in, vec);
            }
        }

        void apply_zero_constraints(PetscVector &vec) const {
            for (const auto &bc : dirichlet_bcs_) {
                bc->apply_zero(vec);
            }
        }

        static DeviceView<PetscMatrix, 2> assembly_view_device(PetscMatrix &mat) {
            return DeviceView<PetscMatrix, 2>(mat, utopia::GLOBAL_ADD);
        }

        static DeviceView<PetscVector, 1> assembly_view_device(PetscVector &vec) {
            return DeviceView<PetscVector, 1>(vec, utopia::GLOBAL_ADD);
        }

        static DeviceView<const PetscVector, 1> assembly_view_device(const PetscVector &vec) {
            return DeviceView<const PetscVector, 1>(vec);
        }

        void global_to_local(const PetscVector &global, PetscVector &local) const {
            mesh_->global_to_local(global, local);
        }

        const Mesh &mesh() const { return *mesh_; }

        Mesh &mesh() { return *mesh_; }

        PetscCommunicator &comm() { return mesh_->comm(); }

        const PetscCommunicator &comm() const { return mesh_->comm(); }

        void set_mesh(const std::shared_ptr<Mesh> &mesh) {
            mesh_ = mesh;
            // elements_ = mesh_->elements_ptr();
        }

        void set_subspace_id(const SizeType &i) { subspace_id_ = i; }

        const ViewDevice &view_device() const { return *this; }

        FunctionSpace<Mesh, 1, UniVarElem_> subspace(const SizeType &i) const {
            FunctionSpace<Mesh, 1, UniVarElem_> space(mesh_, subspace_id_ + i);
            // space.set_dirichlet_conditions(dirichlet_bcs_);
            assert(i < NComponents);

            // std::cout << i << " " << subspace_id_ << " " << mesh_->n_components() << std::endl;

            // assert(i + subspace_id_ < mesh_->n_components());
            return space;
        }

        template <int NVars>
        void subspace(const SizeType &i, FunctionSpace<Mesh, NVars, UniVarElem_> &space) const {
            space.set_mesh(mesh_);
            space.set_subspace_id(subspace_id_ + i);
        }

        template <int NVars>
        FunctionSpace<Mesh, NVars, UniVarElem_> vector_subspace(const SizeType &i) const {
            FunctionSpace<Mesh, NVars, UniVarElem_> space(mesh_, subspace_id_ + i);
            // space.set_dirichlet_conditions(dirichlet_bcs_);
            assert(i + NVars < NComponents);
            // assert(subspace_id_ + i < mesh_->n_components());
            return space;
        }

        std::unique_ptr<FunctionSpace> uniform_refine() const {
            auto fine_space = utopia::make_unique<FunctionSpace>(mesh_->uniform_refine(), subspace_id_);

            // const std::size_t n = dirichlet_bcs_.size();
            // fine_space->dirichlet_bcs_.resize(n);

            // for(std::size_t i = 0; i < n; ++i) {
            //     fine_space->dirichlet_bcs_[i] = std::make_shared<DirichletBC>(*fine_space);
            //     fine_space->dirichlet_bcs_[i]->init_from(*dirichlet_bcs_[i]);
            // }

            return std::move(fine_space);
        }

        inline void create_interpolation(const FunctionSpace &target, PetscMatrix &I) const {
            mesh_->create_interpolation(target.mesh(), I);
        }

        void read(Input &in) override {
            PetscCommunicator comm;
            if (!mesh_) {
                allocate_mesh(comm);
            }

            // FIXME
            mesh_->simplex(true);

            mesh_->read(in);
            mesh_->set_num_fields(1);

            // FIXME Important to specifiy for all num_fields (What does this field do exactly?)
            SizeType num_comp[1] = {NComponents_};

            // int order = 1;
            // in.get("order", order);

            // Important to specifiy all mesh dim 0, 1, 2, 3
            SizeType num_dofs[4] = {NComponents_, 0, 0, 0};

            if (Elem::Order == 2) {
                num_dofs[1] = NComponents_;
            }

            mesh_->create_section(num_comp, num_dofs);

            // dmplex.set_field_name(0, "u");

            mesh_->set_up();
        }

        inline void elem(const SizeType &idx, Elem &mve) const {
            // FIXME (make it for all elem types)
            auto &e = mve.univar_elem();
            e.idx(idx);

            // DMPolytopeType ct;
            // DMPlexGetCellType(mesh_->raw_type(), idx, &ct);

            // std::cout << " DMPolytopeType: " << ct << std::endl;

            DMPlexComputeCellGeometryAffineFEM(mesh_->raw_type(),
                                               idx,
                                               &e.translation()[0],
                                               &e.jacobian().raw_type()[0],
                                               &e.jacobian_inverse().raw_type()[0],
                                               &e.measure());

            // reference size is changed from 2 to 1
            e.jacobian() *= 2.0;
            e.jacobian_inverse() *= 0.5;
            e.measure() *= 2.0;
        }

        // inline SizeType component(const SizeType &idx) const {
        //     const SizeType nc = mesh_->n_components();
        //     return nc == 1 ? 0 : idx % nc;
        // }

        inline Range element_range() const { return mesh_->element_range(); }

        bool write(const Path &path, const PetscVector &x) const { return mesh_->write(path, x); }

        inline bool empty() const { return static_cast<bool>(mesh_); }

        FunctionSpace() : subspace_id_(0) {}

        FunctionSpace(const PetscCommunicator &, const SizeType subspace_id = 0) : subspace_id_(subspace_id) {
            // allocate_mesh(comm);
        }

        FunctionSpace(const std::shared_ptr<Mesh> &mesh, const SizeType subspace_id = 0)
            : mesh_(mesh), subspace_id_(subspace_id) {
            // elements_ = mesh_->elements_ptr();
        }

        void describe() const {
            assert(mesh_);
            mesh_->describe();
        }

        template <class Quadrature>
        ShapeFunction<FunctionSpace, Quadrature> shape(const Quadrature &q) {
            return ShapeFunction<FunctionSpace, Quadrature>(*this, q);
        }

        template <class Quadrature>
        PhysicalGradient<FunctionSpace, Quadrature> shape_grad(const Quadrature &q) {
            return PhysicalGradient<FunctionSpace, Quadrature>(*this, q);
        }

        template <class Quadrature>
        PhysicalPoint<FunctionSpace, Quadrature> points(const Quadrature &q) {
            return PhysicalPoint<FunctionSpace, Quadrature>(*this, q);
        }

        template <class Quadrature>
        ShapeFunction<typename ViewDevice::Elem::Side, typename Quadrature::ViewDevice> side_shape_device(
            const Quadrature &q) {
            return ShapeFunction<typename ViewDevice::Elem::Side, typename Quadrature::ViewDevice>(q.view_device());
        }

        template <class Quadrature>
        PhysicalPoint<typename ViewDevice::Elem::Side, typename Quadrature::ViewDevice> side_points_device(
            const Quadrature &q) {
            return PhysicalPoint<typename ViewDevice::Elem::Side, typename Quadrature::ViewDevice>(q.view_device());
        }

        template <class Quadrature>
        Differential<typename ViewDevice::Elem::Side, typename Quadrature::ViewDevice> side_differential_device(
            const Quadrature &q) {
            return Differential<typename ViewDevice::Elem::Side, typename Quadrature::ViewDevice>(q.view_device());
        }

        template <class Quadrature>
        Differential<FunctionSpace, Quadrature> differential(const Quadrature &q) {
            return Differential<FunctionSpace, Quadrature>(*this, q);
        }

    private:
        std::shared_ptr<Mesh> mesh_;
        // std::shared_ptr<typename Mesh::Elements> elements_;
        std::vector<std::shared_ptr<DirichletBC>> dirichlet_bcs_;
        SizeType subspace_id_;

        void allocate_mesh(const Comm &comm) { mesh_ = std::make_shared<Mesh>(comm); }
    };  // namespace utopia

    template <class Point_, class IntArray_, int NComponents_, class UniVarElem_>
    using PetscDMPlexFunctionSpace = FunctionSpace<PetscDMPlex<Point_, IntArray_>, NComponents_, UniVarElem_>;

    template <class Point_, class IntArray_, int NComponents_, class UniVarElem_>
    const int FunctionSpace<PetscDMPlex<Point_, IntArray_>, NComponents_, UniVarElem_>::NComponents;

}  // namespace utopia

#endif  // UTOPIA_PETSC_DM_PLEX_FUNCTION_SPACE_HPP
