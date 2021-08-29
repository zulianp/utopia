#ifndef UTOPIA_PETSC_FUNCTION_SPACE_HPP
#define UTOPIA_PETSC_FUNCTION_SPACE_HPP

#include "utopia_fe_Core.hpp"

#include "utopia_DirichletBoundary.hpp"
#include "utopia_FunctionSpaceBase.hpp"

#include "utopia_petsc_StructuredGrid.hpp"

namespace utopia {

    template <>
    class Traits<petsc::FunctionSpace> : public Traits<petsc::StructuredGrid> {
    public:
        using Mesh = petsc::StructuredGrid;
    };

    namespace petsc {

        class FunctionSpace : public FunctionSpaceBase<petsc::StructuredGrid>, public Traits<FunctionSpace> {
        public:
            using Super = utopia::FunctionSpaceBase<petsc::StructuredGrid>;
            using Traits = utopia::Traits<FunctionSpace>;
            using Vector = Traits::Vector;
            using Matrix = Traits::Matrix;
            using Scalar = Traits::Scalar;
            using SizeType = Traits::SizeType;
            using Communicator = Traits::Communicator;
            using Mesh = Traits::Mesh;

            FunctionSpace(const Communicator &comm = Communicator::get_default());
            FunctionSpace(const std::shared_ptr<Mesh> &mesh);
            ~FunctionSpace();

            void read(Input &in) override;

            //////////////////////////////////////////

            void init(const std::shared_ptr<Mesh> &mesh) override;
            std::shared_ptr<Mesh> mesh_ptr() const override;
            const Mesh &mesh() const override;
            Mesh &mesh() override;
            int n_var() const override;

            // Below could be made into a cross-backend interface (given that the same algebra is used)
            bool write(const Path &path, const Vector &x) override;

            const Communicator &comm() const override;

            SizeType n_dofs() const override;
            SizeType n_local_dofs() const override;

            void create_vector(Vector &v) const override;
            void create_matrix(Matrix &m) const override;

            void apply_constraints(Matrix &m, const Scalar diag_value = 1.0) const override;
            void apply_constraints(Vector &v) const override;
            void apply_constraints(Matrix &m, Vector &v) const override;
            void apply_zero_constraints(Vector &vec) const override;

            void add_dirichlet_boundary_condition(const std::string &name,
                                                  const Scalar &value,
                                                  const int component = 0) override;

            bool empty() const override;

            void displace(const Vector &displacement) override;

            const std::string &name() const override;
            void initialize() override;

            void create_field(Field<FunctionSpace> &field);
            void global_to_local(const Vector &global, Vector &local) const;
            void local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const;

            // bool is_non_conforming() const { return false; }
            // std::shared_ptr<Matrix> constraint_matrix() const { return nullptr; }

            // int n_components() const;

            // template <class DofIndex>
            // void dofs(const SizeType &idx, DofIndex &dofs) const {
            //     DofMapping::dofs(*mesh_, subspace_id_, idx, dofs);
            // }

            // template <class DofIndex>
            // void dofs_local(const SizeType &idx, DofIndex &dofs) const {
            //     DofMapping::dofs_local(*mesh_, subspace_id_, idx, dofs);
            // }

            // template <class ElementMatrix, class MatView>
            // void add_matrix(const Elem &e, const ElementMatrix &el_mat, MatView &mat) const {
            //     DofMapping::add_matrix(*mesh_, subspace_id_, e, el_mat, mat);
            // }

            // template <class ElementVector, class VecView>
            // void add_vector(const Elem &e, const ElementVector &el_vec, VecView &vec) const {
            //     DofMapping::add_vector(*mesh_, subspace_id_, e, el_vec, vec);
            // }

            // template <class ElementVector, class VecView>
            // void set_vector(const Elem &e, const ElementVector &el_vec, VecView &vec) const {
            //     DofMapping::set_vector(*mesh_, subspace_id_, e, el_vec, vec);
            // }

            // template <class VectorView, class Values>
            // void local_coefficients(const Elem &e, const VectorView &vec, Values &values) const {
            //     DofMapping::local_coefficients(*mesh_, subspace_id_, e, vec, values);
            // }

            // template <class VectorView, class Values>
            // void local_coefficients(const Elem &e, const VectorView &vec, const SizeType &var, Values &values) const
            // {
            //     DofMapping::local_coefficients_for_var(*mesh_, e, vec, subspace_id_ + var, values);
            // }

            // bool on_boundary(const SizeType &elem_idx) const { return mesh_->on_boundary(elem_idx); }

            // template <class... Args>
            // void emplace_dirichlet_condition(Args &&... args) {
            //     dirichlet_bcs_.push_back(utopia::make_unique<DirichletBC>(*this, std::forward<Args>(args)...));
            // }

            // void apply_constraints(Matrix &mat, Vector &vec) const {
            //     for (const auto &bc : dirichlet_bcs_) {
            //         bc->apply(mat, vec);
            //     }
            // }

            // void apply_constraints(Matrix &mat) const {
            //     for (const auto &bc : dirichlet_bcs_) {
            //         bc->apply(mat);
            //     }
            // }

            // void apply_constraints(Vector &vec) const {
            //     for (const auto &bc : dirichlet_bcs_) {
            //         bc->apply(vec);
            //     }
            // }

            // void copy_at_constrained_dofs(const Vector &in, Vector &vec) const {
            //     for (const auto &bc : dirichlet_bcs_) {
            //         bc->copy(in, vec);
            //     }
            // }

            // void apply_zero_constraints(Vector &vec) const {
            //     for (const auto &bc : dirichlet_bcs_) {
            //         bc->apply_zero(vec);
            //     }
            // }

            // void build_constraints_markers(Vector &vec) const {
            //     for (const auto &bc : dirichlet_bcs_) {
            //         bc->apply_val(vec, 1.0);
            //     }
            // }

            // //////////////////////////////////////////

            // void create_local_vector(Vector &vec) const { mesh_->create_local_vector(vec); }

            // void global_to_local(const Vector &global, Vector &local) const { mesh_->global_to_local(global, local);
            // }

            // const Mesh &mesh() const { return *mesh_; }

            // Mesh &mesh() { return *mesh_; }

            // PetscCommunicator &comm() { return mesh_->comm(); }

            // const PetscCommunicator &comm() const { return mesh_->comm(); }

            // Range dof_range() const { return mesh_->dof_range(); }

            // SizeType n_dofs() const { return mesh_->n_nodes() * NComponents; }

            // UTOPIA_INLINE_FUNCTION constexpr SizeType n_local_dofs() const { return mesh_->n_local_dofs(); }

            // void set_mesh(const std::shared_ptr<Mesh> &mesh) {
            //     mesh_ = mesh;
            //     elements_ = mesh_->elements_ptr();
            // }

            // void set_subspace_id(const SizeType &i) { subspace_id_ = i; }

            // void set_dirichlet_conditions(const std::vector<std::shared_ptr<DirichletBC>> &conds) {
            //     dirichlet_bcs_ = conds;
            // }

            // void reset_bc() { dirichlet_bcs_.clear(); }

            // inline Range element_range() const {
            //     assert(elements_);
            //     return elements_->range();
            // }

            // inline void elem(const SizeType &idx, Elem &e) const {
            //     MakeElem<FunctionSpace, Elem>::apply(*this, idx, e);
            // }

            // const ViewDevice &view_device() const { return *this; }

            // FunctionSpace<Mesh, 1, UniVarElem_> subspace(const SizeType &i) const {
            //     FunctionSpace<Mesh, 1, UniVarElem_> space(mesh_, subspace_id_ + i);
            //     // space.set_dirichlet_conditions(dirichlet_bcs_);
            //     assert(i < NComponents);

            //     // utopia::out() <<i << " " << subspace_id_ << " " << mesh_->n_components() << std::endl;

            //     assert(i + subspace_id_ < mesh_->n_components());
            //     return space;
            // }

            // template <int NVars>
            // void subspace(const SizeType &i, FunctionSpace<Mesh, NVars, UniVarElem_> &space) const {
            //     space.set_mesh(mesh_);
            //     space.set_subspace_id(subspace_id_ + i);
            // }

            // template <int NVars>
            // FunctionSpace<Mesh, NVars, UniVarElem_> vector_subspace(const SizeType &i) const {
            //     FunctionSpace<Mesh, NVars, UniVarElem_> space(mesh_, subspace_id_ + i);
            //     // space.set_dirichlet_conditions(dirichlet_bcs_);
            //     assert(i + NVars < NComponents);
            //     assert(subspace_id_ + i < mesh_->n_components());
            //     return space;
            // }

            // std::unique_ptr<FunctionSpace> uniform_refine() const {
            //     auto fine_space = utopia::make_unique<FunctionSpace>(mesh_->uniform_refine(), subspace_id_);

            //     const std::size_t n = dirichlet_bcs_.size();
            //     fine_space->dirichlet_bcs_.resize(n);

            //     for (std::size_t i = 0; i < n; ++i) {
            //         fine_space->dirichlet_bcs_[i] = std::make_shared<DirichletBC>(*fine_space);
            //         fine_space->dirichlet_bcs_[i]->init_from(*dirichlet_bcs_[i]);
            //     }

            //     return fine_space;
            // }

            // template <class F>
            // void sample(Vector &v, F f, const int /*c*/ = 0) {
            //     auto r = v.range();
            //     // auto n = r.extent() * NComponents;
            //     assert(!v.empty());

            //     // Write<Vector> w(v, utopia::AUTO);

            //     // Point p;
            //     // for(auto i = r.begin(); i < r.end(); ++i) {
            //     //     this->mesh().node(i/NComponents, p);
            //     //     v.set(i, f(p));
            //     // }

            //     // FIXME
            //     {
            //         auto space_view = view_device();
            //         auto v_view = utopia::view_device(v);

            //         Device::parallel_for(
            //             this->element_range(), UTOPIA_LAMBDA(const SizeType &i) {
            //                 Elem e;
            //                 space_view.elem(i, e);

            //                 NodeIndex nodes;
            //                 space_view.mesh().nodes(i, nodes);
            //                 const SizeType n_nodes = nodes.size();

            //                 Point p;
            //                 for (SizeType i = 0; i < n_nodes; ++i) {
            //                     auto idx = nodes[i] * mesh_->n_components() + subspace_id_;
            //                     e.node(i, p);

            //                     if (r.inside(idx)) {
            //                         v_view.set(idx, f(p));
            //                     }
            //                 }
            //             });
            //     }
            // }

            // inline void create_interpolation(const FunctionSpace &target, Matrix &I) const {
            //     mesh_->create_interpolation(target.mesh(), I);
            // }

            // void read(Input &in) override {
            //     PetscCommunicator comm;
            //     if (!mesh_) {
            //         allocate_mesh(comm);
            //     }

            //     mesh_->read(in);
            //     elements_ = mesh_->elements_ptr();
            // }

            // bool write(const Path &path, const Vector &x) const { return mesh_->write(path, x); }

            // inline bool empty() const { return static_cast<bool>(mesh_); }

            // inline SizeType component(const SizeType &idx) const {
            //     const SizeType nc = mesh_->n_components();
            //     return nc == 1 ? 0 : idx % nc;
            // }

            // FunctionSpace() : subspace_id_(0) {}

            // FunctionSpace(const PetscCommunicator &comm, const SizeType subspace_id = 0) : subspace_id_(subspace_id)
            // {
            //     allocate_mesh(comm);
            // }

            // FunctionSpace(const std::shared_ptr<Mesh> &mesh, const SizeType subspace_id = 0)
            //     : mesh_(mesh), subspace_id_(subspace_id) {
            //     elements_ = mesh_->elements_ptr();
            // }

            // void describe() const {
            //     assert(mesh_);
            //     mesh_->describe();
            // }

        private:
            class Impl;
            std::unique_ptr<Impl> impl_;
            // std::shared_ptr<Mesh> mesh_;
            // std::shared_ptr<typename Mesh::Elements> elements_;
            // std::vector<std::shared_ptr<DirichletBC>> dirichlet_bcs_;
            // SizeType subspace_id_;

            // void allocate_mesh(const PetscCommunicator &comm) {
            //     mesh_ =
            //         std::make_shared<Mesh>(comm, is_simplex<UniVarElem_>::value ? DMDA_ELEMENT_P1 : DMDA_ELEMENT_Q1);
            //     mesh_->set_n_components(NComponents_);

            //     assert(subspace_id_ + NComponents_ <= mesh_->n_components());
            // }
        };

    }  // namespace petsc
}  // namespace utopia

#endif  // UTOPIA_PETSC_FUNCTION_SPACE_HPP
