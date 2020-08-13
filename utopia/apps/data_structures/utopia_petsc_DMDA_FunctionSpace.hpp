#ifndef UTOPIA_PETSC_DMDA_FUNCTION_SPACE_HPP
#define UTOPIA_PETSC_DMDA_FUNCTION_SPACE_HPP

#include "utopia_FunctionSpace.hpp"
#include "utopia_MultiVariateElement.hpp"
#include "utopia_petsc_DMDA.hpp"
#include "utopia_petsc_DMDA_DofMap.hpp"

// TOBEREMOVED
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_MakeElem.hpp"
// #include "utopia_petsc_DMDA_FunctionSpace.hpp"

namespace utopia {

    template <class Point_, class IntArray_, int NComponents_, class UniVarElem_>
    class FunctionSpace<PetscDMDA<Point_, IntArray_>, NComponents_, UniVarElem_> : public Configurable {
    public:
        // concrete types
        using Vector = utopia::PetscVector;
        using Matrix = utopia::PetscMatrix;
        using Comm = utopia::PetscCommunicator;

        // from template arg list
        using Point = Point_;
        using Mesh = utopia::PetscDMDA<Point_, IntArray_>;
        using Shape = UniVarElem_;
        using Elem = MultiVariateElem<UniVarElem_, NComponents_>;
        static const int NComponents = NComponents_;
        static const int Dim = Mesh::StaticDim;

        using NodeIndex = typename Mesh::NodeIndex;

        using MemType = typename Elem::MemType;
        using Scalar = typename Mesh::Scalar;
        using SizeType = typename Mesh::SizeType;

        using ViewDevice = FunctionSpace;
        using Device = typename Mesh::Device;

        using DirichletBC = utopia::DirichletBoundaryCondition<FunctionSpace>;
        using DofMapping = utopia::DofMapping<Mesh, UniVarElem_, NComponents_>;
        static const int NDofs = DofMapping::NDofs;

        template <int NSubVars>
        using Subspace = FunctionSpace<Mesh, NSubVars, UniVarElem_>;

        //////////////////////////////////////////

        inline static constexpr int n_components() { return NComponents; }

        template <class DofIndex>
        void dofs(const SizeType &idx, DofIndex &dofs) const {
            DofMapping::dofs(*mesh_, subspace_id_, idx, dofs);
        }

        template <class DofIndex>
        void dofs_local(const SizeType &idx, DofIndex &dofs) const {
            DofMapping::dofs_local(*mesh_, subspace_id_, idx, dofs);
        }

        template <class ElementMatrix, class MatView>
        void add_matrix(const Elem &e, const ElementMatrix &el_mat, MatView &mat) const {
            DofMapping::add_matrix(*mesh_, subspace_id_, e, el_mat, mat);
        }

        template <class ElementVector, class VecView>
        void add_vector(const Elem &e, const ElementVector &el_vec, VecView &vec) const {
            DofMapping::add_vector(*mesh_, subspace_id_, e, el_vec, vec);
        }

        template <class ElementVector, class VecView>
        void set_vector(const Elem &e, const ElementVector &el_vec, VecView &vec) const {
            DofMapping::set_vector(*mesh_, subspace_id_, e, el_vec, vec);
        }

        template <class VectorView, class Values>
        void local_coefficients(const Elem &e, const VectorView &vec, Values &values) const {
            DofMapping::local_coefficients(*mesh_, subspace_id_, e, vec, values);
        }

        template <class VectorView, class Values>
        void local_coefficients(const Elem &e, const VectorView &vec, const SizeType &var, Values &values) const {
            DofMapping::local_coefficients_for_var(*mesh_, e, vec, subspace_id_ + var, values);
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

        bool on_boundary(const SizeType &elem_idx) const { return mesh_->on_boundary(elem_idx); }

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

        void build_constraints_markers(PetscVector &vec) const {
            for (const auto &bc : dirichlet_bcs_) {
                bc->apply_val(vec, 1.0);
            }
        }

        //////////////////////////////////////////

        void create_matrix(PetscMatrix &mat) const { mesh_->create_matrix(mat); }

        void create_vector(PetscVector &vec) const { mesh_->create_vector(vec); }

        void create_local_vector(PetscVector &vec) const { mesh_->create_local_vector(vec); }

        void global_to_local(const PetscVector &global, PetscVector &local) const {
            mesh_->global_to_local(global, local);
        }

        const Mesh &mesh() const { return *mesh_; }

        Mesh &mesh() { return *mesh_; }

        PetscCommunicator &comm() { return mesh_->comm(); }

        const PetscCommunicator &comm() const { return mesh_->comm(); }

        Range dof_range() const { return mesh_->dof_range(); }

        SizeType n_dofs() const { return mesh_->n_nodes() * NComponents; }

        UTOPIA_INLINE_FUNCTION constexpr SizeType n_local_dofs() const { return mesh_->n_local_dofs(); }

        void set_mesh(const std::shared_ptr<Mesh> &mesh) {
            mesh_ = mesh;
            elements_ = mesh_->elements_ptr();
        }

        void set_subspace_id(const SizeType &i) { subspace_id_ = i; }

        void set_dirichlet_conditions(const std::vector<std::shared_ptr<DirichletBC>> &conds) {
            dirichlet_bcs_ = conds;
        }

        void reset_bc() { dirichlet_bcs_.clear(); }

        inline Range element_range() const {
            assert(elements_);
            return elements_->range();
        }

        inline void elem(const SizeType &idx, Elem &e) const { MakeElem<FunctionSpace, Elem>::apply(*this, idx, e); }

        const ViewDevice &view_device() const { return *this; }

        FunctionSpace<Mesh, 1, UniVarElem_> subspace(const SizeType &i) const {
            FunctionSpace<Mesh, 1, UniVarElem_> space(mesh_, subspace_id_ + i);
            // space.set_dirichlet_conditions(dirichlet_bcs_);
            assert(i < NComponents);

            // utopia::out() <<i << " " << subspace_id_ << " " << mesh_->n_components() << std::endl;

            assert(i + subspace_id_ < mesh_->n_components());
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
            assert(subspace_id_ + i < mesh_->n_components());
            return space;
        }

        std::unique_ptr<FunctionSpace> uniform_refine() const {
            auto fine_space = utopia::make_unique<FunctionSpace>(mesh_->uniform_refine(), subspace_id_);

            const std::size_t n = dirichlet_bcs_.size();
            fine_space->dirichlet_bcs_.resize(n);

            for (std::size_t i = 0; i < n; ++i) {
                fine_space->dirichlet_bcs_[i] = std::make_shared<DirichletBC>(*fine_space);
                fine_space->dirichlet_bcs_[i]->init_from(*dirichlet_bcs_[i]);
            }

            return fine_space;
        }

        template <class F>
        void sample(Vector &v, F f, const int /*c*/ = 0) {
            auto r = v.range();
            // auto n = r.extent() * NComponents;
            assert(!v.empty());

            // Write<Vector> w(v, utopia::AUTO);

            // Point p;
            // for(auto i = r.begin(); i < r.end(); ++i) {
            //     this->mesh().node(i/NComponents, p);
            //     v.set(i, f(p));
            // }

            // FIXME
            {
                auto space_view = view_device();
                auto v_view = utopia::view_device(v);

                Device::parallel_for(
                    this->element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        Elem e;
                        space_view.elem(i, e);

                        NodeIndex nodes;
                        space_view.mesh().nodes(i, nodes);
                        const SizeType n_nodes = nodes.size();

                        Point p;
                        for (SizeType i = 0; i < n_nodes; ++i) {
                            auto idx = nodes[i] * mesh_->n_components() + subspace_id_;
                            e.node(i, p);

                            if (r.inside(idx)) {
                                v_view.set(idx, f(p));
                            }
                        }
                    });
            }
        }

        inline void create_interpolation(const FunctionSpace &target, PetscMatrix &I) const {
            mesh_->create_interpolation(target.mesh(), I);
        }

        void read(Input &in) override {
            PetscCommunicator comm;
            if (!mesh_) {
                allocate_mesh(comm);
            }

            mesh_->read(in);
            elements_ = mesh_->elements_ptr();
        }

        bool write(const Path &path, const PetscVector &x) const { return mesh_->write(path, x); }

        inline bool empty() const { return static_cast<bool>(mesh_); }

        inline SizeType component(const SizeType &idx) const {
            const SizeType nc = mesh_->n_components();
            return nc == 1 ? 0 : idx % nc;
        }

        FunctionSpace() : subspace_id_(0) {}

        FunctionSpace(const PetscCommunicator &comm, const SizeType subspace_id = 0) : subspace_id_(subspace_id) {
            allocate_mesh(comm);
        }

        FunctionSpace(const std::shared_ptr<Mesh> &mesh, const SizeType subspace_id = 0)
            : mesh_(mesh), subspace_id_(subspace_id) {
            elements_ = mesh_->elements_ptr();
        }

        void describe() const {
            assert(mesh_);
            mesh_->describe();
        }

    private:
        std::shared_ptr<Mesh> mesh_;
        std::shared_ptr<typename Mesh::Elements> elements_;
        std::vector<std::shared_ptr<DirichletBC>> dirichlet_bcs_;
        SizeType subspace_id_;

        void allocate_mesh(const PetscCommunicator &comm) {
            mesh_ = std::make_shared<Mesh>(comm, is_simplex<UniVarElem_>::value ? DMDA_ELEMENT_P1 : DMDA_ELEMENT_Q1);
            mesh_->set_n_components(NComponents_);

            assert(subspace_id_ + NComponents_ <= mesh_->n_components());
        }
    };

    template <class Point_, class IntArray_, int NComponents_, class UniVarElem_>
    using PetscDMDAFunctionSpace = FunctionSpace<PetscDMDA<Point_, IntArray_>, NComponents_, UniVarElem_>;

    template <class Point_, class IntArray_, int NComponents_, class UniVarElem_>
    const int FunctionSpace<PetscDMDA<Point_, IntArray_>, NComponents_, UniVarElem_>::NComponents;

}  // namespace utopia

#endif  // UTOPIA_PETSC_DMDA_FUNCTION_SPACE_HPP
