#ifndef UTOPIA_PETSC_DM_PLEX_FUNCTION_SPACE_HPP
#define UTOPIA_PETSC_DM_PLEX_FUNCTION_SPACE_HPP

#include "utopia_FunctionSpace.hpp"
#include "utopia_MultiVariateElement.hpp"
#include "utopia_petsc_DMPlex.hpp"

// TOBEREMOVED
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_MakeElem.hpp"
#include "utopia_petsc_dma_FunctionSpace.hpp"

namespace utopia {

    template <class Point_, class IntArray_, int NComponents_, class UniVarElem_>
    class FunctionSpace<PetscDMPlex<Point_, IntArray_>, NComponents_, UniVarElem_> : public Configurable {
    public:
        // concrete types
        using Vector = utopia::PetscVector;
        using Matrix = utopia::PetscMatrix;
        using Comm = utopia::PetscCommunicator;

        // from template arg list
        using Point = Point_;
        using Mesh = utopia::PetscDMPlex<Point_, IntArray_>;
        static const int NComponents = NComponents_;
        using Elem = MultiVariateElem<UniVarElem_, NComponents_>;
        // static const int Dim = Mesh::StaticDim;

        using NodeIndex = typename Mesh::NodeIndex;

        using Scalar = typename Mesh::Scalar;
        using SizeType = typename Mesh::SizeType;

        using ViewDevice = FunctionSpace;
        using Device = typename Mesh::Device;

        using DirichletBC = utopia::DirichletBoundaryCondition<FunctionSpace>;

        template <int NSubVars>
        using Subspace = FunctionSpace<Mesh, NSubVars, UniVarElem_>;

        //////////////////////////////////////////

        inline static constexpr int n_components() { return NComponents; }

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

            mesh_->read(in);
            // elements_ = mesh_->elements_ptr();
        }

        inline void elem(const SizeType &idx, Elem &e) const {
            // FIXME
            e.univar_elem().idx(idx);

            ArrayView<SizeType, 3> nodes;
            mesh_->nodes_local(idx, nodes);

            const SizeType n = nodes.size();

            Point p0, p1, p2;

            if (n == 3) {
                // mesh_->point(nodes[0], p0);
                // mesh_->point(nodes[1], p1);
                // mesh_->point(nodes[2], p2);

                e.set(p0, p1, p2);
            } else {
                assert(false);
            }
        }

        bool write(const Path &path, const PetscVector &x) const { return mesh_->write(path, x); }

        inline bool empty() const { return static_cast<bool>(mesh_); }

        FunctionSpace() : subspace_id_(0) {}

        FunctionSpace(const PetscCommunicator &comm, const SizeType subspace_id = 0) : subspace_id_(subspace_id) {
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

    private:
        std::shared_ptr<Mesh> mesh_;
        // std::shared_ptr<typename Mesh::Elements> elements_;
        // std::vector<std::shared_ptr<DirichletBC>> dirichlet_bcs_;
        SizeType subspace_id_;

        void allocate_mesh(const Comm &comm) { mesh_ = std::make_shared<Mesh>(comm); }
    };  // namespace utopia

    template <class Point_, class IntArray_, int NComponents_, class UniVarElem_>
    using PetscDMPlexFunctionSpace = FunctionSpace<PetscDMPlex<Point_, IntArray_>, NComponents_, UniVarElem_>;

    template <class Point_, class IntArray_, int NComponents_, class UniVarElem_>
    const int FunctionSpace<PetscDMPlex<Point_, IntArray_>, NComponents_, UniVarElem_>::NComponents;

}  // namespace utopia

#endif  // UTOPIA_PETSC_DM_PLEX_FUNCTION_SPACE_HPP
