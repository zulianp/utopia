#ifndef UTOPIA_PETSC_DMA_FUNCTIONSPACE_HPP
#define UTOPIA_PETSC_DMA_FUNCTIONSPACE_HPP

#include "utopia_PetscDM.hpp"
#include "utopia_MultiVariateElement.hpp"
#include "utopia_AssemblyView.hpp"
#include "utopia_Input.hpp"

namespace utopia {

    template<class Mesh, class Elem, int NComponents>
    class DofMap {};

    template<int Dim, class Elem_, int NComponents>
    class DofMap<PetscDM<Dim>, Elem_, NComponents> {
    public:
        static const int NDofs = NComponents * Elem_::NNodes;

        using SizeType  = typename PetscDM<Dim>::SizeType;
        using NodeIndex = typename PetscDM<Dim>::NodeIndex;
        using DofIndexNonConst = utopia::ArrayView<SizeType, NDofs>;

        template<class DofIndex>
        static void dofs_local(const PetscDM<Dim> &mesh, const SizeType &var_offset, const SizeType &idx, DofIndex &dofs)
        {
            if(mesh.n_components() == 1) {
                assert(var_offset == 0);
                assert(NComponents == 1);

                assert(idx < mesh.n_elements());

                NodeIndex nodes;
                mesh.nodes_local(idx, nodes);
                dofs = nodes;

            } else {
                assert(idx < mesh.n_elements());
                NodeIndex nodes;
                mesh.nodes_local(idx, nodes);
                tensorize_dofs(mesh, var_offset, nodes, dofs);
            }
        }

        template<class DofIndex>
        static void dofs_local_for_var(
            const PetscDM<Dim> &mesh,
            const SizeType &idx,
            const SizeType &var,
            DofIndex &dofs
        )
        {
            if(mesh.n_components() == 1) {
                assert(var == 0);
                assert(NComponents == 1);

                assert(idx < mesh.n_elements());

                NodeIndex nodes;
                mesh.nodes_local(idx, nodes);
                dofs = nodes;

            } else {
                assert(idx < mesh.n_elements());
                NodeIndex nodes;
                mesh.nodes_local(idx, nodes);
                get_var_dofs(mesh, var, nodes, dofs);
            }
        }

        template<class DofIndex>
        static void dofs(const PetscDM<Dim> &mesh, const SizeType &var_offset, const SizeType &idx, DofIndex &dofs)
        {
            if(mesh.n_components() == 1) {
                assert(var_offset == 0);
                assert(NComponents == 1);

                NodeIndex nodes;
                mesh.nodes(idx, nodes);
                dofs = nodes;

            } else {
                assert(idx < mesh.n_elements());
                NodeIndex nodes;
                mesh.nodes(idx, nodes);
                tensorize_dofs(mesh, var_offset, nodes, dofs);
            }
        }

        template<class DofIndex>
        static void tensorize_dofs(
            const PetscDM<Dim> &mesh,
            const SizeType &var_offset,
            const NodeIndex &nodes,
            DofIndex &dofs)
        {
            UTOPIA_UNUSED(mesh);

            // assert(nodes.size() * (mesh.n_components() - var_offset) == dofs.size());
            // assert(NComponents == (mesh.n_components() - var_offset));

            const SizeType n_nodes = nodes.size();
            const SizeType n_components = mesh.n_components();

            SizeType j = 0;
            for(SizeType c = 0; c < NComponents; ++c) {
                const SizeType offset_c = c + var_offset;

                for(SizeType i = 0; i < n_nodes; ++i) {
                    dofs[j++] = nodes[i] * n_components + offset_c;
                }
            }
        }

        template<class DofIndex>
        static void get_var_dofs(
            const PetscDM<Dim> &mesh,
            const SizeType &var,
            const NodeIndex &nodes,
            DofIndex &dofs)
        {
            UTOPIA_UNUSED(mesh);

            assert(NComponents == (mesh.n_components() - var));

            const SizeType n_nodes = nodes.size();

            const SizeType n_components = mesh.n_components();

            SizeType j = 0;
            for(SizeType i = 0; i < n_nodes; ++i) {
                dofs[j++] = nodes[i] * n_components + var;
            }

        }

        template<class Elem, class ElementMatrix, class MatView>
        static void add_matrix(
            const PetscDM<Dim> &mesh,
            const SizeType &var_offset,
            const Elem &e,
            const ElementMatrix &el_mat,
            MatView &mat)
        {
            if(mesh.n_components() == 1) {
                assert(var_offset == 0);
                assert(NComponents == 1);

                const SizeType n_dofs = e.nodes().size();
                const auto &dofs = e.nodes();

                for(SizeType i = 0; i < n_dofs; ++i) {
                    for(SizeType j = 0; j < n_dofs; ++j) {
                        mat.atomic_add(dofs[i], dofs[j], el_mat(i, j));
                    }
                }
            } else {
                DofIndexNonConst indices;
                dofs(mesh, var_offset, e.idx(), indices);

                // const SizeType n_dofs = indices.size();
                // for(SizeType i = 0; i < n_dofs; ++i) {
                //     for(SizeType j = 0; j < n_dofs; ++j) {
                //         mat.atomic_add(indices[i], indices[j], el_mat(i, j));
                //     }
                // }

                //Potentially breaks
                mat.atomic_add_matrix(indices, indices, &el_mat(0, 0));
            }
       }

       template<class Elem, class ElementVector, class VecView>
       static void add_vector(
        const PetscDM<Dim> &mesh,
        const SizeType &var_offset,
        const Elem &e,
        const ElementVector &el_vec,
        VecView &vec)
       {
            if(mesh.n_components() == 1) {
                assert(var_offset == 0);
                assert(NComponents == 1);

                const SizeType n_dofs = e.nodes().size();
                const auto &dofs = e.nodes();

                for(SizeType i = 0; i < n_dofs; ++i) {
                    vec.atomic_add(dofs[i], el_vec(i));
                }

            } else {
                DofIndexNonConst indices;
                dofs(mesh, var_offset, e.idx(), indices);

                const SizeType n_dofs = indices.size();
                for(SizeType i = 0; i < n_dofs; ++i) {
                    vec.atomic_add(indices[i], el_vec(i));
                }
            }
        }


        template<class Elem, class ElementVector, class VecView>
        static void set_vector(
         const PetscDM<Dim> &mesh,
         const SizeType &var_offset,
         const Elem &e,
         const ElementVector &el_vec,
         VecView &vec)
        {
             if(mesh.n_components() == 1) {
                 assert(var_offset == 0);
                 assert(NComponents == 1);

                 const SizeType n_dofs = e.nodes().size();
                 const auto &dofs = e.nodes();

                 for(SizeType i = 0; i < n_dofs; ++i) {
                     vec.atomic_set(dofs[i], el_vec(i));
                 }

             } else {
                 DofIndexNonConst indices;
                 dofs(mesh, var_offset, e.idx(), indices);

                 const SizeType n_dofs = indices.size();
                 for(SizeType i = 0; i < n_dofs; ++i) {
                     vec.atomic_set(indices[i], el_vec(i));
                 }
             }
         }

       template<class Elem, class VectorView, class Values>
       static void local_coefficients(
           const PetscDM<Dim> &mesh,
           const SizeType &var_offset,
           const Elem &e,
           const VectorView &vec,
           Values &values)
       {
            DofIndexNonConst dofs;
            dofs_local(mesh, var_offset, e.idx(), dofs);
            const SizeType n = dofs.size();

            for(SizeType i = 0; i < n; ++i) {
                assert(dofs[i] < mesh.n_nodes() * mesh.n_components());
                values[i] = vec.get(dofs[i]);
            }
       }

       template<class Elem, class VectorView, class Values>
       static void local_coefficients_for_var(
           const PetscDM<Dim> &mesh,
           const Elem &e,
           const VectorView &vec,
           const SizeType &var,
           Values &values)
       {
            DofIndexNonConst dofs;
            dofs_local_for_var(mesh, e.idx(), var, dofs);
            const SizeType n = dofs.size();

            for(SizeType i = 0; i < n; ++i) {
                assert(dofs[i] < mesh.n_nodes() * mesh.n_components());
                values[i] = vec.get(dofs[i]);
            }
       }
    };

    template<class Elem_, int NComponents_>
    class FunctionSpace<PetscDM<Elem_::Dim>, NComponents_, Elem_> : public Configurable {
    public:
        static const int Dim = Elem_::Dim;
        static const std::size_t UDim = Dim;
        static const int NComponents = NComponents_;
        using Mesh = utopia::PetscDM<Dim>;
        using Elem = MultiVariateElem<Elem_, NComponents>;
        using Shape = Elem_;
        using MemType = typename Elem::MemType;
        using Scalar = typename Mesh::Scalar;
        using SizeType = typename Mesh::SizeType;
        using Point = typename Mesh::Point;

        using ViewDevice = FunctionSpace;
        using Device = typename Mesh::Device;
        // using DofIndex = typename Mesh::NodeIndex;
        using Vector = utopia::PetscVector;
        using Matrix = utopia::PetscMatrix;
        using Comm = utopia::PetscCommunicator;
        using DirichletBC = utopia::DirichletBoundaryCondition<FunctionSpace>;
        using DofMap      = utopia::DofMap<PetscDM<Dim>, Elem_, NComponents>;
        static const int NDofs = DofMap::NDofs;

        template<int NSubVars>
        using Subspace = FunctionSpace<PetscDM<Elem_::Dim>, NSubVars, Elem_>;

        void read(Input &in) override;
        bool write(const Path &path, const PetscVector &x) const;

        FunctionSpace(const std::shared_ptr<Mesh> &mesh, const SizeType &subspace_id = 0)
        : mesh_(mesh), subspace_id_(subspace_id)
        {}

        FunctionSpace(Mesh &mesh, const SizeType &subspace_id = 0)
        : mesh_(utopia::make_ref(mesh)), subspace_id_(subspace_id)
        {}

        ///shallow copy
        FunctionSpace(const FunctionSpace &other)
        : mesh_(other.mesh_), dirichlet_bcs_(other.dirichlet_bcs_), subspace_id_(other.subspace_id_)
        {}

        FunctionSpace()
        : subspace_id_(0)
        {}

        inline static constexpr int n_components() { return NComponents; }

        template<class Quadrature>
        ShapeFunction<FunctionSpace, Quadrature> shape(const Quadrature &q)
        {
            return ShapeFunction<FunctionSpace, Quadrature>(*this, q);
        }

        template<class Quadrature>
        PhysicalGradient<FunctionSpace, Quadrature> shape_grad(const Quadrature &q)
        {
            return PhysicalGradient<FunctionSpace, Quadrature>(*this, q);
        }


        template<class Quadrature>
        Differential<FunctionSpace, Quadrature> differential(const Quadrature &q)
        {
            return Differential<FunctionSpace, Quadrature>(*this, q);
        }

        template<class... Args>
        void build(
            const PetscCommunicator     &comm,
            const std::array<SizeType, UDim> &dims,
            const std::array<Scalar, UDim>   &box_min,
            const std::array<Scalar, UDim>   &box_max,
            const SizeType &subspace_id = 0)
        {
            mesh_ = std::make_shared<Mesh>(comm, dims, box_min, box_max, NComponents);
            subspace_id_ = subspace_id;
        }

        FunctionSpace<PetscDM<Elem_::Dim>, 1, Elem_> subspace(const SizeType &i) const
        {
            FunctionSpace<PetscDM<Elem_::Dim>, 1, Elem_> space(mesh_, subspace_id_ + i);
            // space.set_dirichlet_conditions(dirichlet_bcs_);
            assert(i < NComponents);
            assert(i + subspace_id_ < mesh_->n_components());
            return space;
        }


        template<int NVars>
        void subspace(const SizeType &i, FunctionSpace<PetscDM<Elem_::Dim>, NVars, Elem_> &space) const
        {
            space.set_mesh(mesh_);
            space.set_subspace_id(subspace_id_ + i);
        }

        template<int NVars>
        FunctionSpace<PetscDM<Elem_::Dim>, NVars, Elem_> vector_subspace(const SizeType &i) const
        {
            FunctionSpace<PetscDM<Elem_::Dim>, NVars, Elem_> space(mesh_, subspace_id_ + i);
            // space.set_dirichlet_conditions(dirichlet_bcs_);
            assert(i + NVars < NComponents);
            assert(subspace_id_ + i < mesh_->n_components());
            return space;
        }

        template<class Fun>
        void each_element(Fun fun)
        {
            mesh_->each_element(fun);
        }

        void elem(const SizeType &idx, Elem &e) const;


        bool is_boundary_dof(const SizeType &idx) const
        {
            return mesh_->is_boundary(idx);
        }

        SideSet::BoundaryIdType boundary_id(const SizeType &idx) const
        {
            return mesh_->boundary_id(idx);
        }

        inline SizeType component(const SizeType &idx) const
        {
            return mesh_->n_components() == 1? 0 : idx % mesh_->n_components();
        }

        const ViewDevice &view_device() const
        {
            return *this;
        }

        Range local_element_range() const
        {
            return mesh_->local_element_range();
        }

        void create_matrix(PetscMatrix &mat) const
        {
            mesh_->create_matrix(mat);
        }

        void create_vector(PetscVector &vec) const
        {
            mesh_->create_vector(vec);
        }

        static DeviceView<PetscMatrix, 2> assembly_view_device(PetscMatrix &mat)
        {
            return DeviceView<PetscMatrix, 2>(mat, utopia::GLOBAL_ADD);
        }

        static DeviceView<PetscVector, 1> assembly_view_device(PetscVector &vec)
        {
            return DeviceView<PetscVector, 1>(vec, utopia::GLOBAL_ADD);
        }

        static DeviceView<const PetscVector, 1> assembly_view_device(const PetscVector &vec)
        {
            return DeviceView<const PetscVector, 1>(vec);
        }

        template<class DofIndex>
        void dofs(const SizeType &idx, DofIndex &dofs) const
        {
            DofMap::dofs(*mesh_, subspace_id_, idx, dofs);
        }

        template<class DofIndex>
        void dofs_local(const SizeType &idx, DofIndex &dofs) const
        {
           DofMap::dofs_local(*mesh_, subspace_id_, idx, dofs);
        }

        template<class ElementMatrix, class MatView>
        void add_matrix(
            const Elem &e,
            const ElementMatrix &el_mat,
            MatView &mat) const
        {
            DofMap::add_matrix(*mesh_, subspace_id_, e, el_mat, mat);
        }

        template<class ElementVector, class VecView>
        void add_vector(
            const Elem &e,
            const ElementVector &el_vec,
            VecView &vec) const
        {
            DofMap::add_vector(*mesh_, subspace_id_, e, el_vec, vec);
        }

        template<class ElementVector, class VecView>
        void set_vector(
            const Elem &e,
            const ElementVector &el_vec,
            VecView &vec) const
        {
            DofMap::set_vector(*mesh_, subspace_id_, e, el_vec, vec);
        }

        template<class VectorView, class Values>
        void local_coefficients(
            const Elem &e,
            const VectorView &vec,
            Values &values) const
        {
            DofMap::local_coefficients(*mesh_, subspace_id_, e, vec, values);
        }

        template<class VectorView, class Values>
        void local_coefficients(
            const Elem &e,
            const VectorView &vec,
            const SizeType &var,
            Values &values) const
        {
            DofMap::local_coefficients_for_var(*mesh_, e, vec, subspace_id_ +var, values);
        }

        const Mesh &mesh() const
        {
            return *mesh_;
        }

        Mesh &mesh()
        {
            return *mesh_;
        }

        PetscCommunicator &comm()
        {
            return mesh_->comm();
        }

        const PetscCommunicator &comm() const
        {
            return mesh_->comm();
        }

        inline SizeType n_dofs() const
        {
            return mesh_->n_nodes() * NComponents;
        }

        void set_mesh(const std::shared_ptr<Mesh> &mesh)
        {
            mesh_ = mesh;
        }

        void set_subspace_id(const SizeType &i)
        {
            subspace_id_ = i;
        }

        void set_dirichlet_conditions(const std::vector<std::shared_ptr<DirichletBC>> &conds)
        {
            dirichlet_bcs_ = conds;
        }

        template<class... Args>
        void emplace_dirichlet_condition(Args && ...args)
        {
            dirichlet_bcs_.push_back(utopia::make_unique<DirichletBC>(*this, std::forward<Args>(args)...));
        }

        void apply_constraints(PetscMatrix &mat, PetscVector &vec) const
        {
            for(const auto &bc : dirichlet_bcs_) {
                bc->apply(mat, vec);
            }
        }

        void apply_constraints(PetscMatrix &mat) const
        {
            for(const auto &bc : dirichlet_bcs_) {
                bc->apply(mat);
            }
        }

        void apply_constraints(PetscVector &vec) const
        {
            for(const auto &bc : dirichlet_bcs_) {
                bc->apply(vec);
            }
        }

        void apply_zero_constraints(PetscVector &vec) const
        {
            for(const auto &bc : dirichlet_bcs_) {
                bc->apply_zero(vec);
            }
        }

        inline bool empty() const
        {
            return !static_cast<bool>(mesh_);
        }


        template<class F>
        void sample(Vector &v, F f, const int c = 0)
        {
            auto r = v.range();
            // auto n = r.extent() * NComponents;
            assert(!v.empty());

            // Write<Vector> w(v, utopia::AUTO);

            // Point p;
            // for(auto i = r.begin(); i < r.end(); ++i) {
            //     this->mesh().node(i/NComponents, p);
            //     v.set(i, f(p));
            // }

            //FIXME
            {
                auto space_view = view_device();
                auto v_view = utopia::view_device(v);

                Device::parallel_for(
                    this->local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                {
                    Elem e;
                    space_view.elem(i, e);


                    const SizeType n_nodes = e.n_nodes();

                    Point p;
                    for(SizeType i = 0; i < n_nodes; ++i) {
                        auto idx = e.node_id(i) * mesh_->n_components() + subspace_id_;
                        e.node(i, p);

                        if(r.inside(idx)) {
                            v_view.set(idx, f(p));
                        }
                    }
                });
            }
        }

    private:
        std::shared_ptr<Mesh> mesh_;
        std::vector<std::shared_ptr<DirichletBC>> dirichlet_bcs_;
        SizeType subspace_id_;
    };

    template<class Elem_, int NComponents_>
    const int FunctionSpace<PetscDM<Elem_::Dim>, NComponents_, Elem_>::NComponents;
}

#endif //UTOPIA_PETSC_DMA_FUNCTIONSPACE_HPP
