#ifndef UTOPIA_PETSC_DMDA_DOFMAP_HPP
#define UTOPIA_PETSC_DMDA_DOFMAP_HPP

#include "utopia_AssemblyView.hpp"
#include "utopia_Input.hpp"
#include "utopia_MultiVariateElement.hpp"

// petsc
#include "utopia_petsc_DM.hpp"

namespace utopia {

    template <class Mesh, class Elem_, int NComponents>
    class DofMapping {
        // };

        // template<int Dim, class Elem_, int NComponents>
        // class DofMapping<PetscDM<Dim>, Elem_, NComponents> {
    public:
        //     using Mesh = utopia::PetscStructuredGrid<Dim>;

        static const int NDofs = NComponents * Elem_::NNodes;

        using SizeType = typename Mesh::SizeType;
        using NodeIndex = typename Mesh::NodeIndex;
        using DofIndexNonConst = utopia::ArrayView<SizeType, NDofs>;

        template <class DofIndex>
        static void dofs_local(const Mesh &mesh, const SizeType &var_offset, const SizeType &idx, DofIndex &dofs) {
            if (mesh.n_components() == 1) {
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

        template <class DofIndex>
        static void dofs_local_for_var(const Mesh &mesh, const SizeType &idx, const SizeType &var, DofIndex &dofs) {
            if (mesh.n_components() == 1) {
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

        template <class DofIndex>
        static void dofs(const Mesh &mesh, const SizeType &var_offset, const SizeType &idx, DofIndex &dofs) {
            if (mesh.n_components() == 1) {
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

        template <class DofIndex>
        static void tensorize_dofs(const Mesh &mesh,
                                   const SizeType &var_offset,
                                   const NodeIndex &nodes,
                                   DofIndex &dofs) {
            UTOPIA_UNUSED(mesh);

            // assert(nodes.size() * (mesh.n_components() - var_offset) == dofs.size());
            // assert(NComponents == (mesh.n_components() - var_offset));

            const SizeType n_nodes = nodes.size();
            const SizeType n_components = mesh.n_components();

            SizeType j = 0;
            for (SizeType c = 0; c < NComponents; ++c) {
                const SizeType offset_c = c + var_offset;

                for (SizeType i = 0; i < n_nodes; ++i) {
                    dofs[j++] = nodes[i] * n_components + offset_c;
                }
            }
        }

        template <class DofIndex>
        static void get_var_dofs(const Mesh &mesh, const SizeType &var, const NodeIndex &nodes, DofIndex &dofs) {
            UTOPIA_UNUSED(mesh);

            assert(NComponents == (mesh.n_components() - var));

            const SizeType n_nodes = nodes.size();

            const SizeType n_components = mesh.n_components();

            SizeType j = 0;
            for (SizeType i = 0; i < n_nodes; ++i) {
                dofs[j++] = nodes[i] * n_components + var;
            }
        }

        template <class Elem, class ElementMatrix, class MatView>
        static void add_matrix(const Mesh &mesh,
                               const SizeType &var_offset,
                               const Elem &e,
                               const ElementMatrix &el_mat,
                               MatView &mat) {
            if (mesh.n_components() == 1) {
                assert(var_offset == 0);
                assert(NComponents == 1);

                typename Mesh::NodeIndex dofs;
                mesh.nodes(e.idx(), dofs);

                // Potentially breaks
                mat.atomic_add_matrix(dofs, dofs, &el_mat(0, 0));

            } else {
                DofIndexNonConst indices;
                dofs(mesh, var_offset, e.idx(), indices);

                // Potentially breaks
                mat.atomic_add_matrix(indices, indices, &el_mat(0, 0));
            }
        }

        template <class Elem, class ElementVector, class VecView>
        static void add_vector(const Mesh &mesh,
                               const SizeType &var_offset,
                               const Elem &e,
                               const ElementVector &el_vec,
                               VecView &vec) {
            if (mesh.n_components() == 1) {
                assert(var_offset == 0);
                assert(NComponents == 1);

                typename Mesh::NodeIndex dofs;
                mesh.nodes(e.idx(), dofs);
                vec.atomic_add_vector(dofs, &el_vec(0));

            } else {
                DofIndexNonConst indices;
                dofs(mesh, var_offset, e.idx(), indices);

                vec.atomic_add_vector(indices, &el_vec(0));
            }
        }

        template <class Elem, class ElementVector, class VecView>
        static void set_vector(const Mesh &mesh,
                               const SizeType &var_offset,
                               const Elem &e,
                               const ElementVector &el_vec,
                               VecView &vec) {
            if (mesh.n_components() == 1) {
                assert(var_offset == 0);
                assert(NComponents == 1);

                typename Mesh::NodeIndex dofs;
                mesh.nodes(e.idx(), dofs);

                const SizeType n_dofs = dofs.size();

                for (SizeType i = 0; i < n_dofs; ++i) {
                    vec.atomic_set(dofs[i], el_vec(i));
                }

            } else {
                DofIndexNonConst indices;
                dofs(mesh, var_offset, e.idx(), indices);

                const SizeType n_dofs = indices.size();
                for (SizeType i = 0; i < n_dofs; ++i) {
                    vec.atomic_set(indices[i], el_vec(i));
                }
            }
        }

        template <class Elem, class VectorView, class Values>
        static void local_coefficients(const Mesh &mesh,
                                       const SizeType &var_offset,
                                       const Elem &e,
                                       const VectorView &vec,
                                       Values &values) {
            DofIndexNonConst dofs;
            dofs_local(mesh, var_offset, e.idx(), dofs);
            const SizeType n = dofs.size();

            for (SizeType i = 0; i < n; ++i) {
                assert(dofs[i] < mesh.n_nodes() * mesh.n_components());
                assert(i < values.size());

                values[i] = vec.get(dofs[i]);
            }
        }

        template <class Elem, class VectorView, class Values>
        static void local_coefficients_for_var(const Mesh &mesh,
                                               const Elem &e,
                                               const VectorView &vec,
                                               const SizeType &var,
                                               Values &values) {
            DofIndexNonConst dofs;
            dofs_local_for_var(mesh, e.idx(), var, dofs);
            const SizeType n = dofs.size();

            for (SizeType i = 0; i < n; ++i) {
                assert(dofs[i] < mesh.n_nodes() * mesh.n_components());
                values[i] = vec.get(dofs[i]);
            }
        }
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_DMDA_DOFMAP_HPP
