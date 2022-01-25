#ifndef UTOPIA_PETSC_DOF_MAP_HPP
#define UTOPIA_PETSC_DOF_MAP_HPP

#include "utopia_Input.hpp"

// petsc
#include "utopia_petsc_DMDABase.hpp"
#include "utopia_petsc_StructuredGrid.hpp"

namespace utopia {

    namespace petsc {

        class DofMap {
        public:
            using Mesh = utopia::petsc::StructuredGrid;
            using SizeType = Traits<Mesh>::SizeType;
            using ElemToNodeIndex = DMDABase::ElemToNodeIndex;

            DofMap(const Mesh &mesh)
                : mesh_(mesh),
                  elem_to_node_index_(mesh.dm().elem_to_node_index()),
                  elem_to_local_node_index_(mesh.dm().elem_to_local_node_index()) {}

            inline SizeType n_node_x_element() const { return elem_to_node_index_.extent(1); }
            inline SizeType n_elements() const { return elem_to_node_index_.extent(0); }
            inline SizeType n_components() const { return mesh_.dm().get_dof(); }

            template <class DofIndex>
            void dofs_local(const SizeType &var_offset, const SizeType &idx, DofIndex &dofs) const {
                assert(idx < n_elements());
                const SizeType n = n_node_x_element();

                if (this->n_components() == 1) {
                    assertt(n == dofs.size());
                    for (SizeType k = 0; k < n; ++k) {
                        dofs[k] = elem_to_local_node_index_(idx, k);
                    }

                } else {
                    tensorize_dofs(elem_to_local_node_index_, idx, var_offset, dofs);
                }
            }

            template <class DofIndex>
            void tensorize_dofs(const ElemToNodeIndex &nodes,
                                const SizeType &idx,
                                const SizeType &var_offset,
                                DofIndex &dofs) const {
                const SizeType n_nodes = n_node_x_element();
                const SizeType n_components = this->n_components();

                assert(dofs.size() == n_nodes * n_components);

                SizeType j = 0;
                for (SizeType c = 0; c < n_components; ++c) {
                    const SizeType offset_c = c + var_offset;

                    for (SizeType i = 0; i < n_nodes; ++i) {
                        dofs[j++] = nodes(idx, i) * n_components + offset_c;
                    }
                }
            }

            // template <class DofIndex>
            // static void dofs_local_for_var(const Mesh &mesh, const SizeType &idx, const SizeType &var, DofIndex
            // &dofs) {
            //     if (mesh.n_components() == 1) {
            //         assert(var == 0);
            //         assert(NComponents == 1);

            //         assert(idx < mesh.n_elements());

            //         NodeIndex nodes;
            //         mesh.nodes_local(idx, nodes);
            //         dofs = nodes;

            //     } else {
            //         assert(idx < mesh.n_elements());
            //         NodeIndex nodes;
            //         mesh.nodes_local(idx, nodes);
            //         get_var_dofs(mesh, var, nodes, dofs);
            //     }
            // }

            // template <class DofIndex>
            // static void dofs(const Mesh &mesh, const SizeType &var_offset, const SizeType &idx, DofIndex &dofs) {
            //     if (mesh.n_components() == 1) {
            //         assert(var_offset == 0);
            //         assert(NComponents == 1);

            //         NodeIndex nodes;
            //         mesh.nodes(idx, nodes);
            //         dofs = nodes;

            //     } else {
            //         assert(idx < mesh.n_elements());
            //         NodeIndex nodes;
            //         mesh.nodes(idx, nodes);
            //         tensorize_dofs(mesh, var_offset, nodes, dofs);
            //     }
            // }

            // template <class DofIndex>
            // static void get_var_dofs(const Mesh &mesh, const SizeType &var, const NodeIndex &nodes, DofIndex &dofs) {
            //     UTOPIA_UNUSED(mesh);

            //     assert(NComponents == (mesh.n_components() - var));

            //     const SizeType n_nodes = nodes.size();

            //     const SizeType n_components = mesh.n_components();

            //     SizeType j = 0;
            //     for (SizeType i = 0; i < n_nodes; ++i) {
            //         dofs[j++] = nodes[i] * n_components + var;
            //     }
            // }

            // template <class Elem, class ElementMatrix, class MatView>
            // static void add_matrix(const Mesh &mesh,
            //                        const SizeType &var_offset,
            //                        const Elem &e,
            //                        const ElementMatrix &el_mat,
            //                        MatView &mat) {
            //     if (mesh.n_components() == 1) {
            //         assert(var_offset == 0);
            //         assert(NComponents == 1);

            //         typename Mesh::NodeIndex dofs;
            //         mesh.nodes(e.idx(), dofs);

            //         // Potentially breaks
            //         mat.atomic_add_matrix(dofs, dofs, &el_mat(0, 0));

            //     } else {
            //         DofIndexNonConst indices;
            //         dofs(mesh, var_offset, e.idx(), indices);

            //         // Potentially breaks
            //         mat.atomic_add_matrix(indices, indices, &el_mat(0, 0));
            //     }
            // }

            template <class ElementVector, class VecView>
            void add_vector(const SizeType idx,
                            const SizeType &var_offset,
                            const ElementVector &el_vec,
                            VecView &vec) const {
                const int n_components = this->n_components();
                const SizeType n_nodes = n_node_x_element();

                if (n_components == 1) {
                    for (SizeType k = 0; k < n_nodes; ++k) {
                        vec.atomic_add(elem_to_node_index_(idx, k), el_vec(k));
                    }

                } else {
                    SizeType j = 0;
                    for (SizeType c = 0; c < n_components; ++c) {
                        const SizeType offset_c = c + var_offset;

                        for (SizeType i = 0; i < n_nodes; ++i, ++j) {
                            const SizeType dof = elem_to_node_index_(idx, i) * n_components + offset_c;
                            const auto value = el_vec[j];
                            vec.atomic_add(dof, value);
                        }
                    }
                }
            }

            // template <class Elem, class ElementVector, class VecView>
            // static void set_vector(const Mesh &mesh,
            //                        const SizeType &var_offset,
            //                        const Elem &e,
            //                        const ElementVector &el_vec,
            //                        VecView &vec) {
            //     if (mesh.n_components() == 1) {
            //         assert(var_offset == 0);
            //         assert(NComponents == 1);

            //         typename Mesh::NodeIndex dofs;
            //         mesh.nodes(e.idx(), dofs);

            //         const SizeType n_dofs = dofs.size();

            //         for (SizeType i = 0; i < n_dofs; ++i) {
            //             vec.atomic_set(dofs[i], el_vec(i));
            //         }

            //     } else {
            //         DofIndexNonConst indices;
            //         dofs(mesh, var_offset, e.idx(), indices);

            //         const SizeType n_dofs = indices.size();
            //         for (SizeType i = 0; i < n_dofs; ++i) {
            //             vec.atomic_set(indices[i], el_vec(i));
            //         }
            //     }
            // }

            // template <class Elem, class VectorView, class Values>
            // static void local_coefficients(const Mesh &mesh,
            //                                const SizeType &var_offset,
            //                                const Elem &e,
            //                                const VectorView &vec,
            //                                Values &values) {
            //     DofIndexNonConst dofs;
            //     dofs_local(mesh, var_offset, e.idx(), dofs);
            //     const SizeType n = dofs.size();

            //     for (SizeType i = 0; i < n; ++i) {
            //         assert(SizeType(dofs[i]) < SizeType(mesh.n_nodes() * mesh.n_components()));
            //         assert(i < SizeType(values.size()));

            //         values[i] = vec.get(dofs[i]);
            //     }
            // }

            // template <class Elem, class VectorView, class Values>
            // static void local_coefficients_for_var(const Mesh &mesh,
            //                                        const Elem &e,
            //                                        const VectorView &vec,
            //                                        const SizeType &var,
            //                                        Values &values) {
            //     DofIndexNonConst dofs;
            //     dofs_local_for_var(mesh, e.idx(), var, dofs);
            //     const SizeType n = dofs.size();

            //     for (SizeType i = 0; i < n; ++i) {
            //         assert(SizeType(dofs[i]) < SizeType(mesh.n_nodes() * mesh.n_components()));
            //         values[i] = vec.get(dofs[i]);
            //     }
            // }

        private:
            const Mesh &mesh_;
            ElemToNodeIndex elem_to_node_index_;
            ElemToNodeIndex elem_to_local_node_index_;
        };

    }  // namespace petsc
}  // namespace utopia

#endif  // UTOPIA_PETSC_DOF_MAP_HPP
