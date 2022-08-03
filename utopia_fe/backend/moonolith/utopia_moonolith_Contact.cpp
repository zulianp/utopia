#include "utopia_moonolith_Contact.hpp"

// Utils from utopia
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_NormalizeRows.hpp"
#include "utopia_ZeroRowsToIdentity.hpp"

// Moonolith
#include "moonolith_affine_transform.hpp"
#include "moonolith_assign_functions.hpp"
#include "moonolith_contact.hpp"
#include "moonolith_elem_quad.hpp"
#include "moonolith_elem_segment.hpp"
#include "moonolith_elem_shape.hpp"
#include "moonolith_elem_triangle.hpp"
#include "moonolith_matlab_scripter.hpp"
#include "moonolith_par_contact.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_sparse_matrix.hpp"

// Integrate
// #include "utopia_moonolith_permutations.hpp"

// #include <memory>
// #include <vector>

namespace utopia {
    namespace moonolith {

        class Contact::Impl {
        public:
            bool write() {
                if (coupling_matrix && !empty(*coupling_matrix)) {
                    rename("b", *coupling_matrix);
                    utopia::write("B.m", *coupling_matrix);
                }

                if (mass_matrix && !empty(*mass_matrix)) {
                    rename("d", *mass_matrix);
                    utopia::write("D.m", *mass_matrix);
                }

                if (basis_transform && !empty(*basis_transform)) {
                    rename("q", *basis_transform);
                    utopia::write("Q.m", *basis_transform);
                }

                if (basis_transform_inv && !empty(*basis_transform_inv)) {
                    rename("q_inv", *basis_transform_inv);
                    utopia::write("Q_inv.m", *basis_transform_inv);
                }

                if (complete_transformation && !empty(*complete_transformation)) {
                    rename("t", *complete_transformation);
                    utopia::write("T.m", *complete_transformation);
                }

                if (!empty(inv_mass_vector)) {
                    rename("imv", inv_mass_vector);
                    utopia::write("inv_mass_vector.m", inv_mass_vector);
                }

                if (!empty(weighted_gap)) {
                    rename("wg", weighted_gap);

                    utopia::write("weighted_gap.m", weighted_gap);
                }

                if (!empty(gap)) {
                    rename("g", gap);
                    utopia::write("gap.m", gap);
                }

                if (!empty(weighted_normal)) {
                    rename("wn", weighted_normal);
                    utopia::write("weighted_normal.m", weighted_normal);
                }

                if (!empty(normal)) {
                    rename("n", normal);
                    utopia::write("normal.m", normal);
                }

                if (!empty(is_contact)) {
                    rename("ic", is_contact);
                    utopia::write("is_contact.m", is_contact);
                }

                if (!empty(is_glue)) {
                    rename("ig", is_glue);
                    utopia::write("is_glue.m", is_glue);
                }

                return true;
            }

            std::shared_ptr<Matrix> coupling_matrix, mass_matrix;
            std::shared_ptr<Matrix> basis_transform, basis_transform_inv;
            std::shared_ptr<Matrix> orthogonal_trafo, complete_transformation;

            Vector inv_mass_vector;
            Vector weighted_gap, gap;
            Vector weighted_normal, normal;
            Vector is_contact;
            Vector is_glue;
        };

        template <int Dim>
        class Contact::ImplD : public Contact::Impl {
        public:
            // void finalize(const moonolith::ContactBuffers<double> &buffers,
            //               const SpaceT<Dim> &element_wise_space,
            //               const SpaceT<Dim> &node_wise_space) {
            //     utopia::convert(buffers.B, element_wise.B);
            //     utopia::convert(buffers.D, element_wise.D);
            //     utopia::convert(buffers.Q, element_wise.Q);
            //     utopia::convert(buffers.basis_transform_inv, element_wise.basis_transform_inv);

            //     convert_tensor(buffers.gap, element_wise.weighted_gap);
            //     convert_tensor(buffers.normal, element_wise.weighted_normal);
            //     convert_tensor(buffers.is_glue, element_wise.is_glue);
            //     convert_tensor(buffers.is_contact, element_wise.is_contact);

            //     convert_to_node_wise(element_wise_space, element_wise, node_wise_space, *node_wise);
            // }

            //         static void normalize_and_build_orthgonal_trafo(const SpaceT<Dim> &node_wise_space,
            //                                                         ConvertContactTensors &node_wise) {
            //             node_wise.orthogonal_trafo =
            //                 local_sparse(node_wise_space.dof_map().n_local_dofs(),
            //                 node_wise_space.dof_map().n_local_dofs(), Dim);

            //             moonolith::HouseholderTransformation<double, Dim> H;
            //             auto &normal = node_wise.normal;

            //             auto r = range(normal);

            //             ReadAndWrite<Vector> rw_normal(normal);
            //             Read<Vector> r_is_c(node_wise.is_contact);
            //             Write<Matrix> w_ot(node_wise.orthogonal_trafo, utopia::LOCAL);

            //             moonolith::Vector<double, Dim> n;
            //             for (auto i = r.begin(); i < r.end(); i += Dim) {
            //                 const bool is_contact = node_wise.is_contact.get(i);

            //                 if (is_contact) {
            //                     for (int d = 0; d < Dim; ++d) {
            //                         n[d] = normal.get(i + d);
            //                     }

            //                     auto len = length(n);
            //                     assert(len > 0.0);

            //                     n /= len;

            //                     for (int d = 0; d < Dim; ++d) {
            //                         normal.set(i + d, n[d]);
            //                     }

            //                     n.x -= 1.0;

            //                     len = length(n);

            //                     if (approxeq(len, 0.0, 1e-15)) {
            //                         H.identity();
            //                     } else {
            //                         n /= len;
            //                         H.init(n);
            //                     }

            //                     assert(approxeq(std::abs(measure(H)), 1.0, 1e-8));

            //                     for (int d1 = 0; d1 < Dim; ++d1) {
            //                         for (int d2 = 0; d2 < Dim; ++d2) {
            //                             node_wise.orthogonal_trafo.set(i + d1, i + d2, H(d1, d2));
            //                         }
            //                     }

            //                 } else {
            //                     for (int d1 = 0; d1 < Dim; ++d1) {
            //                         node_wise.orthogonal_trafo.set(i + d1, i + d1, 1.0);
            //                     }
            //                 }
            //             }
            //         }

            //         static void convert_to_node_wise(const SpaceT<Dim> &elem_wise_space,
            //                                          ConvertContactTensors &elem_wise,
            //                                          const SpaceT<Dim> &node_wise_space,
            //                                          ConvertContactTensors &node_wise) {
            //             Matrix perm, vector_perm;

            //             make_permutation(elem_wise_space, node_wise_space, perm);

            //             make_vector_permutation(Dim, elem_wise_space, node_wise_space, vector_perm);

            //             node_wise.weighted_gap = perm * elem_wise.weighted_gap;
            //             node_wise.is_glue = perm * elem_wise.is_glue;
            //             node_wise.is_contact = perm * elem_wise.is_contact;

            //             node_wise.weighted_normal = vector_perm * elem_wise.weighted_normal;

            //             Matrix B_x = perm * elem_wise.B * transpose(perm);
            //             Matrix D_x = perm * elem_wise.D * transpose(perm);
            //             Matrix Q_x = perm * elem_wise.Q * transpose(perm);

            //             Matrix basis_transform_inv_x = perm * elem_wise.basis_transform_inv * transpose(perm);

            //             normalize_rows(Q_x, 1e-15);
            //             normalize_rows(basis_transform_inv_x, 1e-15);

            //             {
            //                 each_transform(node_wise.is_contact,
            //                                node_wise.is_contact,
            //                                [&node_wise](const SizeType i, const double value) -> double {
            //                                    if (value > 0.0) {
            //                                        return 1.0;
            //                                    } else {
            //                                        return 0.0;
            //                                    }
            //                                });

            //                 each_transform(
            //                     node_wise.is_glue, node_wise.is_glue, [&node_wise](const SizeType i, const double
            //                     value)
            //                     -> double {
            //                         if (value > 0.0) {
            //                             return 1.0;
            //                         } else {
            //                             return 0.0;
            //                         }
            //                     });
            //             }

            //             node_wise.inv_mass_vector = sum(D_x, 1);

            //             e_pseudo_inv(node_wise.inv_mass_vector, node_wise.inv_mass_vector, 1e-15);
            //             Matrix D_inv_x = diag(node_wise.inv_mass_vector);

            //             Matrix T_temp_x = D_inv_x * B_x;
            //             Matrix T_x = Q_x * T_temp_x;

            //             tensorize(B_x, Dim, node_wise.B);
            //             tensorize(D_x, Dim, node_wise.D);
            //             tensorize(T_x, Dim, node_wise.T);
            //             tensorize(Q_x, Dim, node_wise.Q);
            //             tensorize(basis_transform_inv_x, Dim, node_wise.basis_transform_inv);

            //             tensorize(Dim, node_wise.inv_mass_vector);
            //             tensorize(Dim, node_wise.is_glue);

            //             normalize_rows(node_wise.T, 1e-15);

            //             assert(check_op(node_wise.T));

            //             node_wise.T += local_identity(local_size(node_wise.T));

            //             node_wise.gap = node_wise.Q * e_mul(node_wise.inv_mass_vector, node_wise.weighted_gap);
            //             node_wise.normal = node_wise.Q * e_mul(node_wise.inv_mass_vector, node_wise.weighted_normal);

            //             normalize_and_build_orthgonal_trafo(node_wise_space, node_wise);
            //             node_wise.complete_transformation = node_wise.T * node_wise.orthogonal_trafo;

            //             {
            //                 static const double LARGE_VALUE = 1e6;
            //                 Write<Vector> w(node_wise.gap);
            //                 each_read(node_wise.is_contact, [&node_wise](const SizeType i, const double value) {
            //                     if (value == 0.0) {
            //                         node_wise.gap.set(i, LARGE_VALUE);
            //                     }
            //                 });
            //             }

            //             // node_wise.write();
            //         }

            //         ConvertContactTensors element_wise;
            //         std::shared_ptr<ConvertContactTensors> node_wise;

            //         ConvertContactBuffers(MPI_Comm comm) { node_wise = std::make_shared<ConvertContactTensors>(); }

            //         using MeshT = moonolith::Mesh<double, Dim>;
            //         using SpaceT = moonolith::FunctionSpace<MeshT>;

            //         static bool apply(const ContactParams &params,
            //                           const std::shared_ptr<ElementBlackList> &black_list,
            //                           libMesh::MeshBase &mesh,
            //                           libMesh::DofMap &dof_map,
            //                           ConvertContactBuffers &contact_data) {
            //             using AlogrithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim,
            //             LibMeshFunctionSpaceAdapter>; using Adapter = typename AlogrithmT::Adapter;

            //             moonolith::Communicator comm = mesh.comm().get();

            //             auto mesh_ptr = std::make_shared<MeshT>(comm);
            //             SpaceT space(mesh_ptr);

            //             extract_trace_space(mesh, dof_map, params.variable_number, space);
            //             // moonolith::MatlabScripter script;
            //             // mesh_ptr->draw(script);
            //             // script.save("contact.m");

            //             SpaceT elem_wise_space;
            //             space.separate_dofs(elem_wise_space);

            //             moonolith::ParContact<double, Dim> par_contact(comm, Dim == 2);

            //             if (par_contact.assemble(
            //                     params.contact_pair_tags, elem_wise_space, params.side_set_search_radius,
            //                     params.is_glue)) {
            //                 contact_data.finalize(par_contact.buffers, elem_wise_space, space);
            //                 return true;
            //             } else {
            //                 return false;
            //             }
            //         }
            //     };
        };

        //     bool Contact::assemble(libMesh::MeshBase &mesh,
        //                                            libMesh::DofMap &dof_map,
        //                                            const ContactParams &params) {
        //         Chrono overall_time;
        //         overall_time.start();

        //         ConvertContactBuffers contact_data(mesh.comm().get());

        //         const int spatial_dim = mesh.spatial_dimension();

        //         has_contact_ = false;
        //         if (spatial_dim == 2) {
        //             has_contact_ = ConvertContactAlgorithm<2>::apply(params, black_list_, mesh, dof_map,
        //             contact_data);
        //         } else if (spatial_dim == 3) {
        //             has_contact_ = ConvertContactAlgorithm<3>::apply(params, black_list_, mesh, dof_map,
        //             contact_data);
        //         }

        //         if (has_contact_) {
        //             contact_tensors_ = contact_data.node_wise;
        //         } else {
        //             // init default
        //             init_no_contact(mesh, dof_map);
        //         }

        //         has_glue_ = has_contact_ && !params.is_glue->empty();

        //         if (has_glue_) {
        //             double n_glued = sum(contact_tensors_->is_glue);

        //             if (n_glued == 0.) {
        //                 has_glue_ = false;
        //             }
        //         }

        //         overall_time.stop();
        //         std::cout << "Contact::assemble: " << overall_time << std::endl;
        //         return has_contact_;
        //     }

        void Contact::transform(const Vector &in, Vector &out) {
            assert(impl_);
            out = transpose((*impl_->complete_transformation)) * in;
        }

        void Contact::inverse_transform(const Vector &in, Vector &out) {
            assert(impl_);
            out = (*impl_->complete_transformation) * in;
        }

        void Contact::transform(const Matrix &in, Matrix &out) {
            assert(impl_);

            const auto &T = (*impl_->complete_transformation);
            out = transpose(T) * in * T;
        }

        const Contact::Vector &Contact::gap() const {
            assert(impl_);
            return impl_->gap;
        }

        Contact::Vector &Contact::gap() {
            assert(impl_);
            return impl_->gap;
        }

        Contact::Contact() : params_(utopia::make_unique<Params>()) {}

        //     bool Contact::init_no_contact(const libMesh::MeshBase &mesh, const libMesh::DofMap
        //     &dof_map) {
        //         if (!contact_tensors_) {
        //             contact_tensors_ = std::make_shared<ConvertContactTensors>();
        //         }

        //         auto n_local_dofs = dof_map.n_local_dofs();

        //         contact_tensors_->gap = local_values(n_local_dofs, 100000000);
        //         contact_tensors_->weighted_gap = local_values(n_local_dofs, 100000000);
        //         contact_tensors_->normal = local_zeros(n_local_dofs);

        //         contact_tensors_->inv_mass_vector = local_values(n_local_dofs, 1.);
        //         contact_tensors_->is_contact = local_zeros(n_local_dofs);

        //         contact_tensors_->T = local_identity(n_local_dofs, n_local_dofs);
        //         contact_tensors_->orthogonal_trafo = local_identity(n_local_dofs, n_local_dofs);
        //         contact_tensors_->complete_transformation = local_identity(n_local_dofs, n_local_dofs);
        //         has_contact_ = false;
        //         has_glue_ = false;
        //         return true;
        //     }

        //     void Contact::remove_mass(const Vector &in, Vector &out) const {
        //         if (!empty(contact_tensors_->Q)) {
        //             out = contact_tensors_->Q * e_mul(contact_tensors_->inv_mass_vector, in);
        //             return;
        //         }

        //         out = e_mul(contact_tensors_->inv_mass_vector, in);
        //     }

        void Contact::read(Input &in) {
            params_->read(in);
            // in.get("black-list", [this](Input &in) {
            //     black_list_ = std::make_shared<ElementBlackList>(true);
            //     black_list_->read(in);
            // });
        }

    }  // namespace moonolith
}  // namespace utopia
