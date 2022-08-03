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

        void Contact::Params::read(Input &in) {
            std::set<int> temp;

            in.get("radius", this->search_radius);

            this->is_glue = std::make_shared<::moonolith::IsGlue>();

            in.get("pairs", [this, &temp](Input &array_node) {
                array_node.get_all([this, &temp](Input &node) {
                    int master = -1, slave = -1;
                    node.get("master", master);
                    node.get("slave", slave);

                    bool is_glued = false;
                    node.get("glue", is_glued);

                    assert(master != -1);
                    assert(slave != -1);
                    temp.insert(master);
                    temp.insert(slave);

                    this->contact_pair_tags.push_back({master, slave});
                    this->glued.push_back(is_glued);

                    if (is_glued) {
                        this->is_glue->insert(master, slave);
                    }
                });
            });

            in.get("search_radius", [this](Input &in) {
                in.get("default", this->search_radius);

                this->side_set_search_radius = std::make_shared<::moonolith::SearchRadius<double>>(this->search_radius);

                in.get("sides", [this](Input &array_node) {
                    array_node.get_all([this](Input &node) {
                        int id = -1;
                        double value = this->search_radius;

                        node.get("id", id);
                        node.get("value", value);

                        this->side_set_search_radius->insert(id, value);
                    });
                });
            });

            if (!this->side_set_search_radius) {
                this->side_set_search_radius = std::make_shared<::moonolith::SearchRadius<double>>(this->search_radius);
            }
        }

        void Contact::Params::describe(std::ostream &os) const {
            os << "search_radius: " << search_radius << "\n";
            os << "variable_number: " << variable_number << "\n";
            os << "use_biorthogonal_basis: " << use_biorthogonal_basis << "\n";
            os << "master, slave:\n";
            for (const auto &p : contact_pair_tags) {
                os << p.first << ", " << p.second << "\n";
            }

            if (side_set_search_radius) {
                os << "search-radius: " << std::endl;
                side_set_search_radius->describe(os);
            }

            os << std::endl;
        }

        class Contact::Impl {
        public:
            virtual ~Impl() = default;
            static constexpr Scalar LARGE_VALUE = 100000;

            virtual bool init(const FunctionSpace &, const Params &) = 0;

            bool init_no_contact(const FunctionSpace &space) {
                space.create_vector(gap);
                gap.set(LARGE_VALUE);

                space.create_vector(weighted_gap);
                weighted_gap.set(LARGE_VALUE);

                space.create_vector(normals);
                normals.set(0.);

                space.create_vector(inv_mass_vector);
                inv_mass_vector.set(1.);

                space.create_vector(is_contact);
                is_contact.set(0);

                complete_transformation = std::make_shared<Matrix>();
                orthogonal_transformation = std::make_shared<Matrix>();

                auto ml = square_matrix_layout(layout(gap));
                complete_transformation->identity(ml, 1);
                orthogonal_transformation->identity(ml, 1);

                has_contact = false;
                // has_glue_ = false;
                return true;
            }

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

                if (!empty(weighted_normals)) {
                    rename("wn", weighted_normals);
                    utopia::write("weighted_normals.m", weighted_normals);
                }

                if (!empty(normals)) {
                    rename("n", normals);
                    utopia::write("normals.m", normals);
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
            std::shared_ptr<Matrix> orthogonal_transformation, complete_transformation;

            Vector inv_mass_vector;
            Vector weighted_gap, gap;
            Vector weighted_normals, normals;
            Vector is_contact;
            Vector is_glue;

            bool has_contact{false};
        };

        template <int Dim>
        class Contact::ImplD : public Contact::Impl {
        public:
            using MoonolithMesh_t = ::moonolith::Mesh<Scalar, Dim>;
            using MoonolithSpace_t = ::moonolith::FunctionSpace<MoonolithMesh_t>;

            bool init(const FunctionSpace &space, const Params &params) override {
                std::shared_ptr<MoonolithSpace_t> space_d;
                if (space.mesh().manifold_dimension() == Dim - 1) {
                    space_d = space.raw_type<Dim>();
                } else {
                    assert(false && "IMPLEMENT ME");
                    return false;
                }

                MoonolithSpace_t elem_wise_space;
                space_d->separate_dofs(elem_wise_space);

                ::moonolith::ParContact<double, Dim> par_contact(space_d->mesh().comm(), Dim == 2);

                if (par_contact.assemble(
                        params.contact_pair_tags, elem_wise_space, params.side_set_search_radius, params.is_glue)) {
                    assert(false);  // TODO
                    // contact_data.finalize(par_contact.buffers, elem_wise_space, space);
                    return true;
                } else {
                    return false;
                }

                // if (impl_->has_contact) {
                //     // contact_tensors_ = contact_data.node_wise;
                //     assert(false);
                // } else {
                //     // init default
                //     init_no_contact(mesh, dof_map);
                // }

                // has_glue_ = impl_->has_contact && !params.is_glue->empty();

                // if (has_glue_) {
                //     double n_glued = sum(contact_tensors_->is_glue);

                //     if (n_glued == 0.) {
                //         has_glue_ = false;
                //     }
                // }

                assert(false);
                return false;
            }
            // void finalize(const moonolith::ContactBuffers<double> &buffers,
            //               const SpaceT<Dim> &element_wise_space,
            //               const SpaceT<Dim> &node_wise_space) {
            //     utopia::convert(buffers.B, element_wise.B);
            //     utopia::convert(buffers.D, element_wise.D);
            //     utopia::convert(buffers.Q, element_wise.Q);
            //     utopia::convert(buffers.basis_transform_inv, element_wise.basis_transform_inv);

            //     convert_tensor(buffers.gap, element_wise.weighted_gap);
            //     convert_tensor(buffers.normals, element_wise.weighted_normals);
            //     convert_tensor(buffers.is_glue, element_wise.is_glue);
            //     convert_tensor(buffers.is_contact, element_wise.is_contact);

            //     convert_to_node_wise(element_wise_space, element_wise, node_wise_space, *node_wise);
            // }

            //         static void normalize_and_build_orthgonal_trafo(const SpaceT<Dim> &node_wise_space,
            //                                                         ConvertContactTensors &node_wise) {
            //             node_wise.orthogonal_transformation =
            //                 local_sparse(node_wise_space.dof_map().n_local_dofs(),
            //                 node_wise_space.dof_map().n_local_dofs(), Dim);

            //             moonolith::HouseholderTransformation<double, Dim> H;
            //             auto &normals = node_wise.normals;

            //             auto r = range(normals);

            //             ReadAndWrite<Vector> rw_normals(normals);
            //             Read<Vector> r_is_c(node_wise.is_contact);
            //             Write<Matrix> w_ot(node_wise.orthogonal_transformation, utopia::LOCAL);

            //             moonolith::Vector<double, Dim> n;
            //             for (auto i = r.begin(); i < r.end(); i += Dim) {
            //                 const bool is_contact = node_wise.is_contact.get(i);

            //                 if (is_contact) {
            //                     for (int d = 0; d < Dim; ++d) {
            //                         n[d] = normals.get(i + d);
            //                     }

            //                     auto len = length(n);
            //                     assert(len > 0.0);

            //                     n /= len;

            //                     for (int d = 0; d < Dim; ++d) {
            //                         normals.set(i + d, n[d]);
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
            //                             node_wise.orthogonal_transformation.set(i + d1, i + d2, H(d1, d2));
            //                         }
            //                     }

            //                 } else {
            //                     for (int d1 = 0; d1 < Dim; ++d1) {
            //                         node_wise.orthogonal_transformation.set(i + d1, i + d1, 1.0);
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

            //             node_wise.weighted_normals = vector_perm * elem_wise.weighted_normals;

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
            //             node_wise.normals = node_wise.Q * e_mul(node_wise.inv_mass_vector,
            //             node_wise.weighted_normals);

            //             normalize_and_build_orthgonal_trafo(node_wise_space, node_wise);
            //             node_wise.complete_transformation = node_wise.T * node_wise.orthogonal_transformation;

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

        bool Contact::assemble(const FunctionSpace &space) {
            Chrono overall_time;
            overall_time.start();

            auto &&mesh = space.mesh();
            const int spatial_dim = mesh.spatial_dimension();

            if (spatial_dim == 2) {
                impl_ = utopia::make_unique<ImplD<2>>();
                impl_->init(space, *params_);
            } else if (spatial_dim == 3) {
                impl_ = utopia::make_unique<ImplD<3>>();
                impl_->init(space, *params_);
            }

            overall_time.stop();
            std::cout << "Contact::assemble: " << overall_time << std::endl;
            return impl_->has_contact;
        }

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

        Contact::Vector &Contact::normals() {
            assert(impl_);
            return impl_->normals;
        }

        Contact::Vector &Contact::is_contact() {
            assert(impl_);
            return impl_->is_contact;
        }

        std::shared_ptr<Contact::Matrix> Contact::mass_matrix() { return impl_->mass_matrix; }

        std::shared_ptr<Contact::Matrix> Contact::orthogonal_transformation() {
            return impl_->orthogonal_transformation;
        }

        std::shared_ptr<Contact::Matrix> Contact::complete_transformation() { return impl_->complete_transformation; }

        Contact::Contact() : params_(utopia::make_unique<Params>()) {}

        Contact::~Contact() {}

        void Contact::set_params(const Params &params) {
            assert(false);  // TODO
        }

        void Contact::set_banned_nodes(const std::shared_ptr<IndexArray> &banned_nodes) {
            assert(false);  // TODO
        }

        void Contact::read(Input &in) {
            params_->read(in);
            // in.get("black-list", [this](Input &in) {
            //     black_list_ = std::make_shared<ElementBlackList>(true);
            //     black_list_->read(in);
            // });
        }

        void Contact::describe(std::ostream &os) const {}

    }  // namespace moonolith
}  // namespace utopia
