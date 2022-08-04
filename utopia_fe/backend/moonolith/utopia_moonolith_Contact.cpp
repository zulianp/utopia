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
#include "moonolith_function_space.hpp"
#include "moonolith_matlab_scripter.hpp"
#include "moonolith_par_contact.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_sparse_matrix.hpp"

#include "utopia_moonolith_ConvertTensor.hpp"

// Integrate
// #include "utopia_moonolith_permutations.hpp"

// #include <memory>
// #include <vector>

namespace utopia {
    namespace moonolith {

        // template <int Dim>
        // void make_tensorize_permutation(const int dim,
        //                                 const MoonolithFunctionSpace<Dim> &from,
        //                                 const MoonolithFunctionSpace<Dim> &to,
        //                                 Matrix &mat) {
        //     using Scalar = Traits<Matrix>::Scalar;
        //     using SizeType = Traits<Matrix>::SizeType;

        //     std::vector<SizeType> irows(1), icols(1);
        //     std::vector<Scalar> vals(1, 1.0);

        //     auto max_nnz = from.dof_map().max_nnz();
        //     assert(max_nnz > 0);
        //     mat = local_sparse(to.dof_map().n_local_dofs() * dim, from.dof_map().n_local_dofs(), max_nnz);

        //     std::size_t n_elems = from.dof_map().n_elements();

        //     assert(n_elems == to.dof_map().n_elements());

        //     Write<Matrix> w(mat, utopia::GLOBAL_INSERT);

        //     for (std::size_t e = 0; e < n_elems; ++e) {
        //         const auto &from_dofs = from.dof_map().dofs(e);
        //         const auto &to_dofs = to.dof_map().dofs(e);

        //         const auto n_from = from_dofs.size();
        //         const auto n_to = to_dofs.size();

        //         assert(n_from == n_to);

        //         for (std::size_t i = 0; i < n_to; ++i) {
        //             for (int d = 0; d < dim; ++d) {
        //                 irows[0] = to_dofs[i] + d;
        //                 icols[0] = from_dofs[i];
        //                 mat.set_matrix(irows, icols, vals);
        //             }
        //         }
        //     }
        // }

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

            static void tensorize(const Matrix &T_x, const SizeType n_var, Matrix &T) {
                auto max_nnz = utopia::max_row_nnz(T_x);
                T.sparse(layout(T_x), max_nnz, max_nnz);

                assert(!empty(T));
                assert(T.row_range().extent() % n_var == 0);

                Write<Matrix> w(T);
                each_read(T_x, [&](const SizeType i, const SizeType j, const Scalar value) {
                    for (SizeType k = 0; k < n_var; ++k) {
                        T.set(i + k, j + k, value);
                    }
                });
            }

            static void tensorize(const SizeType n_var, Vector &t) {
                ReadAndWrite<Vector> w(t);
                auto r = range(t);

                assert(!empty(t));
                assert(r.extent() % n_var == 0);

                for (auto i = r.begin(); i < r.end(); i += n_var) {
                    const auto value = t.get(i);

                    for (SizeType k = 1; k < n_var; ++k) {
                        t.set(i + k, value);
                    }
                }
            }

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
            bool has_glue{false};
        };

        template <int Dim>
        class Contact::ImplD : public Contact::Impl {
        public:
            using MoonolithMesh_t = ::moonolith::Mesh<Scalar, Dim>;
            using MoonolithSpace_t = ::moonolith::FunctionSpace<MoonolithMesh_t>;

            void make_permutation(const MoonolithSpace_t &from, const MoonolithSpace_t &to, Matrix &mat) {
                // using Scalar = Traits<Matrix>::Scalar;
                // using SizeType = Traits<Matrix>::SizeType;

                std::vector<SizeType> irows(1), icols(1);
                std::vector<Scalar> vals(1, 1.0);

                // auto max_nnz = from.dof_map().max_nnz(); assert(max_nnz > 0);
                auto max_nnz = to.dof_map().max_nnz();
                assert(max_nnz > 0);

                Comm comm(from.mesh().comm().get());
                auto vl = layout(comm, to.dof_map().n_local_dofs(), from.dof_map().n_local_dofs());
                mat.sparse(square_matrix_layout(vl), max_nnz, max_nnz);

                std::size_t n_elems = from.dof_map().n_elements();

                assert(n_elems == to.dof_map().n_elements());

                Write<Matrix> w(mat, utopia::GLOBAL_INSERT);

                for (std::size_t e = 0; e < n_elems; ++e) {
                    const auto &from_dofs = from.dof_map().dofs(e);
                    const auto &to_dofs = to.dof_map().dofs(e);

                    const auto n_from = from_dofs.size();
                    const auto n_to = to_dofs.size();

                    assert(n_from == n_to);

                    for (std::size_t i = 0; i < n_to; ++i) {
                        irows[0] = to_dofs[i];
                        icols[0] = from_dofs[i];
                        mat.set_matrix(irows, icols, vals);
                    }
                }
            }

            void make_vector_permutation(const int dim,
                                         const MoonolithSpace_t &from,
                                         const MoonolithSpace_t &to,
                                         Matrix &mat) {
                // using Scalar = Traits<Matrix>::Scalar;
                // using SizeType = Traits<Matrix>::SizeType;

                std::vector<SizeType> irows(1), icols(1);
                std::vector<Scalar> vals(1, 1.0);

                // auto max_nnz =  from.dof_map().max_nnz(); assert(max_nnz > 0);
                auto max_nnz = to.dof_map().max_nnz();
                assert(max_nnz > 0);

                Comm comm(from.mesh().comm().get());
                auto vl = layout(comm, to.dof_map().n_local_dofs(), from.dof_map().n_local_dofs() * dim);
                mat.sparse(square_matrix_layout(vl), max_nnz, max_nnz);

                std::size_t n_elems = from.dof_map().n_elements();

                assert(n_elems == to.dof_map().n_elements());

                Write<Matrix> w(mat, utopia::GLOBAL_INSERT);

                for (std::size_t e = 0; e < n_elems; ++e) {
                    const auto &from_dofs = from.dof_map().dofs(e);
                    const auto &to_dofs = to.dof_map().dofs(e);

                    const auto n_from = from_dofs.size();
                    const auto n_to = to_dofs.size();

                    assert(n_from == n_to);

                    for (std::size_t i = 0; i < n_to; ++i) {
                        for (int d = 0; d < dim; ++d) {
                            irows[0] = to_dofs[i] + d;
                            icols[0] = from_dofs[i] * dim + d;
                            mat.set_matrix(irows, icols, vals);
                        }
                    }
                }
            }

            bool init(const FunctionSpace &space, const Params &params) override {
                std::shared_ptr<MoonolithSpace_t> space_d;
                if (space.mesh().manifold_dimension() == Dim - 1) {
                    space_d = space.raw_type<Dim>();
                } else {
                    assert(false && "IMPLEMENT ME");
                    return false;
                }

                MoonolithSpace_t element_wise_space;
                space_d->separate_dofs(element_wise_space);

                ::moonolith::ParContact<double, Dim> par_contact(space_d->mesh().comm(), Dim == 2);

                if (par_contact.assemble(
                        params.contact_pair_tags, element_wise_space, params.side_set_search_radius, params.is_glue)) {
                    this->has_contact = true;
                    finalize(par_contact.buffers, element_wise_space, *space_d);
                } else {
                    this->has_contact = false;
                    this->has_glue = false;
                }

                if (this->has_contact) {
                    has_glue = has_contact && !params.is_glue->empty();

                    if (has_glue) {
                        double n_glued = sum(is_glue);

                        if (n_glued == 0.) {
                            has_glue = false;
                        }
                    }

                } else {
                    // init default
                    this->init_no_contact(space);
                }

                return true;
            }

            void normalize_and_build_orthgonal_trafo(const MoonolithSpace_t &node_wise_space) {
                Comm comm(node_wise_space.mesh().comm().get());

                auto vl = layout(comm, node_wise_space.dof_map().n_local_dofs(), node_wise_space.dof_map().n_dofs());

                this->orthogonal_transformation->sparse(square_matrix_layout(vl), Dim, 0);

                ::moonolith::HouseholderTransformation<double, Dim> H;
                auto &normals = this->normals;

                auto r = range(normals);

                ReadAndWrite<Vector> rw_normals(normals);
                Read<Vector> r_is_c(this->is_contact);
                Write<Matrix> w_ot(*this->orthogonal_transformation, utopia::LOCAL);

                ::moonolith::Vector<double, Dim> n;
                for (auto i = r.begin(); i < r.end(); i += Dim) {
                    const bool is_contact = this->is_contact.get(i);

                    if (is_contact) {
                        for (int d = 0; d < Dim; ++d) {
                            n[d] = normals.get(i + d);
                        }

                        auto len = length(n);
                        assert(len > 0.0);

                        n /= len;

                        for (int d = 0; d < Dim; ++d) {
                            normals.set(i + d, n[d]);
                        }

                        n.x -= 1.0;

                        len = length(n);

                        if (approxeq(len, 0.0, 1e-15)) {
                            H.identity();
                        } else {
                            n /= len;
                            H.init(n);
                        }

                        assert(approxeq(std::abs(measure(H)), 1.0, 1e-8));

                        for (int d1 = 0; d1 < Dim; ++d1) {
                            for (int d2 = 0; d2 < Dim; ++d2) {
                                this->orthogonal_transformation->set(i + d1, i + d2, H(d1, d2));
                            }
                        }

                    } else {
                        for (int d1 = 0; d1 < Dim; ++d1) {
                            this->orthogonal_transformation->set(i + d1, i + d1, 1.0);
                        }
                    }
                }
            }

            void finalize(::moonolith::ContactBuffers<double> &buffers,
                          const MoonolithSpace_t &element_wise_space,
                          const MoonolithSpace_t &node_wise_space) {
                Matrix B, D, Q, Q_inv;
                Vector wg, g, wn, ig, ic;

                utopia::convert(buffers.B.get(), B);
                utopia::convert(buffers.D.get(), D);
                utopia::convert(buffers.Q.get(), Q);
                utopia::convert(buffers.Q_inv.get(), Q_inv);

                utopia::convert(buffers.gap.get(), wg);
                utopia::convert(buffers.normal.get(), wn);
                utopia::convert(buffers.is_glue.get(), ig);
                utopia::convert(buffers.is_contact.get(), ic);

                Matrix perm, vector_perm;
                make_permutation(element_wise_space, node_wise_space, perm);
                make_vector_permutation(Dim, element_wise_space, node_wise_space, vector_perm);

                this->weighted_gap = perm * wg;
                this->is_glue = perm * ig;
                this->is_contact = perm * ic;

                this->weighted_normals = vector_perm * wn;

                Matrix B_x = perm * B * transpose(perm);
                Matrix D_x = perm * D * transpose(perm);
                Matrix Q_x = perm * Q * transpose(perm);

                Matrix basis_transform_inv_x = perm * Q_inv * transpose(perm);

                normalize_rows(Q_x, 1e-15);
                normalize_rows(basis_transform_inv_x, 1e-15);

                // zeros and ones
                {
                    auto ic_view = local_view_device(this->is_contact);

                    parallel_for(
                        local_range_device(this->is_contact), UTOPIA_LAMBDA(const SizeType i) {
                            if (ic_view.get(i) > 0.0) {
                                ic_view.set(i, 1);
                            }
                        })
                }

                // zeros and ones
                {
                    auto ig_view = local_view_device(this->is_glue);

                    parallel_for(
                        local_range_device(this->is_glue), UTOPIA_LAMBDA(const SizeType i) {
                            if (ig_view.get(i) > 0.0) {
                                ig_view.set(i, 1);
                            }
                        })
                }

                this->inv_mass_vector = sum(D_x, 1);

                e_pseudo_inv(this->inv_mass_vector, this->inv_mass_vector, 1e-15);
                Matrix D_inv_x = diag(this->inv_mass_vector);

                Matrix T_temp_x = D_inv_x * B_x;
                Matrix T_x = Q_x * T_temp_x;

                this->tensorize(B_x, Dim, *this->coupling_matrix);
                this->tensorize(D_x, Dim, *this->mass_matrix);
                this->tensorize(T_x, Dim, *this->complete_transformation);
                this->tensorize(Q_x, Dim, *this->basis_transform);
                this->tensorize(Q_inv, Dim, *this->basis_transform_inv);

                this->tensorize(Dim, this->inv_mass_vector);
                this->tensorize(Dim, this->is_glue);

                normalize_rows(*this->complete_transformation, 1e-15);

                assert(check_op(*this->complete_transformation));

                this->complete_transformation->shift_diag(1);

                this->gap = (*this->basis_transform) * e_mul(this->inv_mass_vector, this->weighted_gap);
                this->normals = (*this->basis_transform) * e_mul(this->inv_mass_vector, this->weighted_normals);

                normalize_and_build_orthgonal_trafo(node_wise_space);
                (*this->complete_transformation) *= *this->orthogonal_transformation;

                {
                    auto g_view = local_view_device(this->gap);
                    auto ic_view = local_view_device(this->is_contact);
                    parallel_for(
                        local_range_device(this->gap), UTOPIA_LAMBDA(const SizeType i) {
                            if (ic_view.get(i) == 0.0) {
                                g_view.set(i, LARGE_VALUE);
                            }
                        });
                }
            }

            // static void convert_to_node_wise(const MoonolithSpace_t &element_wise_space,
            //                                  ConvertContactTensors &elem_wise,
            //                                  const MoonolithSpace_t &node_wise_space,
            //                                  ConvertContactTensors &node_wise)

            //         ConvertContactBuffers(MPI_Comm comm) { node_wise = std::make_shared<ConvertContactTensors>();
            //         }

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

            //             SpaceT element_wise_space;
            //             space.separate_dofs(element_wise_space);

            //             moonolith::ParContact<double, Dim> par_contact(comm, Dim == 2);

            //             if (par_contact.assemble(
            //                     params.contact_pair_tags, element_wise_space, params.side_set_search_radius,
            //                     params.is_glue)) {
            //                 contact_data.finalize(par_contact.buffers, element_wise_space, space);
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

        void Contact::set_params(const Params &params) { *params_ = params; }

        void Contact::set_banned_nodes(const std::shared_ptr<IndexArray> & /*banned_nodes*/) {
            assert(false);  // TODO
        }

        void Contact::read(Input &in) {
            params_->read(in);
            // in.get("black-list", [this](Input &in) {
            //     black_list_ = std::make_shared<ElementBlackList>(true);
            //     black_list_->read(in);
            // });
        }

        void Contact::describe(std::ostream &os) const { params_->describe(os); }

    }  // namespace moonolith
}  // namespace utopia
