#include "utopia_stk_FunctionSpace.hpp"

#include "utopia_Options.hpp"

#include "utopia_DirichletBoundary.hpp"
#include "utopia_FEVar.hpp"

#include "utopia_stk_Commons.hpp"
#include "utopia_stk_DofMap.hpp"
#include "utopia_stk_MeshIO.hpp"

// All stk includes
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <memory>
#include <unordered_set>

namespace utopia {
    namespace stk {

        class FunctionSpace::Impl {
        public:
            using Bucket_t = ::stk::mesh::Bucket;
            using BucketVector_t = ::stk::mesh::BucketVector;
            using Entity_t = ::stk::mesh::Entity;
            using Selector_t = ::stk::mesh::Selector;
            using IOBroker_t = ::stk::io::StkMeshIoBroker;
            using VectorField_t = ::stk::mesh::Field<Scalar, ::stk::mesh::Cartesian>;
            using MatrixField_t = ::stk::mesh::Field<Scalar, ::stk::mesh::Cartesian, ::stk::mesh::Cartesian>;
            using MatrixField3x3_t = ::stk::mesh::Field<Scalar, ::stk::mesh::Cartesian3d, ::stk::mesh::Cartesian3d>;

            std::string name{"main"};
            std::shared_ptr<Mesh> mesh;
            DirichletBoundary dirichlet_boundary;
            std::vector<FEVar> variables;
            std::shared_ptr<DofMap> dof_map;
            bool verbose{false};
            bool print_map{false};

            void node_eval(const BucketVector_t &node_buckets,
                           std::function<void(const SizeType idx, const Scalar *)> fun) {
                using Bucket_t = ::stk::mesh::Bucket;

                auto &meta_data = mesh->meta_data();
                auto &bulk_data = mesh->bulk_data();

                auto *coords = meta_data.coordinate_field();
                assert(coords);

                const int dim = this->mesh->spatial_dimension();

                auto &&local_to_global = dof_map->local_to_global();

                if (local_to_global.empty()) {
                    for (const auto &ib : node_buckets) {
                        const Bucket_t &b = *ib;
                        const Bucket_t::size_type length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            auto node = b[k];
                            auto moonolith_index = utopia::stk::convert_entity_to_index(node);

                            const Scalar *points = (const Scalar *)::stk::mesh::field_data(*coords, node);

                            auto idx = utopia::stk::convert_entity_to_index(node);
                            fun(idx, points);
                        }
                    }

                } else {
                    auto n_local_to_global = local_to_global.size();
                    for (const auto &ib : node_buckets) {
                        const Bucket_t &b = *ib;
                        const Bucket_t::size_type length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            auto node = b[k];
                            auto moonolith_index = utopia::stk::convert_entity_to_index(node);

                            const Scalar *points = (const Scalar *)::stk::mesh::field_data(*coords, node);

                            auto idx = utopia::stk::convert_entity_to_index(node);

                            assert(idx < n_local_to_global);

                            fun(local_to_global.global_node_id(idx), points);
                        }
                    }
                }
            }

            bool write_restart(const Path &path) {
                auto &meta_data = mesh->meta_data();
                auto &bulk_data = mesh->bulk_data();

                try {
                    Impl::IOBroker_t io_broker(mesh->comm().raw_comm());
                    io_broker.set_bulk_data(bulk_data);

                    auto out_id = io_broker.create_output_mesh(path.to_string(), ::stk::io::WRITE_RESTART);

                    for (auto &v : variables) {
                        ::stk::mesh::FieldBase *field = ::stk::mesh::get_field_by_name(v.name, meta_data);
                        io_broker.add_field(out_id, *field, v.name);
                    }

                    io_broker.process_output_request(out_id, 0);
                    return true;

                } catch (const std::exception &ex) {
                    utopia::err() << "Mesh::write(\"" << path.to_string() << "\") error: " << ex.what() << '\n';
                    assert(false);
                    return false;
                }
            }

            bool check_variables() {
                auto &meta_data = mesh->meta_data();

                for (auto &v : variables) {
                    ::stk::mesh::FieldBase *field = ::stk::mesh::get_field_by_name(v.name, meta_data);
                    assert(field && "must be defined");
                    if (!field) {
                        utopia::err() << "Variable with " << v.name << " does not exist!\n";
                        return false;
                    }
                }

                return true;
            }

            void register_variables() {
                auto &meta_data = mesh->meta_data();
                meta_data.enable_late_fields();

                for (auto &v : variables) {
                    register_variable(v);
                }
            }

            void register_variable(const FEVar &v) {
                auto &meta_data = mesh->meta_data();
                auto &&part = meta_data.universal_part();

                ::stk::mesh::FieldBase *has_field = ::stk::mesh::get_field_by_name(v.name, meta_data);

                if (!has_field) {
                    if (v.n_components == 1) {
                        auto &field =
                            meta_data.declare_field<::stk::mesh::Field<Scalar>>(::stk::topology::NODE_RANK, v.name, 1);
                        ::stk::mesh::put_field_on_mesh(field, part, 1, nullptr);
                    } else {
                        auto &field =
                            meta_data.declare_field<Impl::VectorField_t>(::stk::topology::NODE_RANK, v.name, 1);
                        ::stk::mesh::put_field_on_mesh(field, part, v.n_components, nullptr);
                    }
                } else {
                    // TODO check if definition is correct
                }
            }

            void nodal_field_to_local_vector(Vector &v) { nodal_field_to_local_vector(this->variables, v); }

            void nodal_field_to_local_vector(const std::vector<FEVar> &variables, Vector &v) {
                auto &meta_data = mesh->meta_data();
                auto &bulk_data = mesh->bulk_data();
                auto &node_buckets = utopia::stk::universal_nodes(bulk_data);
                // const int n_var = dof_map->n_var();

                int n_var = 0;

                for (auto &v : variables) {
                    n_var += v.n_components;
                }

                auto v_view = local_view_device(v);

                int offset = 0;
                for (auto &v : variables) {
                    ::stk::mesh::FieldBase *field = ::stk::mesh::get_field_by_name(v.name, meta_data);
                    assert(field);

                    if (!field) {
                        Utopia::Abort("utopia::stk::FunctionSpace: is trying to access undefiend field");
                    }

                    for (const auto &ib : node_buckets) {
                        const Bucket_t &b = *ib;
                        const Bucket_t::size_type length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            Entity_t node = b[k];
                            auto idx = utopia::stk::convert_entity_to_index(node);

                            if (bulk_data.in_receive_ghost(node)) {
                                idx = dof_map->shift_aura_idx(idx);
                            }

                            const Scalar *values = (const Scalar *)::stk::mesh::field_data(*field, node);
                            const int n_comp = ::stk::mesh::field_scalars_per_entity(*field, node);

                            for (int d = 0; d < n_comp; ++d) {
                                assert(offset + d < n_var);
                                v_view.set(idx * n_var + offset + d, values[d]);
                            }
                        }
                    }

                    offset += v.n_components;
                }
            }

            void local_vector_to_nodal_field(const Vector &v) {
                auto &meta_data = mesh->meta_data();
                auto &bulk_data = mesh->bulk_data();
                auto &node_buckets = utopia::stk::universal_nodes(bulk_data);
                const int n_var = dof_map->n_var();

                auto v_view = local_view_device(v);

                int offset = 0;
                for (auto &v : variables) {
                    ::stk::mesh::FieldBase *field = ::stk::mesh::get_field_by_name(v.name, meta_data);

                    for (const auto &ib : node_buckets) {
                        const Bucket_t &b = *ib;
                        const Bucket_t::size_type length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            Entity_t node = b[k];
                            auto idx = utopia::stk::convert_entity_to_index(node);

                            if (bulk_data.in_receive_ghost(node)) {
                                idx = dof_map->shift_aura_idx(idx);
                            }

                            Scalar *values = (Scalar *)::stk::mesh::field_data(*field, node);
                            int n_comp = ::stk::mesh::field_scalars_per_entity(*field, node);

                            for (int d = 0; d < n_comp; ++d) {
                                assert(offset + d < n_var);
                                values[d] = v_view.get(idx * n_var + offset + d);
                            }
                        }
                    }

                    offset += v.n_components;
                }
            }

            void copy_local_vector_to_nodal_field(const Vector &v, ::stk::mesh::FieldBase &field) const {
                auto &bulk_data = mesh->bulk_data();
                auto &node_buckets = utopia::stk::universal_nodes(bulk_data);

                auto v_view = local_view_device(v);

                for (const auto &ib : node_buckets) {
                    const Bucket_t &b = *ib;
                    const Bucket_t::size_type length = b.size();

                    for (Bucket_t::size_type k = 0; k < length; ++k) {
                        Entity_t node = b[k];
                        auto idx = utopia::stk::convert_entity_to_index(node);

                        if (bulk_data.in_receive_ghost(node)) {
                            idx = dof_map->shift_aura_idx(idx);
                        }

                        Scalar *points = (Scalar *)::stk::mesh::field_data(field, node);
                        int n_comp = ::stk::mesh::field_scalars_per_entity(field, node);

                        for (int d = 0; d < n_comp; ++d) {
                            points[d] = v_view.get(idx * n_comp + d);
                        }
                    }
                }
            }

            void local_vector_to_nodal_field(const std::string &name, const Vector &v) {
                auto &meta_data = mesh->meta_data();
                ::stk::mesh::FieldBase *field = ::stk::mesh::get_field_by_name(name, meta_data);
                assert(field);
                copy_local_vector_to_nodal_field(v, *field);
            }

            void read_meta(Input &in) {
                int n_var = 0;

                if (!Options()
                         .add_option("name", name, "Unique name of the function space")
                         .add_option(
                             "n_var", n_var, "Number of variables per node (instead of specifiying all of them).")
                         .add_option("boundary_conditions", dirichlet_boundary, "Boundary conditions.")
                         .add_option("verbose", verbose, "Verbose output.")
                         .add_option("print_map", print_map, "Print local to global mapping.")
                         .parse(in)) {
                    return;
                }

                in.get("mesh", [this](Input &mesh_node) {
                    std::string type;
                    mesh_node.get("type", type);

                    if (type == "cube") {
                        dirichlet_boundary.convert_user_space_names(SideSet::Cube());
                    }
                });

                in.get("variables", [this](Input &in) {
                    in.get_all([&](Input &in) {
                        FEVar v;
                        v.read(in);
                        variables.push_back(v);
                    });
                });

                if (variables.empty()) {
                    FEVar v;
                    // At least one variable
                    v.n_components = std::max(1, n_var);
                    variables.push_back(v);
                }

                // int counted_vars = 0;

                // for (auto &v : variables) {
                //     counted_vars += v.n_components;
                // }

                // assert(n_var == 0 || n_var == counted_vars);

                // n_var = counted_vars;
                // dof_map->set_n_var(n_var);

                count_and_set_variables();
            }

            void count_and_set_variables() { dof_map->set_n_var(count_variables()); }

            int count_variables() {
                int counted_vars = 0;

                for (auto &v : variables) {
                    counted_vars += v.n_components;
                }

                return counted_vars;
            }

            template <class Fun>
            void apply_varying_constraints(const BucketVector_t &buckets,
                                           Fun fun,
                                           const int component,
                                           Vector &v) const {
                using Bucket_t = ::stk::mesh::Bucket;

                auto &meta_data = mesh->meta_data();
                // auto &bulk_data = mesh().bulk_data();

                auto &&local_to_global = dof_map->local_to_global();

                const int nv = dof_map->n_var();
                const int spatial_dim = mesh->spatial_dimension();

                auto *coords = meta_data.coordinate_field();
                assert(coords);

                if (local_to_global.empty()) {
                    auto v_view = local_view_device(v);

                    for (auto *b_ptr : buckets) {
                        auto &b = *b_ptr;
                        const Bucket_t::size_type length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            auto node = b[k];
                            auto idx = utopia::stk::convert_entity_to_index(node);

                            Scalar p[3] = {0., 0., 0.};
                            const Scalar *points = (const Scalar *)::stk::mesh::field_data(*coords, node);

                            for (int d = 0; d < spatial_dim; d++) {
                                p[d] = points[d];
                            }

                            const Scalar value = fun(p[0], p[1], p[2]);
                            v_view.set(idx * nv + component, value);
                        }
                    }
                } else {
                    Write<Vector> w(v, utopia::GLOBAL_INSERT);

                    for (auto *b_ptr : buckets) {
                        auto &b = *b_ptr;
                        const Bucket_t::size_type length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            auto node = b[k];
                            auto local_idx = utopia::stk::convert_entity_to_index(node);
                            assert(local_idx < local_to_global.size());
                            Scalar p[3] = {0., 0., 0.};

                            const Scalar *points = (const Scalar *)::stk::mesh::field_data(*coords, node);

                            for (int d = 0; d < spatial_dim; d++) {
                                p[d] = points[d];
                            }

                            const Scalar value = fun(p[0], p[1], p[2]);
                            v.c_set(local_to_global(local_idx, component), value);
                        }
                    }
                }
            }

            void apply_uniform_constraints(const BucketVector_t &buckets,
                                           const Scalar value,
                                           const int component,
                                           Vector &v) const {
                using Bucket_t = ::stk::mesh::Bucket;

                auto &meta_data = mesh->meta_data();
                // auto &bulk_data = mesh().bulk_data();

                auto &&local_to_global = dof_map->local_to_global();

                const int nv = dof_map->n_var();
                const int spatial_dim = mesh->spatial_dimension();

                if (local_to_global.empty()) {
                    auto v_view = local_view_device(v);

                    for (auto *b_ptr : buckets) {
                        auto &b = *b_ptr;
                        const Bucket_t::size_type length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            auto node = b[k];
                            auto idx = utopia::stk::convert_entity_to_index(node);
                            v_view.set(idx * nv + component, value);
                        }
                    }
                } else {
                    Write<Vector> w(v, utopia::GLOBAL_INSERT);

                    for (auto *b_ptr : buckets) {
                        auto &b = *b_ptr;
                        const Bucket_t::size_type length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            auto node = b[k];
                            auto local_idx = utopia::stk::convert_entity_to_index(node);
                            assert(local_idx < local_to_global.size());
                            v.c_set(local_to_global(local_idx, component), value);
                        }
                    }
                }
            }
        };

        void FunctionSpace::create_node_to_element_matrix(Matrix &matrix) const {
            using Bucket_t = ::stk::mesh::Bucket;
            using BucketVector_t = ::stk::mesh::BucketVector;
            using Entity_t = ::stk::mesh::Entity;
            using Scalar_t = Traits<utopia::stk::Mesh>::Scalar;
            using Size_t = Traits<utopia::stk::Mesh>::SizeType;
            using MetaData_t = ::stk::mesh::MetaData;
            using BulkData_t = ::stk::mesh::BulkData;

            //
            if (impl_->dof_map->empty()) {
                impl_->dof_map->init(*this->mesh_ptr(), impl_->print_map);
            }

            auto &bulk_data = this->mesh().bulk_data();
            SizeType num_nodes_x_element = 0;
            {
                auto &elem_buckets = utopia::stk::local_elements(bulk_data);

                {
                    for (const auto &ib : elem_buckets) {
                        const Bucket_t &b = *ib;
                        const Bucket_t::size_type length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            Entity_t elem = b[k];
                            num_nodes_x_element = bulk_data.num_nodes(elem);
                            break;
                        }

                        if (length) {
                            break;
                        }
                    }
                }
            }

            auto nodal_layout = layout(comm(), n_local_dofs(), n_dofs());
            auto elemental_layout = layout(comm(),
                                           this->mesh().n_local_elements() * num_nodes_x_element,
                                           this->mesh().n_elements() * num_nodes_x_element);

            auto ml = matrix_layout(elemental_layout, nodal_layout);

            utopia::out() << "create_node_to_element_matrix\n";
            utopia::out() << "nle: " << this->mesh().n_local_elements() << " nge: " << this->mesh().n_elements()
                          << "\n";
            utopia::out() << ml.size(0) << ", " << ml.size(1) << "\n";

            matrix.sparse(ml, num_nodes_x_element, num_nodes_x_element);

            Write<PetscMatrix> w(matrix, utopia::GLOBAL_ADD);

            const SizeType n_dofs = this->n_dofs();
            // const SizeType nn = n_dofs / this->n_var();

            auto &&local_to_global = this->dof_map().local_to_global();

            IndexArray row_idx(num_nodes_x_element, -1);
            IndexArray col_idx(num_nodes_x_element, -1);
            ScalarArray val(num_nodes_x_element * num_nodes_x_element, 1.0);

            const bool is_block = matrix.is_block();

            auto rr = matrix.row_range();

            auto &meta_data = this->mesh().meta_data();
            const BucketVector_t &elem_buckets =
                bulk_data.get_buckets(::stk::topology::ELEMENT_RANK, meta_data.locally_owned_part());

            Size_t elem_idx = 0;
            for (const auto &ib : elem_buckets) {
                const Bucket_t &b = *ib;
                const Bucket_t::size_type length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    // get the current node entity and extract the id to fill it into the field
                    Entity_t elem = b[k];
                    // const Size_t elem_idx = utopia::stk::convert_entity_to_index(elem) - n_local_nodes;

                    // auto g_eid = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(elem));

                    for (int i = 0; i < num_nodes_x_element; i++) {
                        row_idx[i] = rr.begin() + elem_idx * num_nodes_x_element + i;
                    }

                    const Size_t n_nodes = bulk_data.num_nodes(elem);
                    UTOPIA_UNUSED(n_nodes);

                    assert(num_nodes_x_element == n_nodes);

                    auto node_ids = bulk_data.begin_nodes(elem);

                    if (local_to_global.empty()) {
                        for (Size_t i = 0; i < num_nodes_x_element; ++i) {
                            col_idx[i] = utopia::stk::convert_entity_to_index(node_ids[i]);
                            assert(col_idx[i] < this->n_dofs());
                            assert(col_idx[i] >= 0);
                        }
                    } else {
                        for (Size_t i = 0; i < num_nodes_x_element; ++i) {
                            col_idx[i] = local_to_global.block(utopia::stk::convert_entity_to_index(node_ids[i]));
                            assert(col_idx[i] < this->n_dofs());
                            assert(col_idx[i] >= 0);
                        }
                    }

                    matrix.set_matrix(row_idx, col_idx, val);

                    ++elem_idx;
                }
            }
        }

        FunctionSpace::FunctionSpace(const Comm &comm) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = std::make_shared<Mesh>(comm);
            impl_->dof_map = std::make_shared<DofMap>();
        }

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = mesh;
            impl_->dof_map = std::make_shared<DofMap>();
        }

        void FunctionSpace::copy_meta_info_from(const FunctionSpace &other) {
            impl_->variables = other.impl_->variables;
            dof_map().set_n_var(other.n_var());
        }

        void FunctionSpace::initialize(const bool valid_local_id_mode) {
            impl_->register_variables();

            // if (this->mesh().has_aura()) {
            //     this->mesh().create_edges();
            // }

            // impl_->dof_map->init(*this->mesh_ptr(), impl_->print_map);
            impl_->dof_map->set_valid_local_id_mode(valid_local_id_mode);
            impl_->dof_map->init(*this->mesh_ptr(), impl_->print_map);

            if (impl_->verbose) {
                std::stringstream ss;
                describe(ss);
                comm().synched_print(ss.str());
            }
        }

        int FunctionSpace::add_variable(const FEVar &var) {
            int num_var = impl_->variables.size();
            impl_->variables.push_back(var);
            impl_->register_variable(var);
            impl_->count_and_set_variables();
            return num_var;
        }

        FunctionSpace::~FunctionSpace() = default;

        bool FunctionSpace::write(const Path &path, const Vector &x) {
            if (comm().size() == 1) {
                impl_->local_vector_to_nodal_field(x);
            } else {
                Vector x_local;
                global_to_local(x, x_local);
                impl_->local_vector_to_nodal_field(x_local);
            }

            return impl_->write_restart(path);
        }

        void FunctionSpace::init(const std::shared_ptr<Mesh> &mesh) { impl_->mesh = mesh; }

        void FunctionSpace::update(const SimulationTime<Scalar> &time) { impl_->dirichlet_boundary.update(time); }

        void FunctionSpace::read_meta(Input &in) { impl_->read_meta(in); }

        void FunctionSpace::register_variables() { impl_->register_variables(); }

        void FunctionSpace::read(Input &in) {
            in.get("mesh", *impl_->mesh);

            if (mesh().empty()) return;

            read_meta(in);

            impl_->register_variables();

            // if (this->mesh().has_aura()) {
            //     this->mesh().create_edges();
            // }

            impl_->dof_map->init(*this->mesh_ptr(), impl_->print_map);

            if (impl_->verbose) {
                std::stringstream ss;
                describe(ss);
                comm().synched_print(ss.str());
            }
        }

        bool FunctionSpace::read_with_fields(Input &in, std::vector<std::shared_ptr<Field<FunctionSpace>>> &val) {
            MeshIO io(*impl_->mesh);
            io.import_all_field_data(true);
            in.get("mesh", io);

            if (!io.load()) {
                return false;
            }

            impl_->read_meta(in);
            impl_->register_variables();
            impl_->dof_map->init(*this->mesh_ptr(), impl_->print_map);

            in.get("fields", [&](Input &array_node) {
                array_node.get_all([&](Input &node) {
                    FEVar var;
                    var.read(node);

                    SizeType n_nodes = utopia::stk::count_universal_nodes(mesh().bulk_data());
                    SizeType nn = n_nodes * var.n_components;
                    Vector lv(layout(Comm::self(), nn, nn), 0.0);
                    lv.set_block_size(var.n_components);

                    impl_->nodal_field_to_local_vector({var}, lv);

                    auto gv = std::make_shared<Vector>();

                    gv->zeros(
                        layout(comm(), mesh().n_local_nodes() * var.n_components, mesh().n_nodes() * var.n_components));
                    gv->set_block_size(var.n_components);

                    local_to_global(lv, *gv, OVERWRITE_MODE);

                    auto field = std::make_shared<Field<FunctionSpace>>();
                    field->set_data(gv);
                    field->set_space(make_ref(*this));

                    field->set_name(var.name);

                    val.push_back(field);
                });
            });

            return true;
        }

        bool FunctionSpace::read_with_state(Input &in, Field<FunctionSpace> &field) {
            MeshIO io(*impl_->mesh);
            io.import_all_field_data(true);
            in.get("mesh", io);

            if (!io.load()) {
                return false;
            }

            impl_->read_meta(in);
            impl_->register_variables();
            impl_->dof_map->init(*this->mesh_ptr(), impl_->print_map);

            std::vector<FEVar> fields;
            in.get("fields", [&](Input &array_node) {
                array_node.get_all([&](Input &node) {
                    FEVar var;
                    var.read(node);
                    fields.push_back(var);
                });
            });

            if (fields.empty()) {
                if (!impl_->check_variables()) {
                    return false;
                }

                Vector lv;
                create_local_vector(lv);
                impl_->nodal_field_to_local_vector(lv);

                auto gv = std::make_shared<Vector>();
                create_vector(*gv);
                local_to_global(lv, *gv, OVERWRITE_MODE);
                field.set_data(gv);
                field.set_space(make_ref(*this));

                assert(impl_->variables.size() == 1);

                if (!impl_->variables.empty()) {
                    field.set_name(impl_->variables[0].name);
                    rename(impl_->variables[0].name, *gv);
                    field.set_tensor_size(impl_->variables[0].n_components);
                }

            } else {
                int n_var = 0;

                for (auto &v : fields) {
                    n_var += v.n_components;
                }

                SizeType n_nodes = utopia::stk::count_universal_nodes(mesh().bulk_data());
                SizeType nn = n_nodes * n_var;
                Vector lv(layout(Comm::self(), nn, nn), 0.0);
                lv.set_block_size(n_var);

                impl_->nodal_field_to_local_vector(fields, lv);

                auto gv = std::make_shared<Vector>();

                gv->zeros(layout(comm(), mesh().n_local_nodes() * n_var, mesh().n_nodes() * n_var));
                gv->set_block_size(n_var);

                local_to_global(lv, *gv, OVERWRITE_MODE);
                field.set_data(gv);
                field.set_space(make_ref(*this));

                field.set_name(fields[0].name);
            }

            return true;
        }

        void FunctionSpace::describe(std::ostream &os) const {
            mesh().describe(os);

            if (comm().rank() == 0) {
                impl_->dirichlet_boundary.describe(os);
            }

            os << '\n';
            os << "n_vars: " << n_var() << '\n';
            os << "n_local_dofs: " << n_local_dofs() << '\n';
        }

        std::shared_ptr<Mesh> FunctionSpace::mesh_ptr() const { return impl_->mesh; }

        const Mesh &FunctionSpace::mesh() const {
            assert(impl_->mesh);
            return *impl_->mesh;
        }

        Mesh &FunctionSpace::mesh() {
            assert(impl_->mesh);
            return *impl_->mesh;
        }

        FunctionSpace::SizeType FunctionSpace::n_dofs() const { return impl_->mesh->n_nodes() * n_var(); }
        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const { return impl_->mesh->n_local_nodes() * n_var(); }

        int FunctionSpace::n_var() const { return dof_map().n_var(); }

        void FunctionSpace::set_n_var(const int n_var) { dof_map().set_n_var(n_var); }

        void FunctionSpace::create_vector(Vector &v) const {
            v.zeros(layout(comm(), n_local_dofs(), n_dofs()));
            v.set_block_size(n_var());
        }

        void FunctionSpace::create_field(Field<FunctionSpace> &field) {
            auto gv = std::make_shared<Vector>();
            create_vector(*gv);
            field.set_data(gv);
            field.set_space(make_ref(*this));

            assert(impl_->variables.size() == 1);

            if (!impl_->variables.empty()) {
                field.set_name(impl_->variables[0].name);
                rename(impl_->variables[0].name, *gv);
                field.set_tensor_size(impl_->variables[0].n_components);
            }
        }

        void FunctionSpace::create_nodal_vector_field(const int vector_size, Field<FunctionSpace> &field) {
            auto gv = std::make_shared<Vector>();
            gv->zeros(layout(comm(), mesh().n_local_nodes() * vector_size, mesh().n_nodes() * vector_size));
            gv->set_block_size(vector_size);

            field.set_data(gv);
            field.set_space(make_ref(*this));

            // if (!impl_->variables.empty()) {
            field.set_tensor_size(vector_size);
            // }
        }

        void FunctionSpace::create_local_vector(Vector &v) const {
            SizeType nn = utopia::stk::count_universal_nodes(mesh().bulk_data()) * n_var();
            v.zeros(layout(Comm::self(), nn, nn));
            v.set_block_size(n_var());
        }

        void FunctionSpace::create_matrix(Matrix &m) const {
            if (impl_->dof_map->empty()) {
                impl_->dof_map->init(*this->mesh_ptr(), impl_->print_map);
            }

            auto vl = layout(comm(), n_local_dofs(), n_dofs());
            auto ml = square_matrix_layout(vl);

            // https://petsc.org/release/src/mat/impls/is/matis.c.html#MatSetValues_IS
            // For NON overlapping DD (with petsc)
            // SETUP: set local 2 global after MatSetType(MATIS)
            // 1) ISCreateGeneral ...
            // 2) ISLocalToGlobalMappingCreateIS
            // 3) MatSetLocalToGlobalMapping

            // MatSetValues does to job even for this matrix type

            if (this->n_var() == 1) {
                m.sparse(ml, impl_->dof_map->d_nnz(), impl_->dof_map->o_nnz());
            } else {
                m.block_sparse(ml, impl_->dof_map->d_nnz(), impl_->dof_map->o_nnz(), this->n_var());
            }
        }

        const FunctionSpace::DirichletBoundary &FunctionSpace::dirichlet_boundary() const {
            return impl_->dirichlet_boundary;
        }

        void FunctionSpace::apply_constraints(Matrix &m, const Scalar diag_value) const {
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            using IndexSet = Traits<Matrix>::IndexSet;

            const SizeType nl_dofs = n_local_dofs();
            IndexSet constrains;
            constrains.reserve(nl_dofs);

            const int nv = n_var();

            auto &&local_to_global = dof_map().local_to_global();

            if (local_to_global.empty()) {
                for (auto &bc_ptr : impl_->dirichlet_boundary) {
                    auto &bc = *bc_ptr;
                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
                        auto &buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, *part);

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                // auto idx = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node));
                                auto idx = utopia::stk::convert_entity_to_index(node);

                                const SizeType dof_idx = idx * nv + bc.component;
                                assert(dof_idx < nl_dofs);
                                constrains.push_back(dof_idx);
                            }
                        }
                    }
                }
            } else {
                // auto rr = row_range(m);

                for (auto &bc_ptr : impl_->dirichlet_boundary) {
                    auto &bc = *bc_ptr;

                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
                        auto &buckets =
                            bulk_data.get_buckets(::stk::topology::NODE_RANK, *part & meta_data.locally_owned_part());

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                auto local_idx = utopia::stk::convert_entity_to_index(node);
                                assert(local_idx < local_to_global.size());

                                // if (rr.inside(idx)) {
                                constrains.push_back(local_to_global(local_idx, bc.component));
                                // }
                            }
                        }
                    }
                }
            }

            set_zero_rows(m, constrains, diag_value);
        }

        void FunctionSpace::overwrite_parts(const std::vector<std::string> &parts,
                                            const std::vector<int> &components,
                                            const Vector &source,
                                            Vector &destination) const {
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            auto &&local_to_global = dof_map().local_to_global();

            const int nv = n_var();

            if (local_to_global.empty()) {
                auto source_view = local_view_device(source);
                auto destination_view = local_view_device(destination);

                for (auto &part_name : parts) {
                    auto *part = meta_data.get_part(part_name);
                    if (part) {
                        auto &buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, *part);

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                auto idx = utopia::stk::convert_entity_to_index(node);

                                for (int c : components) {
                                    auto k = idx * nv + c;
                                    destination_view.set(k, source_view.get(k));
                                }
                            }
                        }
                    }
                }
            } else {
                Write<Vector> w(destination, utopia::GLOBAL_INSERT);
                Read<Vector> r(source);

                for (auto &part_name : parts) {
                    auto *part = meta_data.get_part(part_name);
                    if (part) {
                        auto &buckets =
                            bulk_data.get_buckets(::stk::topology::NODE_RANK, *part & meta_data.locally_owned_part());

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                auto local_idx = utopia::stk::convert_entity_to_index(node);
                                assert(local_idx < local_to_global.size());

                                for (int c : components) {
                                    auto k = local_to_global(local_idx, c);
                                    destination.c_set(k, source.get(k));
                                }
                            }
                        }
                    }
                }
            }
        }

        void FunctionSpace::apply_constraints(Vector &v) const {
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            auto &&local_to_global = dof_map().local_to_global();

            const int nv = n_var();
            const int spatial_dim = mesh().spatial_dimension();

            auto *coords = meta_data.coordinate_field();
            assert(coords);

            if (local_to_global.empty()) {
                assert(comm().size() == 1);

                auto v_view = local_view_device(v);

                for (auto &bc_ptr : impl_->dirichlet_boundary) {
                    auto &bc = *bc_ptr;
                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
#ifdef UTOPIA_WITH_TINY_EXPR
                        DirichletBoundary::VaryingCondition *varying_bc = nullptr;
#endif
                        double value = 0;
                        const bool is_uniform = bc.is_uniform();

                        if (is_uniform) {
                            value = static_cast<DirichletBoundary::UniformCondition &>(bc).value();
                        } else
#ifdef UTOPIA_WITH_TINY_EXPR
                            if (!(varying_bc = dynamic_cast<DirichletBoundary::VaryingCondition *>(&bc)))
#endif
                        {
                            continue;
                        }

                        auto &buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, *part);

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                auto idx = utopia::stk::convert_entity_to_index(node);

#ifdef UTOPIA_WITH_TINY_EXPR
                                if (varying_bc) {
                                    Scalar p[3] = {0., 0., 0.};
                                    const Scalar *points = (const Scalar *)::stk::mesh::field_data(*coords, node);

                                    for (int d = 0; d < spatial_dim; d++) {
                                        p[d] = points[d];
                                    }
                                    value = varying_bc->eval(p[0], p[1], p[2]);
                                }
#else
                                Utopia::Abort("Varying boundary conditions require UTOPIA_WITH_TINY_EXPR=ON!");
#endif

                                v_view.set(idx * nv + bc.component, value);
                            }
                        }
                    }
                }
            } else {
                Write<Vector> w(v, utopia::GLOBAL_INSERT);

                for (auto &bc_ptr : impl_->dirichlet_boundary) {
                    auto &bc = *bc_ptr;
                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
#ifdef UTOPIA_WITH_TINY_EXPR
                        DirichletBoundary::VaryingCondition *varying_bc = nullptr;
#endif
                        double value = 0;
                        const bool is_uniform = bc.is_uniform();
                        if (is_uniform) {
                            value = static_cast<DirichletBoundary::UniformCondition &>(bc).value();
                        } else
#ifdef UTOPIA_WITH_TINY_EXPR
                            if (!(varying_bc = dynamic_cast<DirichletBoundary::VaryingCondition *>(&bc)))
#endif
                        {
                            continue;
                        }

                        auto &buckets =
                            bulk_data.get_buckets(::stk::topology::NODE_RANK, *part & meta_data.locally_owned_part());

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                auto local_idx = utopia::stk::convert_entity_to_index(node);
                                assert(local_idx < local_to_global.size());
#ifdef UTOPIA_WITH_TINY_EXPR
                                if (varying_bc) {
                                    Scalar p[3] = {0., 0., 0.};

                                    const Scalar *points = (const Scalar *)::stk::mesh::field_data(*coords, node);

                                    for (int d = 0; d < spatial_dim; d++) {
                                        p[d] = points[d];
                                    }

                                    value = varying_bc->eval(p[0], p[1], p[2]);
                                }
#else
                                Utopia::Abort("Varying boundary conditions require UTOPIA_WITH_TINY_EXPR=ON!");
#endif

                                v.c_set(local_to_global(local_idx, bc.component), value);
                            }
                        }
                    }
                }
            }

            for (auto &bc_ptr : impl_->dirichlet_boundary) {
                if (bc_ptr->is_uniform()) continue;

                auto ow_bc = std::dynamic_pointer_cast<DirichletBoundary::OverwriteCondition>(bc_ptr);

                if (ow_bc && ow_bc->vector) {
                    overwrite_parts({ow_bc->name}, {ow_bc->component}, *ow_bc->vector, v);
                }
            }
        }

        void FunctionSpace::apply_constraints_time_derivative(Vector &v) const {
            // Default is zero
            apply_zero_constraints(v);

#ifdef UTOPIA_WITH_TINY_EXPR
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();
            auto &&local_to_global = dof_map().local_to_global();

            for (auto &bc_ptr : impl_->dirichlet_boundary) {
                auto &bc = *bc_ptr;
                if (!bc.has_time_derivative()) continue;
                auto *part = meta_data.get_part(bc.name);
                if (!part) continue;

                DirichletBoundary::UniformCondition *uniform_bc = nullptr;
                DirichletBoundary::VaryingCondition *varying_bc = nullptr;

                if ((uniform_bc = dynamic_cast<DirichletBoundary::UniformCondition *>(&bc))) {
                    double val = uniform_bc->time_derivative();

                    if (local_to_global.empty()) {
                        auto &buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, *part);
                        impl_->apply_uniform_constraints(buckets, val, bc.component, v);
                    } else {
                        auto &buckets =
                            bulk_data.get_buckets(::stk::topology::NODE_RANK, *part & meta_data.locally_owned_part());
                        impl_->apply_uniform_constraints(buckets, val, bc.component, v);
                    }

                } else if ((varying_bc = dynamic_cast<DirichletBoundary::VaryingCondition *>(&bc))) {
                    if (local_to_global.empty()) {
                        auto &buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, *part);
                        impl_->apply_varying_constraints(
                            buckets,
                            [&](const Scalar x, const Scalar y, const Scalar z) -> Scalar {
                                return varying_bc->eval_time_derivative(x, y, z);
                            },
                            bc.component,
                            v);
                    } else {
                        auto &buckets =
                            bulk_data.get_buckets(::stk::topology::NODE_RANK, *part & meta_data.locally_owned_part());
                        impl_->apply_varying_constraints(
                            buckets,
                            [&](const Scalar x, const Scalar y, const Scalar z) -> Scalar {
                                return varying_bc->eval_time_derivative(x, y, z);
                            },
                            bc.component,
                            v);
                    }
                }
            }
#else
            return;
#endif
        }

        void FunctionSpace::set_overwrite_vector(const Vector &v) {
            bool has_ow_bc = false;
            for (auto &bc_ptr : impl_->dirichlet_boundary) {
                if (bc_ptr->is_uniform()) continue;

                if (std::dynamic_pointer_cast<DirichletBoundary::OverwriteCondition>(bc_ptr)) {
                    has_ow_bc = true;
                    break;
                }
            }

            if (!has_ow_bc) return;

            // v.comm().root_print("SUCCESS set_overwrite_vector!");

            auto copy = std::make_shared<Vector>(v);

            for (auto &bc_ptr : impl_->dirichlet_boundary) {
                if (bc_ptr->is_uniform()) continue;

                auto ow_bc = std::dynamic_pointer_cast<DirichletBoundary::OverwriteCondition>(bc_ptr);

                if (ow_bc) {
                    ow_bc->vector = copy;
                }
            }
        }

        void FunctionSpace::copy_at_constrained_nodes(const Vector &source_vector, Vector &dest_vector) const {
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            auto &&local_to_global = dof_map().local_to_global();

            const int nv = n_var();

            if (local_to_global.empty()) {
                assert(comm().size() == 1);

                auto source_view = local_view_device(source_vector);
                auto dest_view = local_view_device(dest_vector);

                for (auto &bc_ptr : impl_->dirichlet_boundary) {
                    auto &bc = *bc_ptr;

                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
                        auto &buckets =
                            bulk_data.get_buckets(::stk::topology::NODE_RANK, *part & meta_data.locally_owned_part());

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                // auto idx = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node));
                                auto idx = utopia::stk::convert_entity_to_index(node);
                                auto index = idx * nv + bc.component;
                                dest_view.set(index, source_view.get(index));
                            }
                        }
                    }
                }
            } else {
                // auto r = range(v);
                // auto v_view = view_device(v);

                Write<Vector> w(dest_vector, utopia::GLOBAL_INSERT);
                Read<Vector> r(source_vector);

                for (auto &bc_ptr : impl_->dirichlet_boundary) {
                    auto &bc = *bc_ptr;
                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
                        auto &buckets =
                            bulk_data.get_buckets(::stk::topology::NODE_RANK, *part & meta_data.locally_owned_part());

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                auto local_idx = utopia::stk::convert_entity_to_index(node);
                                assert(local_idx < local_to_global.size());

                                auto index = local_to_global(local_idx, bc.component);
                                dest_vector.c_set(index, source_vector.get(index));
                            }
                        }
                    }
                }
            }
        }

        void FunctionSpace::apply_zero_constraints(Vector &v) const {
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            auto &&local_to_global = dof_map().local_to_global();

            const int nv = n_var();

            if (local_to_global.empty()) {
                assert(comm().size() == 1);

                auto v_view = local_view_device(v);

                for (auto &bc_ptr : impl_->dirichlet_boundary) {
                    auto &bc = *bc_ptr;

                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
                        auto &buckets =
                            bulk_data.get_buckets(::stk::topology::NODE_RANK, *part & meta_data.locally_owned_part());

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                // auto idx = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node));
                                auto idx = utopia::stk::convert_entity_to_index(node);
                                v_view.set(idx * nv + bc.component, 0.0);
                            }
                        }
                    }
                }
            } else {
                // auto r = range(v);
                // auto v_view = view_device(v);

                Write<Vector> w(v, utopia::GLOBAL_INSERT);

                for (auto &bc_ptr : impl_->dirichlet_boundary) {
                    auto &bc = *bc_ptr;
                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
                        auto &buckets =
                            bulk_data.get_buckets(::stk::topology::NODE_RANK, *part & meta_data.locally_owned_part());

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                auto local_idx = utopia::stk::convert_entity_to_index(node);
                                assert(local_idx < local_to_global.size());

                                v.c_set(local_to_global(local_idx, bc.component), 0.0);
                            }
                        }
                    }
                }
            }
        }

        void FunctionSpace::apply_constraints(Matrix &m, Vector &v) const {
            apply_constraints(m);
            apply_constraints(v);
        }

        bool FunctionSpace::empty() const { return !impl_->mesh || impl_->mesh->empty(); }

        void FunctionSpace::create_boundary_node_list(IndexArray &node_list) const {
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            const SizeType nl_nodes = mesh().n_local_nodes();
            node_list.clear();
            node_list.reserve(nl_nodes);

            const int nv = n_var();

            auto &&local_to_global = dof_map().local_to_global();

            if (local_to_global.empty()) {
                for (auto &bc_ptr : impl_->dirichlet_boundary) {
                    auto &bc = *bc_ptr;
                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
                        auto &buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, *part);

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                // auto idx = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node));
                                SizeType idx = utopia::stk::convert_entity_to_index(node);
                                assert(idx < nl_nodes);
                                node_list.push_back(idx);
                            }
                        }
                    }
                }
            } else {
                // auto rr = row_range(m);

                for (auto &bc_ptr : impl_->dirichlet_boundary) {
                    auto &bc = *bc_ptr;

                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
                        auto &buckets =
                            bulk_data.get_buckets(::stk::topology::NODE_RANK, *part & meta_data.locally_owned_part());

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                auto local_idx = utopia::stk::convert_entity_to_index(node);
                                assert(local_idx < local_to_global.size());
                                node_list.push_back(local_to_global.global_node_id(local_idx));
                            }
                        }
                    }
                }
            }
        }

        void FunctionSpace::add_dirichlet_boundary_condition(const std::string &name,
                                                             const Scalar &value,
                                                             const int component) {
            assert(component < n_var());

            DirichletBoundary::UniformCondition dirichlet_boundary{name, value, component};
            impl_->dirichlet_boundary.add(dirichlet_boundary);
        }

        // void FunctionSpace::displacement_field_from_transform(const std::vector<Scalar> &scale_factors,
        //                                                       Field<FunctionSpace> &displacement) {
        //     if (displacement.empty()) {
        //         this->create_field(displacement);
        //     }

        //     auto &&space = *this;

        //     int n_var = space.n_var();

        //     assert(n_var == mesh().spatial_dimension());

        //     auto d_view = view_device(displacement.data());

        //     Range r = range(displacement.data());
        //     SizeType r_begin = r.begin() / n_var;
        //     SizeType r_end = r.end() / n_var;

        //     const int n_factors = scale_factors.size();

        //     int dim = std::min(mesh().spatial_dimension(), n_var);

        //     auto fun = [=](const SizeType idx, const Scalar *point) {
        //         if (idx < r_begin || idx >= r_end) return;

        //         Scalar p3[3] = {0.0, 0.0, 0.0};
        //         Scalar transformed_p3[3] = {0.0, 0.0, 0.0};

        //         for (int d = 0; d < dim; ++d) {
        //             p3[d] = point[d];
        //         }

        //         for (int i = 0; i < n_factors; ++i) {
        //             transformed_p3[i] = p3[i] * scale_factors[i];
        //         }

        //         for (int d = 0; d < dim; ++d) {
        //             d_view.set(idx * n_var + d, transformed_p3[d]);
        //         }
        //     };

        //     space.node_eval(fun);
        // }

        void FunctionSpace::displace(const Vector &displacement) {
            Vector local_displacement;
            global_to_local(displacement, local_displacement);

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            ::stk::mesh::Selector s_universal = meta_data.universal_part();
            const auto &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, s_universal);
            auto *coords = meta_data.coordinate_field();

            // FIXME use parallel for etc...
            auto u_view = const_local_view_device(local_displacement);
            const int dim = mesh().spatial_dimension();

            assert(n_var() == dim);

            for (const auto &ib : node_buckets) {
                const auto &b = *ib;
                const SizeType length = b.size();

                for (SizeType k = 0; k < length; ++k) {
                    auto node = b[k];
                    auto idx = utopia::stk::convert_entity_to_index(node);

                    if (bulk_data.in_receive_ghost(node)) {
                        idx = dof_map().shift_aura_idx(idx);
                    }

                    Scalar *points = (Scalar *)::stk::mesh::field_data(*coords, node);

                    for (int d = 0; d < dim; ++d) {
                        points[d] += u_view.get(idx * dim + d);
                    }
                }
            }
        }

        void FunctionSpace::node_eval(std::function<void(const SizeType idx, const Scalar *)> fun) {
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            // ::stk::mesh::Selector s_universal = meta_data.universal_part();
            ::stk::mesh::Selector selector = meta_data.locally_owned_part();
            const auto &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, selector);
            impl_->node_eval(node_buckets, fun);
        }

        void FunctionSpace::node_eval(const std::string &part_name,
                                      std::function<void(const SizeType idx, const Scalar *)> fun) {
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            // ::stk::mesh::Selector s_universal = meta_data.universal_part();
            ::stk::mesh::Selector selector = meta_data.locally_owned_part();

            auto *part = meta_data.get_part(part_name);
            if (!part) {
                // FIXME Handle user space names better
                auto converted_part_name = SideSet::Cube::convert(part_name);
                part = meta_data.get_part(converted_part_name);
            }

            if (part) {
                auto &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, *part & selector);
                impl_->node_eval(node_buckets, fun);
            } else {
                Utopia::Abort("Part not defined!");
            }
        }

        const DofMap &FunctionSpace::dof_map() const {
            assert(impl_->dof_map);
            return *impl_->dof_map;
        }

        DofMap &FunctionSpace::dof_map() {
            assert(impl_->dof_map);
            return *impl_->dof_map;
        }

        void FunctionSpace::global_to_local(const Vector &global, Vector &local) const {
            dof_map().global_to_local(global, local);
        }

        void FunctionSpace::local_to_global(const Vector &local, Vector &global, AssemblyMode mode) const {
            if (global.empty()) {
                create_vector(global);
                global.set(0);
            }

            dof_map().local_to_global(local, global, mode);
        }

        void FunctionSpace::nodal_field_to_local_vector(Vector &v) { impl_->nodal_field_to_local_vector(v); }
        void FunctionSpace::local_vector_to_nodal_field(const Vector &v) { impl_->local_vector_to_nodal_field(v); }

        void FunctionSpace::nodal_field_to_global_vector(Vector &v) {
            if (comm().size() > 1) {
                Vector local_v;
                create_local_vector(local_v);
                nodal_field_to_local_vector(local_v);
                local_to_global(local_v, v, OVERWRITE_MODE);
            } else {
                nodal_field_to_local_vector(v);
            }
        }
        void FunctionSpace::global_vector_to_nodal_field(const Vector &v) {
            Vector local_v;
            global_to_local(v, local_v);
            local_vector_to_nodal_field(local_v);
        }

        const std::string &FunctionSpace::name() const { return impl_->name; }

        FunctionSpace::SizeType FunctionSpace::n_variables() const {
            return static_cast<SizeType>(impl_->variables.size());
        }

        const std::string &FunctionSpace::variable_name(const SizeType var_num) const {
            assert(var_num < n_variables());
            return impl_->variables[var_num].name;
        }

        int FunctionSpace::variable_size(const SizeType var_num) const {
            assert(var_num < n_variables());
            return impl_->variables[var_num].n_components;
        }

        void FunctionSpace::register_output_variables(MeshIO &io) {
            for (auto &v : impl_->variables) {
                io.register_output_field(v.name);
            }
        }

        template <typename FieldType>
        void FunctionSpace::declare_new_nodal_field(const std::string &name, const int n_comp) {
            auto &meta_data = mesh().meta_data();
            auto &&part = meta_data.universal_part();

            FEVar v;
            v.name = name;
            v.n_components = n_comp;
            impl_->variables.push_back(v);

            dof_map().set_n_var(n_var() + n_comp);

            meta_data.enable_late_fields();

            if (n_comp == 1) {
                auto &field = meta_data.declare_field<::stk::mesh::Field<Scalar>>(::stk::topology::NODE_RANK, name, 1);
                ::stk::mesh::put_field_on_mesh(field, part, 1, nullptr);
            } else {
                auto &field = meta_data.declare_field<Impl::VectorField_t>(::stk::topology::NODE_RANK, name, 1);
                ::stk::mesh::put_field_on_mesh(field, part, n_comp, nullptr);
            }
        }

        void FunctionSpace::backend_set_elemental_field(const Field<FunctionSpace> &field) {
            auto &meta_data = mesh().meta_data();
            auto &&part = meta_data.universal_part();

            ::stk::mesh::FieldBase *stk_field = ::stk::mesh::get_field_by_name(field.name(), meta_data);

            int tensor_size = field.tensor_size();
            if (!stk_field) {
                meta_data.enable_late_fields();

                if (tensor_size == 1) {
                    stk_field = &meta_data.declare_field<::stk::mesh::Field<Scalar>>(
                        ::stk::topology::ELEMENT_RANK, field.name(), 1);
                    ::stk::mesh::put_field_on_mesh(*stk_field, part, 1, nullptr);
                } else if (tensor_size == mesh().spatial_dimension()) {
                    assert(tensor_size == mesh().spatial_dimension());
                    stk_field =
                        &meta_data.declare_field<Impl::VectorField_t>(::stk::topology::ELEMENT_RANK, field.name(), 1);
                    ::stk::mesh::put_field_on_mesh(*stk_field, part, tensor_size, nullptr);
                } else if (tensor_size == mesh().spatial_dimension() * mesh().spatial_dimension()) {
                    // if (mesh().spatial_dimension() == 3) {
                    //     stk_field = &meta_data.declare_field<Impl::MatrixField3x3_t>(
                    //         ::stk::topology::ELEMENT_RANK, field.name(), 1);
                    //     ::stk::mesh::put_field_on_mesh(*stk_field, part, 3, 3, nullptr);
                    // } else {
                    stk_field =
                        &meta_data.declare_field<Impl::MatrixField_t>(::stk::topology::ELEMENT_RANK, field.name(), 1);
                    ::stk::mesh::put_field_on_mesh(
                        *stk_field, part, mesh().spatial_dimension(), mesh().spatial_dimension(), nullptr);
                    // }
                } else {
                    stk_field = &meta_data.declare_field<::stk::mesh::Field<Scalar>>(
                        ::stk::topology::ELEMENT_RANK, field.name(), 1);
                    ::stk::mesh::put_field_on_mesh(*stk_field, part, tensor_size, nullptr);
                }
            }

            auto &data = field.data();

            auto &bulk_data = mesh().bulk_data();
            auto &elem_buckets = utopia::stk::local_elements(bulk_data);

            auto v_view = local_view_device(data);

            using Bucket_t = ::stk::mesh::Bucket;
            using Entity_t = ::stk::mesh::Entity;

            int elem_idx = 0;
            for (const auto &ib : elem_buckets) {
                const Bucket_t &b = *ib;
                const Bucket_t::size_type length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k, ++elem_idx) {
                    Entity_t elem = b[k];

                    Scalar *points = (Scalar *)::stk::mesh::field_data(*stk_field, elem);
                    int n_comp = ::stk::mesh::field_scalars_per_entity(*stk_field, elem);
                    assert(field.tensor_size() == n_comp);

                    for (int d = 0; d < n_comp; ++d) {
                        points[d] = v_view.get(elem_idx * n_comp + d);
                    }
                }
            }
        }

        void FunctionSpace::backend_set_nodal_field(const Field<FunctionSpace> &field) {
            auto &meta_data = mesh().meta_data();
            auto &&part = meta_data.universal_part();

            ::stk::mesh::FieldBase *stk_field = ::stk::mesh::get_field_by_name(field.name(), meta_data);

            int tensor_size = field.tensor_size();
            if (!stk_field) {
                meta_data.enable_late_fields();

                if (tensor_size == 1) {
                    stk_field = &meta_data.declare_field<::stk::mesh::Field<Scalar>>(
                        ::stk::topology::NODE_RANK, field.name(), 1);
                    ::stk::mesh::put_field_on_mesh(*stk_field, part, 1, nullptr);
                } else {
                    assert(tensor_size == mesh().spatial_dimension());

                    stk_field =
                        &meta_data.declare_field<Impl::VectorField_t>(::stk::topology::NODE_RANK, field.name(), 1);
                    ::stk::mesh::put_field_on_mesh(*stk_field, part, tensor_size, nullptr);
                }
            }

            auto &data = field.data();

            if (data.comm().size() == 1) {
                impl_->copy_local_vector_to_nodal_field(data, *stk_field);
            } else {
                Vector local;
                global_to_local(data, local);
                impl_->copy_local_vector_to_nodal_field(local, *stk_field);
            }
        }

        void FunctionSpace::create_vector(const std::string &field_name, Vector &v) const {
            int nv = -1;
            for (auto &v : impl_->variables) {
                if (v.name == field_name) {
                    nv = v.n_components;
                }
            }

            assert(nv != -1);
            if (nv == -1) {
                Utopia::Abort("Field with name: " + field_name + "does not exist!");
            }

            v.zeros(layout(comm(), mesh().n_local_nodes() * nv, mesh().n_nodes() * nv));
            v.set_block_size(n_var());
        }

        void FunctionSpace::global_vector_to_nodal_field(const std::string &field_name, const Vector &v) {
            Vector local;
            global_to_local(v, local);
            impl_->local_vector_to_nodal_field(field_name, local);
        }

        void FunctionSpace::set_print_map(const bool val) { impl_->print_map = val; }

        template void FunctionSpace::declare_new_nodal_field<Traits<FunctionSpace>::Scalar>(const std::string &,
                                                                                            const int);
        // template void FunctionSpace::declare_new_nodal_field<Traits<FunctionSpace>::SizeType>(const std::string
        // &, const int);
        template void FunctionSpace::declare_new_nodal_field<int>(const std::string &, const int);

    }  // namespace stk
}  // namespace utopia
