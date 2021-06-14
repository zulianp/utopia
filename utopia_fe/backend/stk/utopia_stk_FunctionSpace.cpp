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

            std::string name{"main"};
            std::shared_ptr<Mesh> mesh;
            DirichletBoundary dirichlet_boundary;
            std::vector<FEVar> variables;
            std::shared_ptr<DofMap> dof_map;
            bool verbose{false};

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
                // auto &&part = meta_data.universal_part();

                meta_data.enable_late_fields();

                for (auto &v : variables) {
                    // if (v.n_components == 1) {
                    //     auto &field =
                    //         meta_data.declare_field<::stk::mesh::Field<Scalar>>(::stk::topology::NODE_RANK, v.name,
                    //         1);
                    //     ::stk::mesh::put_field_on_mesh(field, part, 1, nullptr);
                    // } else {
                    //     auto &field =
                    //         meta_data.declare_field<Impl::VectorField_t>(::stk::topology::NODE_RANK, v.name, 1);
                    //     ::stk::mesh::put_field_on_mesh(field, part, v.n_components, nullptr);
                    // }

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
        };

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

        void FunctionSpace::initialize() {
            impl_->register_variables();

            impl_->dof_map->init(this->mesh().bulk_data());

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

        void FunctionSpace::read_meta(Input &in) { impl_->read_meta(in); }

        void FunctionSpace::register_variables() { impl_->register_variables(); }

        void FunctionSpace::read(Input &in) {
            in.get("mesh", *impl_->mesh);
            read_meta(in);

            impl_->register_variables();

            impl_->dof_map->init(this->mesh().bulk_data());

            if (impl_->verbose) {
                std::stringstream ss;
                describe(ss);
                comm().synched_print(ss.str());
            }
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
            impl_->dof_map->init(this->mesh().bulk_data());

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
                impl_->dof_map->init(this->mesh().bulk_data());
            }

            auto vl = layout(comm(), n_local_dofs(), n_dofs());
            auto ml = square_matrix_layout(vl);

            if (this->n_var() == 1) {
                m.sparse(ml, impl_->dof_map->d_nnz(), impl_->dof_map->o_nnz());
            } else {
                m.block_sparse(ml, impl_->dof_map->d_nnz(), impl_->dof_map->o_nnz(), this->n_var());
            }
        }

        void FunctionSpace::apply_constraints(Matrix &m, const Scalar diag_value) {
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
                for (auto &bc : impl_->dirichlet_boundary.conditions) {
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

                for (auto &bc : impl_->dirichlet_boundary.conditions) {
                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
                        auto &buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, *part);

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

        void FunctionSpace::apply_constraints(Vector &v) {
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            auto &&local_to_global = dof_map().local_to_global();

            const int nv = n_var();

            if (local_to_global.empty()) {
                assert(comm().size() == 1);

                auto v_view = local_view_device(v);

                for (auto &bc : impl_->dirichlet_boundary.conditions) {
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
                                v_view.set(idx * nv + bc.component, bc.value);
                            }
                        }
                    }
                }
            } else {
                // auto r = range(v);
                // auto v_view = view_device(v);

                Write<Vector> w(v, utopia::GLOBAL_INSERT);

                for (auto &bc : impl_->dirichlet_boundary.conditions) {
                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
                        auto &buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, *part);

                        for (auto *b_ptr : buckets) {
                            auto &b = *b_ptr;
                            const Bucket_t::size_type length = b.size();

                            for (Bucket_t::size_type k = 0; k < length; ++k) {
                                auto node = b[k];
                                auto local_idx = utopia::stk::convert_entity_to_index(node);
                                assert(local_idx < local_to_global.size());

                                v.c_set(local_to_global(local_idx, bc.component), bc.value);
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

                for (auto &bc : impl_->dirichlet_boundary.conditions) {
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
                                v_view.set(idx * nv + bc.component, 0.0);
                            }
                        }
                    }
                }
            } else {
                // auto r = range(v);
                // auto v_view = view_device(v);

                Write<Vector> w(v, utopia::GLOBAL_INSERT);

                for (auto &bc : impl_->dirichlet_boundary.conditions) {
                    auto *part = meta_data.get_part(bc.name);
                    if (part) {
                        auto &buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, *part);

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

        void FunctionSpace::apply_constraints(Matrix &m, Vector &v) {
            apply_constraints(m);
            apply_constraints(v);
        }

        bool FunctionSpace::empty() const { return !impl_->mesh || impl_->mesh->empty(); }

        void FunctionSpace::add_dirichlet_boundary_condition(const std::string &name,
                                                             const Scalar &value,
                                                             const int component) {
            assert(component < n_var());

            DirichletBoundary::Condition dirichlet_boundary{name, value, component};
            impl_->dirichlet_boundary.conditions.push_back(dirichlet_boundary);
        }

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

                    Scalar *points = (Scalar *)::stk::mesh::field_data(*coords, node);

                    for (int d = 0; d < dim; ++d) {
                        points[d] += u_view.get(idx * dim + d);
                    }
                }
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

        template void FunctionSpace::declare_new_nodal_field<Traits<FunctionSpace>::Scalar>(const std::string &,
                                                                                            const int);
        // template void FunctionSpace::declare_new_nodal_field<Traits<FunctionSpace>::SizeType>(const std::string &,
        // const int);
        template void FunctionSpace::declare_new_nodal_field<int>(const std::string &, const int);

    }  // namespace stk
}  // namespace utopia
