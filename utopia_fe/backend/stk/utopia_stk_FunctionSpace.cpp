#include "utopia_stk_FunctionSpace.hpp"

#include "utopia_Options.hpp"

#include "utopia_stk_Commons.hpp"
#include "utopia_stk_DofMap.hpp"

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

        class DirichletBoundary : public Configurable, public Describable {
        public:
            class Condition : public Configurable, public Describable {
            public:
                std::string name;
                double value;
                int component{0};
                int side{-1};

                Condition() = default;
                Condition(std::string name, double value, const int component)
                    : name(std::move(name)), value(value), component(component) {}

                void read(Input &in) override {
                    in.get("name", name);
                    in.get("value", value);
                    in.get("var", component);
                    in.get("side", side);

                    if (name.empty() && side != -1) {
                        name = "surface_" + std::to_string(side);
                    }
                }

                void describe(std::ostream &os) const override {
                    os << "name:\t" << name << '\n';
                    os << "value:\t" << value << '\n';
                    os << "var:\t" << component << '\n';
                    os << "side:\t" << side << '\n';
                }
            };

            void read(Input &in) override {
                in.get_all([this](Input &in) {
                    conditions.emplace_back();
                    auto &c = conditions.back();
                    c.read(in);
                });
            }

            void describe(std::ostream &os) const override {
                os << "DirichletBoundary:\n";
                for (auto &c : conditions) {
                    c.describe(os);
                }
            }

            std::vector<Condition> conditions;
        };

        class FunctionSpace::Impl {
        public:
            using Bucket_t = ::stk::mesh::Bucket;
            using BucketVector_t = ::stk::mesh::BucketVector;
            using Entity_t = ::stk::mesh::Entity;
            using Selector_t = ::stk::mesh::Selector;
            using IOBroker_t = ::stk::io::StkMeshIoBroker;

            std::shared_ptr<Mesh> mesh;
            DirichletBoundary dirichlet_boundary;
            std::shared_ptr<DofMap> dof_map;

            void vector_to_nodal_field_old(const std::string &name, const Vector &v) {
                auto &meta_data = mesh->meta_data();
                auto &bulk_data = mesh->bulk_data();

                // ::stk::mesh::Selector s_universal = meta_data.universal_part();
                // const BucketVector_t &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, s_universal);

                auto &node_buckets = utopia::stk::local_nodes(bulk_data);
                ::stk::mesh::FieldBase *field = ::stk::mesh::get_field_by_name(name, meta_data);
                assert(field);

                auto v_view = view_device(v);

                auto &&local_to_global = dof_map->local_to_global();
#ifndef NDEBUG
                const int n_var = dof_map->n_var();
#endif  // NDEBUG

                for (const auto &ib : node_buckets) {
                    const Bucket_t &b = *ib;
                    const Bucket_t::size_type length = b.size();

                    for (Bucket_t::size_type k = 0; k < length; ++k) {
                        Entity_t node = b[k];
                        auto idx = utopia::stk::convert_entity_to_index(node);

                        Scalar *points = (Scalar *)::stk::mesh::field_data(*field, node);
                        auto n_comp = ::stk::mesh::field_scalars_per_entity(*field, node);
                        assert(n_comp == n_var);

                        if (!local_to_global.empty()) {
                            for (int d = 0; d < n_comp; ++d) {
                                points[d] = v_view.get(idx * n_comp + d);
                            }
                        } else {
                            for (int d = 0; d < n_comp; ++d) {
                                points[d] = v_view.get(local_to_global(idx, d));
                            }
                        }
                    }
                }
            }

            void vector_to_nodal_field(const std::string &name, const Vector &v) {
                auto &meta_data = mesh->meta_data();
                auto &bulk_data = mesh->bulk_data();

                // ::stk::mesh::Selector s_universal = meta_data.universal_part();
                // const BucketVector_t &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, s_universal);

                auto &node_buckets = utopia::stk::universal_nodes(bulk_data);
                ::stk::mesh::FieldBase *field = ::stk::mesh::get_field_by_name(name, meta_data);
                assert(field);

                auto v_view = local_view_device(v);
#ifndef NDEBUG
                const int n_var = dof_map->n_var();
#endif  // NDEBUG

                for (const auto &ib : node_buckets) {
                    const Bucket_t &b = *ib;
                    const Bucket_t::size_type length = b.size();

                    for (Bucket_t::size_type k = 0; k < length; ++k) {
                        Entity_t node = b[k];
                        auto idx = utopia::stk::convert_entity_to_index(node);

                        Scalar *points = (Scalar *)::stk::mesh::field_data(*field, node);
                        auto n_comp = ::stk::mesh::field_scalars_per_entity(*field, node);
                        assert(n_comp == n_var);

                        for (int d = 0; d < n_comp; ++d) {
                            points[d] = v_view.get(idx * n_comp + d);
                        }
                    }
                }
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

        FunctionSpace::~FunctionSpace() = default;

        bool FunctionSpace::write(const Path &path, const Vector &x) {
            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            meta_data.enable_late_fields();

            // FIXME declare this somewhere else
            ::stk::mesh::Field<Scalar> &field =
                meta_data.declare_field<::stk::mesh::Field<Scalar>>(::stk::topology::NODE_RANK, "output", 1);

            // const std::vector<Scalar> init(this->n_var(), 1);
            // ::stk::mesh::put_field_on_mesh(field, meta_data.universal_part(), this->n_var(), init.data());
            ::stk::mesh::put_field_on_mesh(field, meta_data.universal_part(), this->n_var(), nullptr);

            if (comm().size() == 1) {
                impl_->vector_to_nodal_field("output", x);
            } else {
                Vector x_local;
                global_to_local(x, x_local);
                impl_->vector_to_nodal_field("output", x_local);
            }

            try {
                Impl::IOBroker_t io_broker(comm().raw_comm());
                io_broker.set_bulk_data(bulk_data);

                auto out_id = io_broker.create_output_mesh(path.to_string(), ::stk::io::WRITE_RESTART);
                io_broker.add_field(out_id, field, "output");

                // io_broker.write_output_mesh(out_id);
                io_broker.process_output_request(out_id, 0);
                return true;

            } catch (const std::exception &ex) {
                utopia::err() << "Mesh::write(\"" << path.to_string() << "\") error: " << ex.what() << '\n';
                assert(false);
                return false;
            }
        }

        void FunctionSpace::init(const std::shared_ptr<Mesh> &mesh) { impl_->mesh = mesh; }

        void FunctionSpace::read(Input &in) {
            auto mesh = std::make_shared<Mesh>();

            // mesh->read(in);
            in.get("mesh", *mesh);
            init(mesh);

            int n_var = 0;
            bool verbose = false;

            if (!Options()
                     .add_option("n_var", n_var, "Number of variables per node (instead of specifiying all of them).")
                     .add_option("boundary_conditions", impl_->dirichlet_boundary, "Boundary conditions.")
                     .add_option("verbose", verbose, "Verbose output.")
                     .parse(in)) {
                return;
            }

            int counted_vars = 0;
            in.get("variables", [&counted_vars](Input &in) { in.get_all([&](Input &in) { ++counted_vars; }); });

            assert(n_var == 0 || n_var == counted_vars);

            n_var = counted_vars;

            dof_map().set_n_var(n_var);

            if (verbose) {
                std::stringstream ss;
                describe(ss);

                comm().synched_print(ss.str());
            }
        }

        void FunctionSpace::read_with_state(Input &in, Vector &val) {
            read(in);
            assert(false && "IMPLEMENT ME");
        }

        void FunctionSpace::describe(std::ostream &os) const {
            // using Bucket_t = ::stk::mesh::Bucket;
            // using BucketVector_t = ::stk::mesh::BucketVector;
            // using Entity_t = ::stk::mesh::Entity;

            if (comm().rank() == 0) {
                os << "Parts:\n";
                for (auto ptr : mesh().meta_data().get_parts()) {
                    auto &p = *ptr;
                    if (p.id() != -1) {
                        os << p.name() << ' ' << p.id() << '\n';
                    }
                }

                impl_->dirichlet_boundary.describe(os);
            }

            os << '\n';
            os << "n_vars: " << n_var() << '\n';
            os << "n_local_dofs: " << n_local_dofs() << '\n';

            // auto &meta_data = mesh().meta_data();
            // auto &bulk_data = mesh().bulk_data();

            // const ::stk::mesh::Selector selector = meta_data.locally_owned_part();
            // const auto &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, selector);

            // for (auto *b_ptr : node_buckets) {
            //     const auto &b = *b_ptr;
            //     const auto length = b.size();

            //     for (Bucket_t::size_type k = 0; k < length; ++k) {
            //         const Entity_t elem = b[k];
            //         const auto id = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(elem));
            //         os << id << ' ';
            //     }

            //     os << '\n';
            // }
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

        void FunctionSpace::create_vector(Vector &v) const { v.zeros(layout(comm(), n_local_dofs(), n_dofs())); }

        void FunctionSpace::create_local_vector(Vector &v) const {
            // FIXME
            v.zeros(layout(comm(), n_local_dofs(), n_dofs()));
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

        void FunctionSpace::apply_constraints(Matrix &m) {
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            using IndexSet = Traits<Matrix>::IndexSet;

            IndexSet constrains;
            constrains.reserve(n_local_dofs());

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
                                auto idx = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node));
                                constrains.push_back(idx * nv + bc.component);
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

            set_zero_rows(m, constrains, 1.0);
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
                                auto idx = utopia::stk::convert_stk_index_to_index(bulk_data.identifier(node));
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

        void FunctionSpace::apply_constraints(Matrix &m, Vector &v) {
            apply_constraints(m);
            apply_constraints(v);
        }

        bool FunctionSpace::empty() const { return !impl_->mesh || impl_->mesh->empty(); }

        void FunctionSpace::add_dirichlet_boundary_condition(const std::string &name,
                                                             const Scalar &value,
                                                             const int component) {
            DirichletBoundary::Condition dirichlet_boundary{name, value, component};
            impl_->dirichlet_boundary.conditions.push_back(dirichlet_boundary);
        }

        void FunctionSpace::displace(const Vector &displacement) { assert(false && "IMPLEMENT ME"); }

        const DofMap &FunctionSpace::dof_map() const { return *impl_->dof_map; }
        DofMap &FunctionSpace::dof_map() { return *impl_->dof_map; }

        void FunctionSpace::global_to_local(const Vector &global, Vector &local) const {
            dof_map().global_to_local(global, local);
        }

    }  // namespace stk
}  // namespace utopia
