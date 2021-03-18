#include "utopia_stk_FunctionSpace.hpp"

#include "utopia_Options.hpp"

#include "utopia_stk_Commons.hpp"

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
            class Condition : public Configurable {
            public:
                std::string name;
                double value;
                int set_id{-1};

                Condition() = default;
                Condition(std::string name, double value) : name(std::move(name)), value(value) {}

                void read(Input &in) override {
                    in.get("name", name);
                    in.get("value", value);
                    in.get("set_id", set_id);
                }
            };

            void read(Input &in) override {
                in.get_all([this](Input &in) {
                    conditions.emplace_back();
                    auto &c = conditions.back();
                    c.read(in);
                });
            }

            void describe(std::ostream &os) const override { os << "DirichletBoundary"; }

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
            int n_var{1};
            DirichletBoundary dirichlet_boundary;

            void vector_to_nodal_field(const std::string &name, const Vector &v) {
                auto &meta_data = mesh->meta_data();
                auto &bulk_data = mesh->bulk_data();

                ::stk::mesh::Selector s_universal = meta_data.universal_part();
                const BucketVector_t &node_buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, s_universal);

                ::stk::mesh::FieldBase *field = ::stk::mesh::get_field_by_name(name, meta_data);
                assert(field);

                auto v_view = const_local_view_device(v);
                for (const auto &ib : node_buckets) {
                    const Bucket_t &b = *ib;
                    const Bucket_t::size_type length = b.size();

                    for (Bucket_t::size_type k = 0; k < length; ++k) {
                        Entity_t node = b[k];
                        auto idx = utopia::stk::convert_entity_to_index(node);

                        Scalar *points = (Scalar *)::stk::mesh::field_data(*field, node);
                        auto n_comp = ::stk::mesh::field_scalars_per_entity(*field, node);

                        for (int d = 0; d < n_comp; ++d) {
                            points[d] = v_view.get(idx * n_comp + d);
                        }
                    }
                }
            }
        };

        FunctionSpace::FunctionSpace(const Comm &comm) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = std::make_shared<Mesh>(comm);
        }

        FunctionSpace::FunctionSpace(const std::shared_ptr<Mesh> &mesh) : impl_(utopia::make_unique<Impl>()) {
            impl_->mesh = mesh;
        }

        FunctionSpace::~FunctionSpace() = default;

        bool FunctionSpace::write(const Path &path, const Vector &x) {
            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            meta_data.enable_late_fields();

            // FIXME declare this somewhere else
            ::stk::mesh::Field<Scalar> &field =
                meta_data.declare_field<::stk::mesh::Field<Scalar>>(::stk::topology::NODE_RANK, "output", impl_->n_var);

            const std::vector<Scalar> init(impl_->n_var, 1);
            ::stk::mesh::put_field_on_mesh(field, meta_data.universal_part(), impl_->n_var, init.data());
            impl_->vector_to_nodal_field("output", x);

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
            mesh->read(in);
            init(mesh);

            if (!Options()
                     .add_option("n_var", impl_->n_var, "Number of variables per node.")
                     .add_option("boundary_conditions", impl_->dirichlet_boundary, "Boundary conditions.")
                     .parse(in)) {
                return;
            }
        }

        void FunctionSpace::describe(std::ostream &os) const {
            // impl_->mesh->describe(os);
            os << "Parts:\n";
            for (auto ptr : mesh().meta_data().get_parts()) {
                auto &p = *ptr;
                os << p.name() << ' ' << p.id() << '\n';
            }
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

        FunctionSpace::SizeType FunctionSpace::n_dofs() const { return impl_->mesh->n_nodes() * impl_->n_var; }
        FunctionSpace::SizeType FunctionSpace::n_local_dofs() const {
            return impl_->mesh->n_local_nodes() * impl_->n_var;
        }

        void FunctionSpace::create_vector(Vector &v) const { v.zeros(layout(comm(), n_local_dofs(), n_dofs())); }

        void FunctionSpace::create_local_vector(Vector &v) const {
            // FIXME
            v.zeros(layout(comm(), n_local_dofs(), n_dofs()));
        }

        void FunctionSpace::create_matrix(Matrix &m) const {
            using Bucket_t = ::stk::mesh::Bucket;
            using BucketVector_t = ::stk::mesh::BucketVector;
            using Entity_t = ::stk::mesh::Entity;

            auto vl = layout(comm(), n_local_dofs(), n_dofs());
            auto ml = square_matrix_layout(vl);

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            const ::stk::mesh::Selector selector = meta_data.locally_owned_part();
            const BucketVector_t &node_buckets = bulk_data.get_buckets(::stk::topology::ELEMENT_RANK, selector);

            // Slow but for the moment it will do...
            const SizeType nln = mesh().n_local_nodes();
            std::vector<std::unordered_set<SizeType>> node2node(nln);

            for (auto *b_ptr : node_buckets) {
                const auto &b = *b_ptr;
                const auto length = b.size();

                for (Bucket_t::size_type k = 0; k < length; ++k) {
                    const Entity_t elem = b[k];
                    const auto node_ids = bulk_data.begin_nodes(elem);
                    const Size_t n_nodes = bulk_data.num_nodes(elem);

                    for (Size_t i = 0; i < n_nodes; ++i) {
                        const auto node_i = utopia::stk::convert_entity_to_index(node_ids[i]);

                        for (Size_t j = 0; j < n_nodes; ++j) {
                            const auto node_j = utopia::stk::convert_entity_to_index(node_ids[j]);
                            node2node[node_i].insert(node_j);
                        }
                    }
                }
            }

            IndexSet d_nnz(nln, 0), o_nnz(nln, 0);
            for (SizeType i = 0; i < nln; ++i) {
                d_nnz[i] = node2node[i].size();
            }

            if (impl_->n_var == 1) {
                m.sparse(ml, d_nnz, o_nnz);
            } else {
                m.block_sparse(ml, d_nnz, o_nnz, impl_->n_var);
            }
        }

        void FunctionSpace::apply_constraints(Matrix &m) {
            using Bucket_t = ::stk::mesh::Bucket;

            auto &meta_data = mesh().meta_data();
            auto &bulk_data = mesh().bulk_data();

            using IndexSet = Traits<Matrix>::IndexSet;

            IndexSet constrains;
            constrains.reserve(n_local_dofs());

            for (auto &bc : impl_->dirichlet_boundary.conditions) {
                auto *part = meta_data.get_part(bc.name);
                if (part) {
                    auto &buckets = bulk_data.get_buckets(::stk::topology::NODE_RANK, *part);

                    for (auto *b_ptr : buckets) {
                        auto &b = *b_ptr;
                        const Bucket_t::size_type length = b.size();

                        for (Bucket_t::size_type k = 0; k < length; ++k) {
                            auto node = b[k];
                            auto idx = utopia::stk::convert_entity_to_index(node);
                            constrains.push_back(idx);
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
                            auto idx = utopia::stk::convert_entity_to_index(node);
                            v_view.set(idx, bc.value);
                        }
                    }
                }
            }
        }

        void FunctionSpace::apply_constraints(Matrix &m, Vector &v) {
            apply_constraints(m);
            apply_constraints(v);
        }

        void FunctionSpace::add_dirichlet_boundary_condition(const std::string &name, const Scalar &value) {
            DirichletBoundary::Condition dirichlet_boundary{name, value};
            impl_->dirichlet_boundary.conditions.push_back(dirichlet_boundary);
        }

    }  // namespace stk
}  // namespace utopia
