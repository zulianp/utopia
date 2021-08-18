#include "utopia_stk_Mesh.hpp"

#include "utopia_make_unique.hpp"

#include "utopia_stk_Commons.hpp"

#include "utopia_stk_MeshIO.hpp"

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <cstdio>

namespace utopia {
    namespace stk {

        class Mesh::Impl {
        public:
            using MetaData = ::stk::mesh::MetaData;
            using BulkData = ::stk::mesh::BulkData;
            using IOBroker = ::stk::io::StkMeshIoBroker;

            Comm comm;
            std::shared_ptr<MetaData> meta_data;
            std::shared_ptr<BulkData> bulk_data;

            SizeType n_elements{-1};
            SizeType n_nodes{-1};
            SizeType n_local_elements{-1};
            SizeType n_local_nodes{-1};
            bool verbose{false};
            bool is_generated_cube{false};

            void compute_mesh_stats() {
                std::vector<size_t> entity_counts;
                ::stk::mesh::comm_mesh_counts(*bulk_data, entity_counts);

                n_elements = entity_counts[::stk::topology::ELEMENT_RANK];
                n_nodes = entity_counts[::stk::topology::NODE_RANK];

                n_local_elements = utopia::stk::count_local_elements(*bulk_data);
                n_local_nodes = utopia::stk::count_local_nodes(*bulk_data);
            }
        };

        Mesh::~Mesh() = default;
        Mesh::Mesh(const Comm &comm) : impl_(utopia::make_unique<Impl>()) { impl_->comm = comm; }

        // 1736
        // https://github.com/NaluCFD/Nalu/blob/master/src/Realm.C

        // automatic_decomposition_type (DECOMPOSITION_METHOD)
        // Used only for parallel runs, this indicates how the a single mesh database must be decomposed amongst the MPI
        // processes during initialization. This option should not be used if the mesh has already been decomposed by an
        // external utility. Possible values are: Purpose Generate statistics for the flow field Extract integrated data
        // from the simulation Compare the solution error to a reference solution Extract data using probes Model
        // turbine blades/tower using actuator lines Momentum source term to drive ABL flows to a desired velocity
        // profile
        //         The name of the realm. The name provided here is used in the Time_Integrators section to determine
        //         the time-integration scheme used for this computational domain.
        // The name of the Exodus-II mesh file that defines the computational domain for this realm. Note that only the
        // base name (i.e., without the .NPROCS.IPROC suffix) is provided even for simulations using previously
        // decomposed mesh/restart files. Value Description rcb recursive coordinate bisection rib recursive inertial
        // bisection linear elements in order first n/p to proc 0, next to proc 1. cyclic elements handed out to id %
        // proc_count

        bool Mesh::read(const Path &path) {
            MeshIO io(*this);
            io.set_read_path(path);
            return io.load();
        }

        bool Mesh::write(const Path &path) {
            MeshIO io(*this);
            return io.write(path);
        }

        void Mesh::read(Input &in) {
            MeshIO io(*this);
            io.read(in);
            if (!io.load()) {
                assert(false);
            }

            Scalar rescale = 1.0;
            in.get("rescale", rescale);

            if (rescale != 1.0) {
                scale(rescale);
            }

            in.get("verbose", impl_->verbose);
        }

        void Mesh::unit_cube(const SizeType &nx, const SizeType &ny, const SizeType &nz) {
            MeshIO io(*this);

            char format[100];
            std::sprintf(format, "generated:%dx%dx%d|sideset:xX", int(nx), int(ny), int(nz));
            io.set_read_specification(format);
            if (!io.load()) {
                assert(false);
            }

            set_is_generated_cube(true);
        }

        void Mesh::box(const AABB &box,
                       const SizeType &nx,
                       const SizeType &ny,
                       const SizeType &nz,
                       const std::string &elem_type) {
            InputParameters params;

            const int dim = box.min.size();

            params.set("nx", nx);
            params.set("ny", ny);
            params.set("nz", nz);

            const std::string postfix[4] = {"x", "y", "z", "t"};

            for (int d = 0; d < dim; ++d) {
                params.set("min_" + postfix[d], box.min[d]);
                params.set("max_" + postfix[d], box.max[d]);
            }

            params.set("elem_type", elem_type);

            MeshIO io(*this);
            io.read(params);

            if (!io.load()) {
                assert(false);
                Utopia::Abort("Failed to generate box!");
            }
        }

        void Mesh::describe(std::ostream &os) const {
            if (comm().rank() == 0) {
                os << "Parts:\n";
                for (auto ptr : meta_data().get_parts()) {
                    auto &p = *ptr;
                    if (p.id() != -1) {
                        os << p.name() << ' ' << p.id() << '\n';
                    }
                }

                for (auto &field : impl_->meta_data->get_fields()) {
                    os << field->name() << ", " << field->entity_rank() << ", num states: " << field->number_of_states()
                       << '\n';
                }
            }

            std::stringstream ss;

            ss << "n_local_nodes:\t\t" << n_local_nodes() << "\n";
            ss << "n_aura_nodes:\t\t" << count_aura_nodes(bulk_data()) << "\n";
            ss << "n_shared_nodes:\t\t" << count_shared_nodes(bulk_data()) << "\n";
            ss << "n_universal_nodes:\t" << count_universal_nodes(bulk_data()) << "\n";
            ss << "\n";
            ss << "n_local_elements:\t" << n_local_elements() << "\n";
            ss << "n_aura_elements:\t" << count_aura_elements(bulk_data()) << "\n";
            ss << "n_shared_elements:\t" << count_shared_elements(bulk_data()) << "\n";
            ss << "n_universal_elements:\t" << count_universal_elements(bulk_data()) << "\n";

            comm().synched_print(ss.str());

            // if (impl_->verbose) {
            //     // impl_->bulk_data->dump_all_mesh_info(os);
            // }
        }

        const Mesh::Comm &Mesh::comm() const { return impl_->comm; }

        ::stk::mesh::BulkData &Mesh::bulk_data() const {
            assert(impl_->bulk_data);
            return *impl_->bulk_data;
        }

        ::stk::mesh::MetaData &Mesh::meta_data() const {
            assert(impl_->meta_data);
            return *impl_->meta_data;
        }

        void Mesh::wrap(const std::shared_ptr<::stk::mesh::MetaData> &meta_data,
                        const std::shared_ptr<::stk::mesh::BulkData> &bulk_data) {
            impl_->meta_data = meta_data;
            impl_->bulk_data = bulk_data;

            if (!empty()) {
                impl_->compute_mesh_stats();
            }
        }

        bool Mesh::empty() const { return !static_cast<bool>(impl_->bulk_data); }

        int Mesh::spatial_dimension() const {
            assert(impl_->meta_data);
            return impl_->meta_data->spatial_dimension();
        }

        Mesh::SizeType Mesh::n_elements() const { return impl_->n_elements; }

        Mesh::SizeType Mesh::n_nodes() const { return impl_->n_nodes; }

        Mesh::SizeType Mesh::n_local_elements() const { return impl_->n_local_elements; }

        Mesh::SizeType Mesh::n_local_nodes() const { return impl_->n_local_nodes; }

        void Mesh::displace(const Vector &) {
            Utopia::Abort("IMPLEMENT ME");
            assert(false);
        }

        void Mesh::init() { impl_->compute_mesh_stats(); }

        void Mesh::bounding_box(AABB &output) const {
            ::stk::mesh::Selector s_universal = meta_data().universal_part();
            const auto &node_buckets = bulk_data().get_buckets(::stk::topology::NODE_RANK, s_universal);
            auto *coords = meta_data().coordinate_field();

            const int dim = spatial_dimension();

            std::vector<Scalar> minmax(2 * dim);

            for (const auto &ib : node_buckets) {
                const auto &b = *ib;
                const SizeType length = b.size();

                for (SizeType k = 0; k < length; ++k) {
                    auto node = b[k];
                    Scalar *points = (Scalar *)::stk::mesh::field_data(*coords, node);

                    for (int d = 0; d < dim; ++d) {
                        minmax[d] = std::min(minmax[d], points[d]);
                        minmax[dim + d] = std::min(minmax[dim + d], -points[d]);
                    }
                }
            }

            this->comm().min(2 * dim, &minmax[0]);

            output.min.resize(dim);
            output.max.resize(dim);

            for (int d = 0; d < dim; ++d) {
                output.min[d] = minmax[d];
                output.max[d] = -minmax[dim + d];
            }
        }

        void Mesh::scale(const Scalar &scale_factor) {
            ::stk::mesh::Selector s_universal = meta_data().universal_part();
            const auto &node_buckets = bulk_data().get_buckets(::stk::topology::NODE_RANK, s_universal);
            auto *coords = meta_data().coordinate_field();

            const int dim = spatial_dimension();

            for (const auto &ib : node_buckets) {
                const auto &b = *ib;
                const SizeType length = b.size();

                for (SizeType k = 0; k < length; ++k) {
                    auto node = b[k];
                    Scalar *points = (Scalar *)::stk::mesh::field_data(*coords, node);

                    for (int d = 0; d < dim; ++d) {
                        points[d] *= scale_factor;
                    }
                }
            }
        }

        bool Mesh::has_aura() const { return bulk_data().is_automatic_aura_on(); }

        void Mesh::create_edges() {
            auto &edges_part = meta_data().declare_part(universal_edge_set_name(), ::stk::topology::EDGE_RANK);
            ::stk::mesh::create_edges(bulk_data(), meta_data().universal_part(), &edges_part);

            std::stringstream ss;
            ss << "n_universal_edges:\t"
               << utopia::stk::count_entities(bulk_data(), ::stk::topology::EDGE_RANK, edges_part) << "\n";

            comm().synched_print(ss.str());
        }

        void Mesh::set_is_generated_cube(const bool val) { impl_->is_generated_cube = val; }
    }  // namespace stk
}  // namespace utopia
