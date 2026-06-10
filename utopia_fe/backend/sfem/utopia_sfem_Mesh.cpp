#include "utopia_sfem_Mesh.hpp"

#include "utopia_IOStream.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_make_unique.hpp"

#include "smesh_extractions.hpp"
#include "smesh_mesh.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <utopia_ArrayView.hpp>
#include <vector>

namespace utopia {
    namespace sfem {
        using idx_t = ::smesh::idx_t;
        using geom_t = ::smesh::geom_t;

        class MeshView::Impl {
        public:
            ArrayView<idx_t *> elements;
            std::function<void(ArrayView<idx_t *> &)> destroy_elements;

            ArrayView<geom_t *> points;
            std::function<void(ArrayView<geom_t *> &)> destroy_points;

            enum ::smesh::ElemType element_type { ::smesh::INVALID };
            ptrdiff_t nelements{0};
            ptrdiff_t nnodes{0};
            ::smesh::SharedBuffer<idx_t *> element_buffer;
        };

        MeshView::MeshView() : impl_(std::make_shared<Impl>()) {}
        MeshView::~MeshView() = default;

        int MeshView::element_type() const { return impl_->element_type; }

        MeshView::Impl &MeshView::impl() { return *impl_; }

        ArrayView<void *> MeshView::elements() const { return impl_->elements.reinterpret<void *>(); }
        ArrayView<void *> MeshView::points() const { return impl_->points.reinterpret<void *>(); }

        ptrdiff_t MeshView::n_elements() const { return impl_->nelements; }
        ptrdiff_t MeshView::n_nodes() const { return impl_->nnodes; }

        class Mesh::Impl {
        public:
            Communicator comm;
            std::shared_ptr<::smesh::Mesh> mesh;
            mutable std::vector<Mesh::SizeType> node_mapping_cache;
            bool has_distributed{false};
            Impl() : comm(Communicator::get_default()) {}
            ~Impl() {}
        };

        Mesh::Mesh() : impl_(utopia::make_unique<Impl>()) {}
        Mesh::~Mesh() {}

        const Mesh::Communicator &Mesh::comm() const { return impl_->comm; }
        Mesh::Communicator &Mesh::comm() { return impl_->comm; }

        int Mesh::spatial_dimension() const { return impl_->mesh->spatial_dimension(); }

        bool Mesh::read(const Path &path) {
            auto mesh = std::make_shared<::smesh::Mesh>(::smesh::Communicator::wrap(impl_->comm.get()));
            if (SMESH_SUCCESS != mesh->read(::smesh::Path(path.c_str()))) {
                return false;
            }

            impl_->mesh = std::move(mesh);
            impl_->node_mapping_cache.clear();
            impl_->has_distributed = impl_->comm.size() > 1;
            return true;
        }

        void Mesh::read(Input &in) {
            Path path;
            in.require("path", path);
            if (!this->read(path)) {
                utopia::err() << "Unable to read mesh at " << path << "\n";
            }
        }

        void Mesh::create_vector_nodal(Vector &out, int ncomponents) const {
            auto dist = impl_->has_distributed ? impl_->mesh->distributed() : nullptr;
            const auto n_global_nodes = dist ? dist->n_nodes_global() : impl_->mesh->n_nodes();
            out.zeros(layout(impl_->comm, n_local_nodes() * ncomponents, n_global_nodes * ncomponents));
        }

        void Mesh::write_nodal_field(const Path &path, const Vector &field) {
            // TODO
            assert(false);
        }

        bool Mesh::write(const Path &path) {
            UTOPIA_TRACE_SCOPE("sfem::Mesh::write");
            return SMESH_SUCCESS == impl_->mesh->write(::smesh::Path(path.c_str()));
        }

        void Mesh::describe(std::ostream &os) const {
            if (!this->comm().rank()) {
                os << "sfem::Mesh\n";
            }

            os << "n_local_nodes:    " << n_local_nodes() << "\n";
            os << "n_owned_elements: " << n_owned_elements() << "\n";
        }

        std::shared_ptr<::smesh::Mesh> Mesh::raw_type() const { return impl_->mesh; }

        void Mesh::set_raw_type(const std::shared_ptr<::smesh::Mesh> &mesh) {
            impl_->mesh = mesh;
            impl_->node_mapping_cache.clear();
            impl_->has_distributed = false;
        }

        ArrayView<const Mesh::SizeType> Mesh::node_mapping() const {
            const auto n_nodes = impl_->mesh->n_nodes();
            impl_->node_mapping_cache.resize(n_nodes);

            auto mapping = impl_->mesh->node_mapping();
            if (mapping) {
                const auto *data = mapping->data();
                for (ptrdiff_t i = 0; i < n_nodes; ++i) {
                    impl_->node_mapping_cache[i] = data[i];
                }
            } else if (impl_->has_distributed) {
                auto dist = impl_->mesh->distributed();
                auto dist_mapping = dist->node_mapping();
                const auto *data = dist_mapping->data();
                for (ptrdiff_t i = 0; i < n_nodes; ++i) {
                    impl_->node_mapping_cache[i] = data[i];
                }
            } else {
                for (ptrdiff_t i = 0; i < n_nodes; ++i) {
                    impl_->node_mapping_cache[i] = i;
                }
            }

            return ArrayView<const SizeType>(impl_->node_mapping_cache.data(), n_nodes);
        }

        Mesh::SizeType Mesh::n_owned_elements() const {
            auto dist = impl_->has_distributed ? impl_->mesh->distributed() : nullptr;
            return dist ? dist->n_elements_owned() : impl_->mesh->n_elements();
        }

        Mesh::SizeType Mesh::n_local_nodes() const {
            auto dist = impl_->has_distributed ? impl_->mesh->distributed() : nullptr;
            return dist ? dist->n_nodes_owned() : impl_->mesh->n_nodes();
        }

        std::shared_ptr<MeshView> Mesh::sharp_edges(const Scalar angle_threshold) const {
            auto edges = ::smesh::extract_sharp_edges(*impl_->mesh, angle_threshold);

            auto view = std::make_shared<MeshView>();
            view->impl().element_buffer = edges;
            view->impl().elements = ArrayView<idx_t *>(edges->data(), 2);
            view->impl().nelements = edges->extent(1);
            view->impl().destroy_elements = [](ArrayView<idx_t *> &elements) {
                elements = ArrayView<idx_t *>();
            };

            view->impl().points = ArrayView<geom_t *>(impl_->mesh->points()->data(), impl_->mesh->spatial_dimension());
            view->impl().nnodes = impl_->mesh->n_nodes();
            view->impl().destroy_points = [](ArrayView<geom_t *> &points) { points = ArrayView<geom_t *>(); };
            view->impl().element_type = ::smesh::BEAM2;
            return view;
        }

        std::shared_ptr<MeshView> Mesh::disconnected_faces_from_sharp_edges(MeshView &sharp_edges) const {
            auto &se = sharp_edges.impl();

            auto disconnected_elements = ::smesh::extract_disconnected_faces(*impl_->mesh, se.element_buffer);
            ptrdiff_t n_disconnected_elements = disconnected_elements->size();

            if (n_disconnected_elements) {
                int nxe = impl_->mesh->n_nodes_per_element(0);
                auto selected_elements = ::smesh::create_host_buffer<idx_t>(nxe, n_disconnected_elements);
                auto in = impl_->mesh->elements(0)->data();
                auto out = selected_elements->data();
                auto selected = disconnected_elements->data();
                for (ptrdiff_t e = 0; e < n_disconnected_elements; ++e) {
                    const auto idx = selected[e];
                    for (int v = 0; v < nxe; ++v) {
                        out[v][e] = in[v][idx];
                    }
                }

                auto view = std::make_shared<MeshView>();
                view->impl().element_buffer = selected_elements;
                view->impl().elements = ArrayView<idx_t *>(selected_elements->data(), nxe);
                view->impl().nelements = n_disconnected_elements;
                view->impl().destroy_elements = [](ArrayView<idx_t *> &elements) {
                    elements = ArrayView<idx_t *>();
                };

                view->impl().points =
                    ArrayView<geom_t *>(impl_->mesh->points()->data(), impl_->mesh->spatial_dimension());
                view->impl().nnodes = impl_->mesh->n_nodes();
                view->impl().destroy_points = [](ArrayView<geom_t *> &points) { points = ArrayView<geom_t *>(); };

                view->impl().element_type = impl_->mesh->element_type(0);
                return view;

            } else {
                return nullptr;
            }
        }
        std::shared_ptr<MeshView> Mesh::separate_corners_from_sharp_edges(MeshView &sharp_edges,
                                                                          const bool remove_edges) const {
            auto &se = sharp_edges.impl();
            auto edges = se.element_buffer;
            auto corners = ::smesh::extract_sharp_corners(impl_->mesh->n_nodes(), edges, remove_edges);
            se.element_buffer = edges;
            se.elements = ArrayView<idx_t *>(edges->data(), 2);
            se.nelements = edges->extent(1);

            if (corners->size()) {
                auto elements = ::smesh::create_host_buffer<idx_t>(1, corners->size());
                std::copy(corners->data(), corners->data() + corners->size(), elements->data()[0]);

                auto view = std::make_shared<MeshView>();
                view->impl().element_buffer = elements;
                view->impl().elements = ArrayView<idx_t *>(elements->data(), 1);
                view->impl().nelements = corners->size();
                view->impl().destroy_elements = [](ArrayView<idx_t *> &elements) {
                    elements = ArrayView<idx_t *>();
                };

                view->impl().points =
                    ArrayView<geom_t *>(impl_->mesh->points()->data(), impl_->mesh->spatial_dimension());
                view->impl().nnodes = impl_->mesh->n_nodes();
                view->impl().destroy_points = [](ArrayView<geom_t *> &points) { points = ArrayView<geom_t *>(); };
                view->impl().element_type = ::smesh::NODE1;

                return view;
            } else {
                return nullptr;
            }
        }

    }  // namespace sfem
}  // namespace utopia
