#include "utopia_sfem_Mesh.hpp"

#include "utopia_IOStream.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_make_unique.hpp"

#include "read_mesh.h"
#include "sfem_base.h"
#include "sfem_defs.h"
#include "sfem_mesh_write.h"

#include "crs_graph.h"
#include "extract_sharp_features.h"

#include "mesh_utils.h"

#include <cassert>
#include <cstddef>
#include <functional>
#include <utopia_ArrayView.hpp>

namespace utopia {
    namespace sfem {
        class MeshView::Impl {
        public:
            ArrayView<idx_t *> elements;
            std::function<void(ArrayView<idx_t *> &)> destroy_elements;

            ArrayView<geom_t *> points;
            std::function<void(ArrayView<geom_t *> &)> destroy_points;

            enum ElemType element_type { INVALID };
            ptrdiff_t nelements{0};
            ptrdiff_t nnodes{0};
        };

        MeshView::MeshView() : impl_(std::make_shared<Impl>()) {}
        MeshView::~MeshView() = default;

        int MeshView::element_type() const
        {
            return impl_->element_type;
        }

        MeshView::Impl &MeshView::impl() { return *impl_; }

        ArrayView<void *> MeshView::elements() const
        {
            return impl_->elements.reinterpret<void*>();
        }
        ArrayView<void *> MeshView::points() const
        {
            return impl_->points.reinterpret<void*>();
        }

        ptrdiff_t MeshView::n_elements() const
        {
            return impl_->nelements;
        }
        ptrdiff_t MeshView::n_nodes() const
        {
            return impl_->nnodes;
        }

        class Mesh::Impl {
        public:
            Communicator comm;
            mesh_t mesh;
            Impl() : comm(Communicator::get_default()) { mesh_init(&mesh); }
            ~Impl() { mesh_destroy(&mesh); }
        };

        Mesh::Mesh() : impl_(utopia::make_unique<Impl>()) {}
        Mesh::~Mesh() {}

        const Mesh::Communicator &Mesh::comm() const { return impl_->comm; }
        Mesh::Communicator &Mesh::comm() { return impl_->comm; }

        int Mesh::spatial_dimension() const { return impl_->mesh.spatial_dim; }

        bool Mesh::read(const Path &path) {
            return SFEM_OK == mesh_read(impl_->comm.get(), path.c_str(), &impl_->mesh);
        }

        void Mesh::read(Input &in) {
            Path path;
            in.require("path", path);
            if (!this->read(path)) {
                utopia::err() << "Unable to read mesh at " << path << "\n";
            }
        }

        void Mesh::create_vector_nodal(Vector &out, int ncomponents) const {
            out.zeros(layout(impl_->comm, impl_->mesh.n_owned_nodes * ncomponents, impl_->mesh.nnodes * ncomponents));
        }

        void Mesh::write_nodal_field(const Path &path, const Vector &field) {
            // TODO
            assert(false);
        }

        bool Mesh::write(const Path &path) {
            UTOPIA_TRACE_SCOPE("sfem::Mesh::write");
            return SFEM_OK == mesh_write(path.c_str(), &impl_->mesh);
        }

        void Mesh::describe(std::ostream &os) const {
            // TODO
        }

        void *Mesh::raw_type() const { return (void *)&impl_->mesh; }

        ArrayView<const Mesh::SizeType> Mesh::node_mapping() const {
            return ArrayView<const SizeType>(impl_->mesh.node_mapping, impl_->mesh.nnodes);
        }

        Mesh::SizeType Mesh::n_local_nodes() const { return impl_->mesh.n_owned_nodes; }

        std::shared_ptr<MeshView> Mesh::sharp_edges(const Scalar angle_threshold) const {
            ptrdiff_t n_sharp_edges = 0;
            idx_t *e0 = 0;
            idx_t *e1 = 0;

            {
                count_t *rowptr = 0;
                idx_t *colidx = 0;
                build_crs_graph_for_elem_type(impl_->mesh.element_type,
                                              impl_->mesh.nelements,
                                              impl_->mesh.nnodes,
                                              impl_->mesh.elements,
                                              &rowptr,
                                              &colidx);

                extract_sharp_edges((enum ElemType)impl_->mesh.element_type,
                                    impl_->mesh.nelements,
                                    impl_->mesh.nnodes,
                                    impl_->mesh.elements,
                                    impl_->mesh.points,
                                    // CRS-graph (node to node)
                                    rowptr,
                                    colidx,
                                    angle_threshold,
                                    &n_sharp_edges,
                                    &e0,
                                    &e1);

                free(rowptr);
                free(colidx);
            }

            idx_t **elements = (idx_t **)malloc(2 * sizeof(idx_t *));
            elements[0] = e0;
            elements[1] = e1;

            auto view = std::make_shared<MeshView>();
            view->impl().elements = ArrayView<idx_t *>(elements, 2);
            view->impl().nelements = n_sharp_edges;
            view->impl().destroy_elements = [](ArrayView<idx_t *> &elements) {
                free(elements[0]);
                free(elements[1]);
                free(elements.begin());
                elements = ArrayView<idx_t *>();
            };

            view->impl().points = ArrayView<geom_t *>(impl_->mesh.points, impl_->mesh.spatial_dim);
            view->impl().nnodes = impl_->mesh.nnodes;
            view->impl().destroy_points = [](ArrayView<geom_t *> &points) { points = ArrayView<geom_t *>(); };
            view->impl().element_type = BEAM2;
            return view;
        }

        std::shared_ptr<MeshView> Mesh::disconnected_faces_from_sharp_edges(MeshView &sharp_edges) const {
            auto &se = sharp_edges.impl();

            ptrdiff_t n_disconnected_elements = 0;
            element_idx_t *disconnected_elements = 0;

            extract_disconnected_faces((enum ElemType)impl_->mesh.element_type,
                                       impl_->mesh.nelements,
                                       impl_->mesh.nnodes,
                                       impl_->mesh.elements,
                                       se.nelements,
                                       se.elements[0],
                                       se.elements[1],
                                       &n_disconnected_elements,
                                       &disconnected_elements);

            if (disconnected_elements) {
                int nxe = elem_num_nodes((enum ElemType)impl_->mesh.element_type);
                idx_t **selected_elements = allocate_elements(nxe, n_disconnected_elements);
                select_elements(
                    nxe, n_disconnected_elements, disconnected_elements, impl_->mesh.elements, selected_elements);

                auto view = std::make_shared<MeshView>();
                view->impl().elements = ArrayView<idx_t *>(selected_elements, nxe);
                view->impl().nelements = n_disconnected_elements;
                view->impl().destroy_elements = [](ArrayView<idx_t *> &elements) {
                    free_elements(elements.size(), elements.begin());
                    elements = ArrayView<idx_t *>();
                };

                view->impl().points = ArrayView<geom_t *>(impl_->mesh.points, impl_->mesh.spatial_dim);
                view->impl().nnodes = impl_->mesh.nnodes;
                view->impl().destroy_points = [](ArrayView<geom_t *> &points) { points = ArrayView<geom_t *>(); };

                view->impl().element_type = (enum ElemType)impl_->mesh.element_type;
                return view;

            } else {
                free(disconnected_elements);
                return nullptr;
            }
        }
        std::shared_ptr<MeshView> Mesh::separate_corners_from_sharp_edges(MeshView &sharp_edges,
                                                                          const bool remove_edges) const {
            auto &se = sharp_edges.impl();
            ptrdiff_t n_corners = 0;
            idx_t *corners = nullptr;

            se.nelements = extract_sharp_corners(
                impl_->mesh.nnodes, se.nelements, se.elements[0], se.elements[1], &n_corners, &corners, remove_edges);

            if (n_corners) {
                idx_t **elements = (idx_t **)malloc(sizeof(idx_t *));
                elements[0] = corners;

                auto view = std::make_shared<MeshView>();
                view->impl().elements = ArrayView<idx_t *>(elements, 1);
                view->impl().nelements = n_corners;
                view->impl().destroy_elements = [](ArrayView<idx_t *> &elements) {
                    free_elements(elements.size(), elements.begin());
                    elements = ArrayView<idx_t *>();
                };

                view->impl().points = ArrayView<geom_t *>(impl_->mesh.points, impl_->mesh.spatial_dim);
                view->impl().nnodes = impl_->mesh.nnodes;
                view->impl().destroy_points = [](ArrayView<geom_t *> &points) { points = ArrayView<geom_t *>(); };
                view->impl().element_type = NODE1;

                return view;
            } else {
                return nullptr;
            }
        }

    }  // namespace sfem
}  // namespace utopia
