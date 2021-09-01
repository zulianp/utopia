#include "utopia_petsc_StructuredGrid.hpp"
#include "utopia_mesh_StructuredGrid.hpp"

#include "utopia_petsc_DMDABase.hpp"

namespace utopia {
    namespace petsc {

        StructuredGrid::StructuredGrid(const Communicator &comm) : impl_(utopia::make_unique<DMDABase>(comm)) {}
        StructuredGrid::~StructuredGrid() = default;

        void StructuredGrid::unit_cube(const SizeType nx, const SizeType ny, const SizeType nz) {
            // Add one to be conforming to the other backends
            impl_->unit_cube(nx + 1, ny + 1, nz + 1);
        }

        StructuredGrid::Communicator &StructuredGrid::comm() { return impl_->comm(); }
        const StructuredGrid::Communicator &StructuredGrid::comm() const { return impl_->comm(); }

        std::unique_ptr<StructuredGrid> StructuredGrid::uniform_refine() const {
            auto fine = utopia::make_unique<StructuredGrid>(comm());
            fine->impl_ = std::move(impl_->uniform_refine());
            return fine;
        }

        std::unique_ptr<StructuredGrid> StructuredGrid::clone(const SizeType &n_components) const {
            auto cloned = utopia::make_unique<StructuredGrid>(comm());
            cloned->impl_ = std::move(impl_->clone());
            // cloned->impl_->set_dof(n_components);
            return cloned;
        }

        std::unique_ptr<StructuredGrid> StructuredGrid::clone() const { return clone(this->impl_->get_dof()); }

        void StructuredGrid::read(Input &in) { impl_->read(in); }

        void StructuredGrid::set_field_name(const SizeType &nf, const std::string &name) {
            impl_->set_field_name(nf, name);
        }

        void StructuredGrid::set_field_names(const std::vector<std::string> &names) { impl_->set_field_names(names); }

        // void StructuredGrid::nodes(const SizeType &idx, NodeIndex &nodes) const {
        //     assert(elements_ && "init_elements() has to be called before any access to element data");
        //     nodes = impl_->elements_->nodes(idx);
        // }

        // void StructuredGrid::nodes_local(const SizeType &idx, NodeIndex &nodes) const {
        //     assert(elements_ && "init_elements() has to be called before any access to element data");
        //     nodes = impl_->elements_->nodes_local(idx);
        // }

        // bool StructuredGrid::on_boundary(const SizeType &elem_idx) const {
        //     NodeIndex idx;
        //     nodes(elem_idx, idx);

        //     for (auto i : idx) {
        //         if (this->impl_->is_node_on_boundary(i)) {
        //             return true;
        //         }
        //     }

        //     return false;
        // }

        StructuredGrid::View StructuredGrid::view() const {
            SizeType dof_range_begin, dof_range_end;
            IntArray dims, corners_begin, corners_extent, ghost_corners_begin, ghost_corners_extent;

            impl_->dof_ownership_range(dof_range_begin, dof_range_end);
            const SizeType dim = impl_->get_dimension();
            impl_->get_dims(dims);
            impl_->get_corners(corners_begin, corners_extent);
            impl_->get_ghost_corners(ghost_corners_begin, ghost_corners_extent);
            SizeType n_components = impl_->get_dof();
            auto &box_min = impl_->box_min();
            auto &box_max = impl_->box_max();

            std::vector<Scalar> v_box_min(dim), v_box_max(dim);
            std::vector<SizeType> v_dims(dim), v_corners_begin(dim), v_corners_extent(dim), v_ghost_corners_begin(dim),
                v_ghost_corners_extent(dim);

            for (SizeType d = 0; d < dim; ++d) {
                v_dims[d] = dims[d];
                v_corners_begin[d] = corners_begin[d];
                v_corners_extent[d] = corners_extent[d];
                v_ghost_corners_begin[d] = ghost_corners_begin[d];
                v_ghost_corners_extent[d] = ghost_corners_extent[d];
                v_box_min[d] = box_min[d];
                v_box_max[d] = box_max[d];
            }

            View v(v_dims,
                   v_corners_begin,
                   v_corners_extent,
                   v_ghost_corners_begin,
                   v_ghost_corners_extent,
                   v_box_min,
                   v_box_max,
                   n_components,
                   impl_->elements_x_cell());

            v.set_dof_range_begin(dof_range_begin);
            v.set_dof_range_end(dof_range_end);
            return v;
        }

        void StructuredGrid::describe(std::ostream &os) const {
            auto v = view();
            os << "n_elements      : " << v.n_elements() << '\n';
            os << "n_nodes         : " << v.n_nodes() << '\n';
            os << "dim             : " << v.dim() << '\n';
            os << "elements_x_cell : " << v.elements_x_cell() << '\n';
        }

        StructuredGrid::SizeType StructuredGrid::n_nodes() const {
            IntArray dims;
            impl_->get_dims(dims);

            SizeType ret = 1;
            for (auto d : dims) {
                if (d > 0) ret *= d;
            }
            return ret;
        }

        DMDABase &StructuredGrid::dm() const { return *impl_; }

        StructuredGrid::ElemToNodeIndex StructuredGrid::elem_to_node_index() const {
            return impl_->elem_to_node_index();
        }

        StructuredGrid::ElemToNodeIndex StructuredGrid::elem_to_local_node_index() const {
            return impl_->elem_to_local_node_index();
        }

    }  // namespace petsc
}  // namespace utopia