#ifndef UTOPIA_S_ELEMENT_ADAPTER_HPP
#define UTOPIA_S_ELEMENT_ADAPTER_HPP

#include "utopia_BoxAdapter.hpp"

#include "MortarAssemble.hpp"
#include "libmesh/serial_mesh.h"
#include "moonolith_input_stream.hpp"
#include "moonolith_output_stream.hpp"
#include "moonolith_serializable.hpp"

namespace utopia {

    template <int Dimension>
    class SElementAdapter : public moonolith::Serializable {
    public:
        inline int tag() const { return tag_; }

        const BoxBoxAdapter<Dimension> &bound() const { return bound_; }

        BoxBoxAdapter<Dimension> &bound() { return bound_; }

        void applyRW(moonolith::Stream &stream) {
            stream &bound_;
            stream &element_;
            stream &element_handle_;
        }

        SElementAdapter(libMesh::MeshBase &fe,
                        const libMesh::dof_id_type &element,
                        const long element_handle,
                        const int tag,
                        /*std::vector<long> &map,*/ const libMesh::Real blow_up)
            : fe_(&fe),
              element_(element),
              element_handle_(element_handle),
              tag_(tag),
              dof_map_(nullptr),
              face_id_(nullptr) {
            assert(element < fe.n_elem());

            libMesh::Elem *e = fe.elem(element);

            const int dim = fe.mesh_dimension();

            libMesh::Point o, u, v, n, c, p;

            for (uint side = 0; side < e->n_sides(); ++side) {
                libMesh::UniquePtr<libMesh::Elem> s = e->build_side_ptr(side);

                if (!s->on_boundary()) continue;

                auto side_ptr = e->build_side_ptr(side);

                compute_side_normal(dim, *side_ptr, n);

                if (fix_normal_orientation(*e, side, n)) {
                    std::cerr << "[Warning] face with wrong orientation detected\n" << std::endl;
                }

                assert(n.contract(side_ptr->centroid() - e->centroid()) > 0);

                for (uint i = 0; i < side_ptr->n_nodes(); ++i) {
                    c = side_ptr->point(i);

                    n *= blow_up;
                    p = c;
                    p += n;

                    std::array<double, Dimension> p_a;
                    for (int d = 0; d < Dimension; ++d) {
                        p_a[d] = p(d);
                    }

                    bound_.static_bound() += p_a;
                    bound_.dynamic_bound() += p_a;
                    p = c;
                    n *= 0.5;
                    p -= n;
                    bound_.static_bound() += p_a;
                    bound_.dynamic_bound() += p_a;
                }
            }
        }

        SElementAdapter()
            : fe_(nullptr), element_(-1), element_handle_(-1), tag_(-1), dof_map_(nullptr), face_id_(nullptr) {}

        inline long handle() const { return element_handle_; }

        inline long element() const { return element_; }

        libMesh::Elem *get() {
            assert(fe_);

            // std::cout<<"I AM IN GET"<<std::endl;

            assert(element_ < fe_->n_local_elem());

            return fe_->elem(element_);
        }

        const libMesh::Elem *get() const {
            assert(fe_);
            return fe_->elem(element_);
        }

        inline const libMesh::MeshBase &space() const {
            assert(fe_);
            return *fe_;
        }

        void set_dof_map(std::vector<long> *ptr) { dof_map_ = ptr; }

        void set_face_id(std::vector<long> *ptr2) { face_id_ = ptr2; }

        inline const std::vector<long> &dof_map() const {
            assert(dof_map_);
            return *dof_map_;
        }

        inline const std::vector<long> &dof_map_face() const {
            assert(face_id_);
            return *face_id_;
        }

    private:
        libMesh::MeshBase *fe_;
        libMesh::dof_id_type element_;
        long element_handle_;
        int tag_;
        BoxBoxAdapter<Dimension> bound_;
        std::vector<long> *dof_map_;
        std::vector<long> *face_id_;
    };
}  // namespace utopia

#endif  // UTOPIA_S_ELEMENT_ADAPTER_HPP
