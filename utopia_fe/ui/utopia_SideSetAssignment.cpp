#include "utopia_SideSetAssignment.hpp"
#include "utopia_Input.hpp"
#include "utopia_fe_base.hpp"

#include "MortarAssemble.hpp"
#include "utopia_SideSetAssignment.hpp"

#include "libmesh/elem.h"
#include "libmesh/parallel_mesh.h"

#include <cmath>

namespace utopia {

    using MeshT = libMesh::UnstructuredMesh;

    template <class Mesh>
    class SideSetAssignment<Mesh>::Impl : public Configurable {
    public:
        class Normal2SideSet : public Configurable {
        public:
            void read(Input &in) override {
                in.get("selection", selection);
                double nx = 0, ny = 0, nz = 0;

                in.get("nx", nx);
                in.get("ny", ny);
                in.get("nz", nz);

                n(0) = nx;
                n(1) = ny;
                n(2) = nz;

                in.get("side-set", side_set);

                in.get("tol", tol);

                in.get("block", block);
            }

            inline bool valid() const { return side_set != -1; }

            void describe(std::ostream &os = std::cout) const {
                os << "selection " << selection << std::endl;
                os << "nx        " << n(0) << std::endl;
                os << "ny        " << n(1) << std::endl;
                os << "nz        " << n(2) << std::endl;
                os << "side-set  " << side_set << std::endl;
                os << "block     " << block << std::endl;
            }

            Normal2SideSet() : selection(-1), n(), side_set(-1), tol(1e-8), block(-1) {}

            int selection;
            libMesh::Point n;
            int side_set;
            double tol;
            int block;
        };

        void read(Input &in) override {
            in.get_all([this](Input &in) {
                Normal2SideSet n;
                n.read(in);

                n.describe(std::cout);

                if (n.valid()) {
                    n2ss.push_back(n);
                } else {
                    std::cerr << "[Error] bad format: " << std::endl;
                    n.describe(std::cerr);
                }
            });
        }

        void apply(Mesh &mesh) {
            if (n2ss.empty()) {
                return;
            }

            const int dim = mesh.mesh_dimension();
            auto &bi = mesh.get_boundary_info();

            libMesh::Point normal;
            for (const auto &elem_ptr : mesh.active_local_element_ptr_range()) {
                const std::size_t n_sides = elem_ptr->n_sides();

                for (std::size_t i = 0; i < n_sides; ++i) {
                    if ((elem_ptr->neighbor_ptr(i) != libmesh_nullptr)) {
                        continue;
                    }

                    int side_set = bi.boundary_id(elem_ptr, i);

                    auto side_ptr = elem_ptr->build_side_ptr(i);

                    compute_side_normal(dim, *side_ptr, normal);

                    if (side_set < 0) {
                        side_set = -1;
                    }

                    for (const auto &n : n2ss) {
                        if (n.selection == side_set) {
                            if (n.block != -1 && n.block != elem_ptr->subdomain_id()) continue;

                            const double angle = normal * n.n;

                            if (std::abs(angle - 1) < n.tol) {
                                if (side_set != -1) {
                                    bi.remove_side(elem_ptr, i);
                                }

                                bi.add_side(elem_ptr, i, n.side_set);
                                assert(bi.boundary_id(elem_ptr, i) == n.side_set);

                                std::cout << "changed " << side_set << " to " << n.side_set << std::endl;
                            }
                        }
                    }
                }
            }
        }

    private:
        std::vector<Normal2SideSet> n2ss;
    };

    template <class Mesh>
    void SideSetAssignment<Mesh>::read(Input &in) {
        impl_->read(in);
    }

    template <class Mesh>
    void SideSetAssignment<Mesh>::apply(Mesh &mesh) {
        impl_->apply(mesh);
    }

    template <class Mesh>
    SideSetAssignment<Mesh>::SideSetAssignment() : impl_(utopia::make_unique<Impl>()) {}

    template <class Mesh>
    SideSetAssignment<Mesh>::~SideSetAssignment() {}

    template class SideSetAssignment<MeshT>;

}  // namespace utopia
