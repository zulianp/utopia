#include "utopia_Obstacle.hpp"

#include "moonolith_obstacle.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_make_unique.hpp"

#include "moonolith_affine_transform.hpp"
#include "moonolith_assign_functions.hpp"
#include "moonolith_contact.hpp"
#include "moonolith_elem_quad.hpp"
#include "moonolith_elem_segment.hpp"
#include "moonolith_elem_shape.hpp"
#include "moonolith_elem_triangle.hpp"
#include "moonolith_matlab_scripter.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_sparse_matrix.hpp"

#include "utopia_LibMeshToMoonolithConvertions.hpp"
#include "utopia_TransferUtils.hpp"

namespace utopia {

    void ObstacleParams::read(Input &in) {
        in.get("variable_number", variable_number);
        in.get("gap_negative_bound", gap_negative_bound);
        in.get("gap_positive_bound", gap_positive_bound);
    }

    void Obstacle::transform(const USparseMatrix &in, USparseMatrix &out) {
        out = transpose(output().orthogonal_trafo) * in * output().orthogonal_trafo;
    }

    void Obstacle::transform(const UVector &in, UVector &out) { out = output().orthogonal_trafo * in; }

    void Obstacle::inverse_transform(const UVector &in, UVector &out) { out = output().orthogonal_trafo * in; }

    class Obstacle::Impl : public Configurable {
    public:
        virtual ~Impl() {}
        Impl() {}

        virtual void init(libMesh::MeshBase &mesh) = 0;

        virtual bool assemble(libMesh::MeshBase &mesh,
                              libMesh::DofMap &dof_map,
                              const ObstacleParams &params,
                              Obstacle::Output &output) = 0;

        virtual void read(Input &) override {}
    };

    template <int Dim>
    class ImplD final : public Obstacle::Impl {
    public:
        using MeshT = moonolith::Mesh<double, Dim>;
        using SpaceT = moonolith::FunctionSpace<MeshT>;
        using ObstacleT = moonolith::Obstacle<double, Dim>;

        ImplD() {}
        ~ImplD() override {}

        bool assemble(libMesh::MeshBase &mesh,
                      libMesh::DofMap &dof_map,
                      const ObstacleParams &params,
                      Obstacle::Output &output) override {
            moonolith::Communicator comm = mesh.comm().get();

            auto mesh_ptr = std::make_shared<MeshT>(comm);
            SpaceT space(mesh_ptr);
            extract_trace_space(mesh, dof_map, params.variable_number, space);

            obstacle->set_gap_bounds(params.gap_negative_bound, params.gap_positive_bound);
            obstacle->assemble({}, space);

            finalize_tensors(output);
            return true;
        }

        void init(libMesh::MeshBase &mesh) override {
            moonolith::Communicator comm = mesh.comm().get();
            obstacle_mesh = utopia::make_unique<MeshT>(comm);

            if (mesh.mesh_dimension() == Dim - 1) {
                convert(mesh, *obstacle_mesh);
            } else {
                assert(mesh.mesh_dimension() == Dim);
                extract_surface<Dim>(mesh, *obstacle_mesh, {});
            }

            auto tree = moonolith::RayCastingTree<double, Dim>::New();
            tree->init(*obstacle_mesh);
            obstacle = utopia::make_unique<ObstacleT>(tree);
        }

        void finalize_tensors(Obstacle::Output &out) {
            auto &buffers = obstacle->buffers();

            USparseMatrix mass_matrix_x, trafo_x;  //, inverse_trafo_x;
            UVector gap_x;

            convert_matrix(buffers.mass_matrix, mass_matrix_x);
            convert_matrix(buffers.trafo, trafo_x);
            tensorize(trafo_x, Dim, out.basis_trafo);
            // convert_matrix(buffers.inverse_trafo, inverse_trafo_x);

            convert_tensor(buffers.gap, gap_x);
            convert_tensor(buffers.normal, out.normals);
            convert_tensor(buffers.is_contact, out.is_contact);

            UVector mass_vector = sum(mass_matrix_x, 1);

            tensorize(Dim, mass_vector);
            e_pseudo_inv(mass_vector, out.inverse_mass_vector, 1e-15);

            // FIXME
            out.gap = e_mul(out.inverse_mass_vector, gap_x);
            out.normals = e_mul(out.inverse_mass_vector, out.normals);

            normalize(out.normals);

            build_orthogonal_transformation(out.is_contact, out.normals, out.orthogonal_trafo);

            // UVector test_normals = out.orthogonal_trafo * out.normals;

            write("gap.m", out.gap);

            replace_zeros(out.is_contact, out.gap);
        }

        void replace_zeros(const UVector &is_contact, UVector &gap) {
            auto r = range(gap);

            Read<UVector> r_ic(is_contact);
            Write<UVector> w_g(gap);

            for (auto i = r.begin(); i < r.end(); ++i) {
                if (is_contact.get(i) == 0.0) {
                    gap.set(i, 1e8);
                }
            }
        }

        void normalize(UVector &normal) {
            ReadAndWrite<UVector> rw_normal(normal);
            auto r = normal.range();

            moonolith::Vector<double, Dim> n;
            for (auto i = r.begin(); i < r.end(); i += Dim) {
                for (int d = 0; d < Dim; ++d) {
                    n[d] = normal.get(i + d);
                }

                auto len_n = length(n);

                if (len_n == 0.0) continue;

                n /= len_n;

                for (int d = 0; d < Dim; ++d) {
                    normal.set(i + d, n[d]);
                }
            }
        }

        void build_orthogonal_transformation(const UVector &is_contact, const UVector &normal, USparseMatrix &trafo) {
            trafo.sparse(square_matrix_layout(layout(normal)), Dim, 0);

            moonolith::HouseholderTransformation<double, Dim> H;
            auto r = range(normal);

            Read<UVector> r_normal(normal), r_ic(is_contact);
            Write<USparseMatrix> w_ot(trafo, utopia::LOCAL);

            moonolith::Vector<double, Dim> n;
            for (auto i = r.begin(); i < r.end(); i += Dim) {
                if (is_contact.get(i) == 0.0) {
                    H.identity();
                } else {
                    for (int d = 0; d < Dim; ++d) {
                        n[d] = normal.get(i + d);
                    }

                    auto len = length(n);
                    assert(len > 0.0);

                    n /= len;
                    n.x -= 1.0;

                    len = length(n);

                    if (approxeq(len, 0.0, 1e-15)) {
                        H.identity();
                    } else {
                        n /= len;
                        H.init(n);
                    }
                }

                assert(approxeq(std::abs(measure(H)), 1.0, 1e-8));

                for (int d1 = 0; d1 < Dim; ++d1) {
                    for (int d2 = 0; d2 < Dim; ++d2) {
                        trafo.set(i + d1, i + d2, H(d1, d2));
                    }
                }
            }
        }

        std::unique_ptr<ObstacleT> obstacle;
        std::unique_ptr<MeshT> obstacle_mesh;
    };

    Obstacle::Obstacle() {}
    Obstacle::~Obstacle() {}

    bool Obstacle::assemble(libMesh::MeshBase &mesh, libMesh::DofMap &dof_map, const ObstacleParams &params) {
        if (!impl_) return false;

        return impl_->assemble(mesh, dof_map, params, output_);
    }

    bool Obstacle::init(libMesh::MeshBase &obstacle_mesh) {
        if (obstacle_mesh.spatial_dimension() == 3) {
            impl_ = utopia::make_unique<ImplD<3>>();
            impl_->init(obstacle_mesh);
            return true;
        } else {
            return false;
        }
    }

}  // namespace utopia
