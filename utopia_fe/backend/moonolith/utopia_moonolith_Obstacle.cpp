#include "utopia_moonolith_Obstacle.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_moonolith_ConvertTensor.hpp"

#include "par_moonolith_instance.hpp"

#include "moonolith_affine_transform.hpp"
#include "moonolith_assign_functions.hpp"
#include "moonolith_contact.hpp"
#include "moonolith_elem_quad.hpp"
#include "moonolith_elem_segment.hpp"
#include "moonolith_elem_shape.hpp"
#include "moonolith_elem_triangle.hpp"
#include "moonolith_matlab_scripter.hpp"
#include "moonolith_obstacle.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_sparse_matrix.hpp"

#include <unordered_set>

namespace utopia {
    namespace moonolith {

        void Obstacle::Params::read(Input &in) {
            in.get("variable_number", variable_number);
            in.get("gap_negative_bound", gap_negative_bound);
            in.get("gap_positive_bound", gap_positive_bound);

            in.get("surfaces", [this](Input &in) {
                in.get_all([this](Input &in) {
                    int tag = -1;
                    in.get("tag", tag);
                    if (tag >= 0) {
                        this->tags.insert(tag);
                    }
                });
            });
        }

        class Obstacle::Output {
        public:
            Vector gap;
            Vector normals;
            Matrix orthogonal_trafo;

            Matrix basis_trafo;
            Vector mass_vector;
            Vector inverse_mass_vector;
            Vector is_contact;
        };

        class Obstacle::Impl {
        public:
            virtual ~Impl() = default;
            Impl() = default;

            inline static void tensorize(const Matrix &T_x, const SizeType n_var, Matrix &T) {
                auto max_nnz = utopia::max_row_nnz(T_x);
                T.sparse(layout(T_x), max_nnz, max_nnz);

                assert(!empty(T));
                assert(T.row_range().extent() % n_var == 0);

                Write<Matrix> w(T);
                each_read(T_x, [&](const SizeType i, const SizeType j, const Scalar value) {
                    for (SizeType k = 0; k < n_var; ++k) {
                        T.set(i + k, j + k, value);
                    }
                });
            }

            inline static void tensorize(const SizeType n_var, Vector &t) {
                ReadAndWrite<Vector> w(t);
                auto r = range(t);

                assert(!empty(t));
                assert(r.extent() % n_var == 0);

                for (auto i = r.begin(); i < r.end(); i += n_var) {
                    const auto value = t.get(i);

                    for (SizeType k = 1; k < n_var; ++k) {
                        t.set(i + k, value);
                    }
                }
            }

            template <class Buffers>
            void finalize_tensors(int dim, Buffers &buffers, Output &out) {
                Matrix mass_matrix_x, trafo_x;
                Vector gap_x;

                convert(buffers.mass_matrix.get(), mass_matrix_x);
                convert(buffers.trafo.get(), trafo_x);
                tensorize(trafo_x, dim, out.basis_trafo);

                convert(buffers.gap.get(), gap_x);
                convert(buffers.normal.get(), out.normals);
                convert(buffers.is_contact.get(), out.is_contact);

                assert(gap_x.local_size() == out.normals.local_size());
                assert(gap_x.local_size() == out.is_contact.local_size());

                Vector mass_vector = sum(mass_matrix_x, 1);

                assert(gap_x.local_size() == mass_vector.local_size());

                tensorize(dim, mass_vector);
                e_pseudo_inv(mass_vector, out.inverse_mass_vector, 1e-15);

                // FIXME
                out.gap = e_mul(out.inverse_mass_vector, gap_x);
                out.normals = e_mul(out.inverse_mass_vector, out.normals);

                normalize(out.normals);

                build_orthogonal_transformation(out.is_contact, out.normals, out.orthogonal_trafo);
                replace_zeros(out.is_contact, out.gap);
            }

            void replace_zeros(const Vector &is_contact, Vector &gap) {
                auto r = range(gap);

                Read<Vector> r_ic(is_contact);
                Write<Vector> w_g(gap);

                for (auto i = r.begin(); i < r.end(); ++i) {
                    if (is_contact.get(i) == 0.0) {
                        gap.set(i, 1e8);
                    }
                }
            }

            virtual bool init(const Mesh &mesh) = 0;
            virtual bool assemble(const FunctionSpace &space, const Params &params, Output &output) = 0;
            virtual void normalize(Vector &normal) = 0;
            virtual void build_orthogonal_transformation(const Vector &is_contact,
                                                         const Vector &normal,
                                                         Matrix &trafo) = 0;
        };

        template <int Dim>
        class Obstacle::ImplD : public Obstacle::Impl {
        public:
            using MoonolithMesh_t = ::moonolith::Mesh<Scalar, Dim>;
            using MoonolithSpace_t = ::moonolith::FunctionSpace<MoonolithMesh_t>;
            using MoonolithObstacle_t = ::moonolith::Obstacle<Scalar, Dim>;

            void normalize(Vector &normal) override {
                ReadAndWrite<Vector> rw_normal(normal);
                auto r = normal.range();

                ::moonolith::Vector<Scalar, Dim> n;
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

            void build_orthogonal_transformation(const Vector &is_contact,
                                                 const Vector &normal,
                                                 Matrix &trafo) override {
                trafo.sparse(square_matrix_layout(layout(normal)), Dim, 0);

                ::moonolith::HouseholderTransformation<Scalar, Dim> H;
                auto r = range(normal);

                Read<Vector> r_normal(normal), r_ic(is_contact);
                Write<Matrix> w_ot(trafo, utopia::LOCAL);

                ::moonolith::Vector<Scalar, Dim> n;
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

            bool assemble(const FunctionSpace &space, const Params &params, Output &output) override {
                std::shared_ptr<MoonolithSpace_t> space_d;
                if (space.mesh().manifold_dimension() == Dim - 1) {
                    space_d = space.raw_type<Dim>();
                } else {
                    assert(false && "IMPLEMENT ME");
                    return false;
                }

                obstacle->set_gap_bounds(params.gap_negative_bound, params.gap_positive_bound);
                obstacle->assemble(params.tags, *space_d);

                this->finalize_tensors(Dim, obstacle->buffers(), output);
                return true;
            }

            bool init(const Mesh &mesh) override {
                assert(mesh.spatial_dimension() == Dim);

                auto raw_mesh = mesh.raw_type<Dim>();
                assert(raw_mesh);
                if (!raw_mesh) return false;

                if (mesh.manifold_dimension() == Dim - 1) {
                    obstacle_mesh = raw_mesh;
                } else {
                    auto mesh_ptr = std::make_shared<MoonolithMesh_t>(raw_mesh->comm());
                    assert(false && "IMPLEMENT ME");
                    return false;
                }

                auto tree = ::moonolith::RayCastingTree<Scalar, Dim>::New();
                tree->init(*obstacle_mesh);
                obstacle = utopia::make_unique<MoonolithObstacle_t>(tree);

                return true;
            }

            std::unique_ptr<MoonolithObstacle_t> obstacle;
            std::shared_ptr<MoonolithMesh_t> obstacle_mesh;
        };

        Obstacle::Obstacle() : params_(utopia::make_unique<Params>()), output_(utopia::make_unique<Output>()) {}

        Obstacle::~Obstacle() = default;

        void Obstacle::read(Input &in) { params_->read(in); }

        void Obstacle::describe(std::ostream &os) const {}

        bool Obstacle::init_obstacle(const Mesh &obstacle_mesh) {
            if (obstacle_mesh.spatial_dimension() == 3) {
                impl_ = utopia::make_unique<ImplD<3>>();
                impl_->init(obstacle_mesh);
                return true;
            } else {
                return false;
            }
        }

        void Obstacle::set_params(const Params &params) { *this->params_ = params; }

        bool Obstacle::assemble(const FunctionSpace &space) { return impl_->assemble(space, *params_, *output_); }

        void Obstacle::transform(const Matrix &in, Matrix &out) {
            out = transpose(output().orthogonal_trafo) * in * output().orthogonal_trafo;
        }

        void Obstacle::transform(const Vector &in, Vector &out) { out = transpose(output().orthogonal_trafo) * in; }

        void Obstacle::inverse_transform(const Vector &in, Vector &out) { out = output().orthogonal_trafo * in; }

        const Obstacle::Vector &Obstacle::gap() const { return output().gap; }

        const Obstacle::Vector &Obstacle::is_contact() const { return output().is_contact; }

        const Obstacle::Vector &Obstacle::normals() const { return output().normals; }

        Obstacle::Output &Obstacle::output() { return *output_; }
        const Obstacle::Output &Obstacle::output() const { return *output_; }

    }  // namespace moonolith
}  // namespace utopia
