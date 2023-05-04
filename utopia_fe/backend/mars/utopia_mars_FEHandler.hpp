#ifndef UTOPIA_MARS_FEHANDLER_HPP
#define UTOPIA_MARS_FEHANDLER_HPP

#include "mars.hpp"
#include "mars_base.hpp"
#include "mars_context.hpp"
#include "mars_distributed_dof_management.hpp"
#include "mars_globals.hpp"
#include "utopia_mars_FunctionSpace.hpp"

#include "utopia_mars_Factory_impl.hpp"

namespace utopia {
    namespace mars {

        template <class DMesh, int Degree_ = 1>
        class FEHandler : public IFEHandler {
        public:
            static constexpr ::mars::Integer Degree = Degree_;
            static constexpr int Dim = DMesh::Dim;
            using DofHandler = ::mars::DofHandler<DMesh, Degree, 0>;
            // using DofHandler = ::mars::DofHandler<DMesh, Degree, 1>;
            using FEDofMap = ::mars::FEDofMap<DofHandler>;
            using SPattern =
                ::mars::SparsityPattern<Scalar, LocalSizeType, SizeType, DofHandler, MarsCrsMatrix::size_type>;

            using FE = Traits<FunctionSpace>::FE;
            using KokkosDiscretization = utopia::kokkos::Discretization<FunctionSpace, FE>;
            using Part = KokkosDiscretization::Part;

            static_assert(std::is_same<SizeType, Matrix::CrsMatrixType::global_ordinal_type>::value, "Weird!");

            auto new_crs_matrix() -> MarsCrsMatrix override { return sparsity_pattern->new_crs_matrix(); }

            void describe() const override {
                dof_handler->print_dofs();
                // sparsity_pattern->print_sparsity_pattern();

                // ::mars::print_fe_dof_map(*dof_handler_impl, *fe_dof_map_impl);
            }

            void matrix_apply_constraints(Matrix &m,
                                          const Scalar diag_value,
                                          const std::string side,
                                          const int component) override {
                // BC set constrained rows to zero, except diagonal where you set diag_value
                auto sp = *sparsity_pattern;
                auto mat = m.raw_type()->getLocalMatrixDevice();

                dof_handler->boundary_dof_iterate(
                    MARS_LAMBDA(const ::mars::Integer local_dof) {
                        sp.matrix_apply_constraints(local_dof, mat, diag_value);
                    },
                    side,
                    component);
            }

            void vector_apply_constraints(Vector &v,
                                          const Scalar value,
                                          const std::string side,
                                          const int component) override {
                // auto sp = *sparsity_pattern;
                auto vec = v.raw_type()->getLocalView<::mars::KokkosSpace>(Tpetra::Access::ReadWrite);
                auto sp_dof_handler = get_dof_handler();
                // BC set values to constraint value (i.e., boundary value)
                dof_handler->boundary_dof_iterate(
                    MARS_LAMBDA(const ::mars::Integer local_dof) {
                        sp_dof_handler.vector_apply_constraints(local_dof, vec, value);
                    },
                    side,
                    component);
            }

            void apply_zero_constraints(Vector &v, const std::string side, const int component) override {
                auto vec = v.raw_type()->getLocalView<::mars::KokkosSpace>(Tpetra::Access::ReadWrite);
                auto dof_handler = get_dof_handler();
                // BC set values to constraint value to zero
                dof_handler.boundary_dof_iterate(
                    MARS_LAMBDA(const ::mars::Integer local_dof) {
                        dof_handler.apply_zero_constraints(local_dof, vec);
                    },
                    side,
                    component);
            }

            void copy_at_constrained_nodes(const Vector &in,
                                           Vector &out,
                                           const std::string side,
                                           const int component) override {
                // auto sp = *sparsity_pattern;
                auto in_view = in.raw_type()->getLocalView<::mars::KokkosSpace>(Tpetra::Access::ReadOnly);
                auto out_view = out.raw_type()->getLocalView<::mars::KokkosSpace>(Tpetra::Access::ReadWrite);
                auto sp_dof_handler = get_dof_handler();

                // BC set values to constraint value (i.e., boundary value)
                dof_handler->boundary_dof_iterate(
                    MARS_LAMBDA(const ::mars::Integer local_dof) {
                        // sp.vector_apply_constraints(local_dof, vec, value);
                        auto idx = sp_dof_handler.local_to_owned_index(local_dof);
                        out_view(idx, 0) = in_view(idx, 0);
                    },
                    side,
                    component);
            }

            /* void system_apply_constraints(Matrix &m, Vector &v) override {
                vector_apply_constraints(v);
                matrix_apply_constraints(m, 1.0);
            } */

            SizeType n_local_dofs() override { return dof_handler->get_owned_dof_size(); };
            SizeType n_dofs() override { return dof_handler->get_global_dof_size(); };

            auto factory() -> Factory & override { return ConcreteFactory<DMesh>::instance(); };

            void init(DMesh &mesh_impl, int block_size) {
                dof_handler = std::make_shared<DofHandler>(mesh_impl);  //, mesh->raw_type_context());
                dof_handler->set_block(block_size);
                dof_handler->enumerate_dofs();

                fe_dof_map = std::make_shared<FEDofMap>(build_fe_dof_map(*dof_handler));

                // ensure_sparsity_pattern();
            }

            void ensure_sparsity_pattern() override {
                if (!sparsity_pattern) {
                    sparsity_pattern = std::make_shared<SPattern>(*dof_handler);
                    sparsity_pattern->build_pattern(*fe_dof_map);
                }
            }

            inline SPattern &get_sparsity_pattern() {
                if (!sparsity_pattern) {
                    Utopia::Abort("sparsity_pattern needs to be initialized with \"ensure_sparsity_pattern\"");
                }

                return *sparsity_pattern;
            }

            inline DofHandler &get_dof_handler() { return *dof_handler; }
            inline FEDofMap &get_fe_dof_map() { return *fe_dof_map; }

            inline const SPattern &get_sparsity_pattern() const { return *sparsity_pattern; }
            inline const DofHandler &get_dof_handler() const { return *dof_handler; }
            inline const FEDofMap &get_fe_dof_map() const { return *fe_dof_map; }

            void collect_ghost_layer(const Vector &in, ::mars::ViewVectorType<Scalar> &out) {
                UTOPIA_TRACE_REGION_BEGIN("mars::FEHandler::collect_ghost_layer");

                auto fe_dof_map = this->get_fe_dof_map();
                auto dof_handler = this->get_dof_handler();

                auto x_view = local_view_device(in).raw_type();  // Rank 2 tensor N x 1
                if (out.extent(0) != dof_handler.get_dof_size()) {
                    out = ::mars::ViewVectorType<Scalar>("x_local", dof_handler.get_dof_size());  // Rank 1 tensor
                }

                ::mars::ViewVectorType<Scalar> x_view_rank1(::Kokkos::subview(x_view, ::Kokkos::ALL, 0));
                ::mars::set_locally_owned_data(dof_handler, out, x_view_rank1);
                ::mars::gather_ghost_data(dof_handler, out);
                Kokkos::fence();

                UTOPIA_TRACE_REGION_END("mars::FEHandler::collect_ghost_layer");
            }

            ////////////////////////////////////////////////////////////////////////////////////

            void create(std::vector<KokkosDiscretization::FE_ptr> &fes,
                        int order,
                        const KokkosDiscretization::Part &part = KokkosDiscretization::all()) override {
                fes.clear();
                auto fe = std::make_shared<FE>();
                FEBuilder<FEHandler, FE> builder;
                builder.build(*this, *fe);
                fes.push_back(fe);

                std::cout << "n_cells: " << fe->n_cells() << "\n";

                assert(fe->n_cells() > 0);
            }

            void create_on_boundary(std::vector<KokkosDiscretization::FE_ptr> &fe,
                                    int order,
                                    const KokkosDiscretization::Part &part = KokkosDiscretization::all()) override {
                assert(false && "IMPLEMENT ME!");
                Utopia::Abort();
            }

            ////////////////////////////////////////////////////////////////////////////////////

            void convert_field(const Field<FunctionSpace> &in,
                               std::vector<std::shared_ptr<KokkosDiscretization::FEField>> &out,
                               const KokkosDiscretization::Part &part = KokkosDiscretization::all()) override {
                ::mars::ViewVectorType<Scalar> x;
                collect_ghost_layer(in.data(), x);

                int block_size = in.tensor_size();
                int b_offset_j = 0;

                auto fe_dof_map = this->get_fe_dof_map();
                auto dof_handler = this->get_dof_handler();

                assert(out.size() == 1);
                auto &&field = out[0]->data();
                auto fe = out[0]->fe();
                int n_fun = fe->n_shape_functions();

                size_t num_elem = fe_dof_map.get_fe_dof_map_size();

                if (Size_t(field.extent(0)) < num_elem || Size_t(field.extent(1)) < n_fun * block_size) {
                    field = KokkosDiscretization::FEField::DynRankView("Coefficients", num_elem, n_fun * block_size);
                }

                auto kernel = MARS_LAMBDA(const ::mars::Integer elem_index) {
                    Scalar val = 0;
                    for (int j = 0; j < n_fun; j++) {
                        for (int sub_j = 0; sub_j < block_size; ++sub_j) {
                            auto offset_j = dof_handler.compute_block_index(j, sub_j + b_offset_j);
                            const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, offset_j);

                            if (local_dof_j > -1) {
                                auto x_j = x(local_dof_j);
                                field(elem_index, j * block_size + sub_j) = x_j;
                            }
                        }
                    }
                };

                fe_dof_map.iterate(kernel);

                out[0]->set_tensor_size(in.tensor_size());
                out[0]->set_elem_type(in.elem_type());
            }

            void convert_field(const std::vector<std::shared_ptr<KokkosDiscretization::FEField>> &in,
                               Field<FunctionSpace> &out,
                               const KokkosDiscretization::Part &part = KokkosDiscretization::all()) override {
                assert(false && "IMPLEMENT ME!");
                Utopia::Abort();

                // assert(in.size() == 1);
                // auto &&field = in[0]->data();

                // fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                //     for (int i = 0; i < n_fun; i++) {
                //         for (int sub_i = 0; sub_i < block_size; ++sub_i) {
                //             auto offset_i = dof_handler.compute_block_index(i, sub_i + b_offset_i);
                //             const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, offset_i);
                //             if (local_dof_i < 0 || !dof_handler.is_owned(local_dof_i)) continue;

                //             Scalar val = field(elem_index, i * block_size + sub_i);

                //             const auto owned_dof_i = dof_handler.local_to_owned_index(local_dof_i);
                //             &y(owned_dof_i, 0) = val;
                //         }
                //     }
                // });
            }

            ////////////////////////////////////////////////////////////////////////////////////

            void global_to_local(const Vector &vector,
                                 std::vector<KokkosDiscretization::VectorAccumulator> &element_vectors,
                                 const KokkosDiscretization::Part &part = KokkosDiscretization::all(),
                                 const int comp = 0) override {
                assert(false && "IMPLEMENT ME!");
                Utopia::Abort();
            }

            // Local to global

            void local_to_global(const std::vector<KokkosDiscretization::MatrixAccumulator> &acc,
                                 AssemblyMode mode,
                                 Matrix &mat,
                                 const KokkosDiscretization::Part &part = KokkosDiscretization::all()) override {
                auto fe_dof_map = this->get_fe_dof_map();
                auto dof_handler = this->get_dof_handler();
                auto sp = this->get_sparsity_pattern();

                int b_offset_i = 0;
                int b_offset_j = 0;

                int block_size = dof_handler.get_block();
                assert(acc.size() == 1);

                auto &&m = acc[0];

                const int n_fun_i = m.extent(1) / block_size;
                const int n_fun_j = m.extent(2) / block_size;

                using SparseMatrixBuilder = ::mars::SparsityMatrix<SPattern>;

                SparseMatrixBuilder matrix_builder(sp, mat.raw_type()->getLocalMatrixDevice());

                assert(block_size <= dof_handler.get_block());

                // int rank = mat.comm().rank();

                // printf("%d) %d\n", rank, m.extent(0));

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    // printf("%d) idx: %d\n", rank, (int)elem_index);

                    // if (elem_index >= m.extent(0)) return;

                    assert(elem_index < m.extent(0));

                    for (int i = 0; i < n_fun_i; i++) {
                        for (int sub_i = 0; sub_i < block_size; ++sub_i) {
                            auto offset_i = dof_handler.compute_block_index(i, b_offset_i + sub_i);
                            const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, offset_i);
                            if (local_dof_i < 0 || !dof_handler.is_owned(local_dof_i)) continue;

                            for (int j = 0; j < n_fun_j; j++) {
                                for (int sub_j = 0; sub_j < block_size; ++sub_j) {
                                    auto offset_j = dof_handler.compute_block_index(j, b_offset_j + sub_j);
                                    const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, offset_j);
                                    if (local_dof_j > -1) {
                                        const Scalar val =
                                            m(elem_index, i * block_size + sub_i, j * block_size + sub_j);

                                        assert(val == val);
                                        matrix_builder.atomic_add_value(local_dof_i, local_dof_j, val);
                                    }
                                }
                            }
                        }
                    }
                });
            }

            void local_to_global(const std::vector<KokkosDiscretization::VectorAccumulator> &acc,
                                 AssemblyMode mode,
                                 Vector &vec,
                                 const KokkosDiscretization::Part &part = KokkosDiscretization::all()) override {
                auto fe_dof_map = this->get_fe_dof_map();
                auto dof_handler = this->get_dof_handler();

                int block_size = dof_handler.get_block();
                assert(acc.size() == 1);

                auto &&v = acc[0];

                const int n_fun = v.extent(1) / block_size;
                auto vec_view = local_view_device(vec).raw_type();

                auto kernel = MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < n_fun; i++) {
                        for (int sub_i = 0; sub_i < block_size; ++sub_i) {
                            auto offset_i = dof_handler.compute_block_index(i, sub_i);
                            const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, offset_i);
                            if (local_dof_i < 0 || !dof_handler.is_owned(local_dof_i)) continue;
                            const auto owned_dof_i = dof_handler.local_to_owned_index(local_dof_i);

                            auto val = v(elem_index, i * block_size + sub_i);
                            Kokkos::atomic_fetch_add(&vec_view(owned_dof_i, 0), val);
                        }
                    }
                };

                fe_dof_map.iterate(kernel);
            }

            void local_to_global(const Comm &comm,
                                 const std::vector<KokkosDiscretization::ScalarAccumulator> &acc,
                                 std::vector<Scalar> &scalars,
                                 const Part &part = KokkosDiscretization::all()) override {
                auto n = fe_dof_map->get_dof_handler().get_mesh().get_chunk_size();

                Scalar summed = 0;
                for (auto &a : acc) {
                    Scalar temp = 0;

                    Kokkos::parallel_reduce(
                        "scalar_local_to_global", n, UTOPIA_LAMBDA(const int i, Scalar &acc) { acc += a(i, 0); }, temp);

                    summed += temp;
                }

                summed = comm.sum(summed);
                scalars[0] = summed;
            }

            void local_to_global_on_boundary(
                const std::vector<KokkosDiscretization::VectorAccumulator> &acc,
                AssemblyMode mode,
                Vector &vec,
                const KokkosDiscretization::Part &part = KokkosDiscretization::all()) override {
                assert(false && "IMPLEMENT ME!");
                Utopia::Abort();
            }

        private:
            std::shared_ptr<SPattern> sparsity_pattern;
            std::shared_ptr<DofHandler> dof_handler;
            std::shared_ptr<FEDofMap> fe_dof_map;
        };

    }  // namespace mars
}  // namespace utopia

#endif
