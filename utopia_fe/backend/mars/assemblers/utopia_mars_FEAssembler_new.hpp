#ifndef UTOPIA_MARS_FEASSEMBLER_NEW_HPP
#define UTOPIA_MARS_FEASSEMBLER_NEW_HPP

#include "utopia_Instance.hpp"
#include "utopia_Traits.hpp"

#include "utopia_FEAssembler.hpp"

#include "utopia_SimulationTime.hpp"
#include "utopia_fe_Core.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_Field.hpp"

#include "utopia_kokkos_TensorAccumulator.hpp"

#include "utopia_kokkos_SubdomainValue.hpp"

#include "utopia_mars_Discretization.hpp"

#include <memory>

namespace utopia {

    namespace mars {

        template <class FE_>
        class FEAssemblerNew : public Describable {
        public:
            using Matrix = Traits<FunctionSpace>::Matrix;
            using Vector = Traits<FunctionSpace>::Vector;
            using Scalar = Traits<FunctionSpace>::Scalar;
            using FE = FE_;
            using Discretization = utopia::Discretization<utopia::mars::FunctionSpace, FE>;
            using Part = typename Discretization::Part;

            inline std::shared_ptr<kokkos::Field<FE>> current_solution() { return current_solution_; }
            std::shared_ptr<Discretization> discretization() { return nullptr; }

            template <class Op>
            bool assemble_matrix_eij(const std::string &name,
                                     AssemblyMode mode,
                                     Op op,
                                     const Part part = Discretization::all()) {
                //         ensure_matrix_accumulator();
                //         auto data = this->matrix_data();

                //         const Scalar a = OVERWRITE_MODE == mode ? 0 : 1;
                //         const Scalar b = SUBTRACT_MODE == mode ? -1 : 1;

                //         this->loop_cell_test_trial(
                //             name, UTOPIA_LAMBDA(const int cell, const int i, const int j) {
                //                 auto value = op(cell, i, j);
                //                 auto &v = data(cell, i + op.offset_test(), j + op.offset_trial());
                //                 v += a * v + b * value;
                //             });

                //         // matrix_accumulator_->describe(std::cout);

                return true;
            }

            //     // Block
            //     template <class Op>
            //     bool assemble_matrix_eij_block(const std::string &name,
            //                                    AssemblyMode mode,
            //                                    Op op,
            //                                    const Part part = Discretization::all()) {
            //         ensure_matrix_accumulator();
            //         auto data = this->matrix_data();

            //         const Scalar a = OVERWRITE_MODE == mode ? 0 : 1;
            //         const Scalar b = SUBTRACT_MODE == mode ? -1 : 1;

            //         this->loop_cell_test_trial(
            //             name, UTOPIA_LAMBDA(const int cell, const int i, const int j) {
            //                 StaticMatrix<Scalar, Op::NComponentsTest, Op::NComponentsTrial> block;
            //                 block.set(0.);

            //                 // Evaluate block operator
            //                 op(cell, i, j, block);

            //                 for (int di = 0; di < Op::NComponentsTest; ++di) {
            //                     for (int dj = 0; dj < Op::NComponentsTrial; ++dj) {
            //                         auto &v = data(cell, op.offset_test() + di, op.offset_trial() + dj);
            //                         v += a * v + b * block(di, dj);
            //                     }
            //                 }
            //             });

            //         return true;
            //     }

            template <class Op>
            bool assemble_apply_ei(const std::string &name,
                                   AssemblyMode mode,
                                   Op op,
                                   const kokkos::Field<FE> &field,
                                   const Part part = Discretization::all()) {
                //         ensure_vector_accumulator();

                //         auto data = this->vector_data();
                //         auto x = field.data();

                //         const Scalar a = OVERWRITE_MODE == mode ? 0 : 1;
                //         const Scalar b = SUBTRACT_MODE == mode ? -1 : 1;

                //         auto &fe = this->fe();
                //         const int n_shape_functions = fe.n_shape_functions();

                //         this->loop_cell_test(
                //             name, UTOPIA_LAMBDA(const int cell, const int i) {
                //                 Scalar val = 0.0;
                //                 for (int j = 0; j < n_shape_functions; ++j) {
                //                     val += op(cell, i, j) * x(cell, j);
                //                 }
                //                 auto &v = data(cell, i);
                //                 v = a * v + b * val;
                //             });

                return true;
            }

            //     template <class Op>
            //     bool assemble_vector_ei_block(const std::string &name,
            //                                   AssemblyMode mode,
            //                                   Op op,
            //                                   const Part part = Discretization::all()) {
            //         ensure_vector_accumulator();
            //         auto data = this->vector_data();

            //         const Scalar a = OVERWRITE_MODE == mode ? 0 : 1;
            //         const Scalar b = SUBTRACT_MODE == mode ? -1 : 1;

            //         this->loop_cell_test(
            //             name, UTOPIA_LAMBDA(const int cell, const int i) {
            //                 StaticVector<Scalar, Op::NComponentsTest> block;
            //                 block.set(0.);

            //                 // Evaluate block operator
            //                 op(cell, i, block);

            //                 for (int di = 0; di < Op::NComponentsTest; ++di) {
            //                     auto &v = data(cell, op.offset_test() + di);
            //                     v = a * v + b * block(di);
            //                 }
            //             });

            //         return true;
            //     }

            //     template <class Op>
            //     bool assemble_vector_ei(const std::string &name,
            //                             AssemblyMode mode,
            //                             Op op,
            //                             const Part part = Discretization::all()) {
            //         ensure_vector_accumulator();

            //         auto data = this->vector_data();

            //         const Scalar a = OVERWRITE_MODE == mode ? 0 : 1;
            //         const Scalar b = SUBTRACT_MODE == mode ? -1 : 1;

            //         this->loop_cell_test(
            //             name, UTOPIA_LAMBDA(const int cell, const int i) {
            //                 const Scalar val = op(cell, i);
            //                 auto &v = data(cell, i);
            //                 v = a * v + b * val;
            //             });

            //         return true;
            //     }

            //     template <class Op>
            //     bool assemble_scalar_e(const std::string &name,
            //                            AssemblyMode mode,
            //                            Op op,
            //                            const Part part = Discretization::all()) {
            //         ensure_scalar_accumulator();

            //         auto data = this->scalar_data();

            //         const Scalar a = OVERWRITE_MODE == mode ? 0 : 1;
            //         const Scalar b = SUBTRACT_MODE == mode ? -1 : 1;

            //         this->loop_cell_test(
            //             name, UTOPIA_LAMBDA(const int &cell) {
            //                 const Scalar val = op(cell);
            //                 auto &v = data(cell);
            //                 v = a * v + b * val;
            //             });

            //         return true;
            //     }

            //     template <class Op>
            //     bool assemble_scalar_e_block(const std::string &name,
            //                                  AssemblyMode mode,
            //                                  Op op,
            //                                  const Part part = Discretization::all()) {
            //         ensure_scalar_accumulator();

            //         auto data = this->scalar_data();

            //         const Scalar a = OVERWRITE_MODE == mode ? 0 : 1;
            //         const Scalar b = SUBTRACT_MODE == mode ? -1 : 1;

            //         this->loop_cell_test(
            //             name, UTOPIA_LAMBDA(const int cell) {
            //                 StaticVector<Scalar, Op::NComponentsTest> block;
            //                 block.set(0.);

            //                 op(cell, block);

            //                 for (int di = 0; di < Op::NComponentsTest; ++di) {
            //                     auto &v = data(cell, op.offset() + di);
            //                     v = a * v + b * block(di);
            //                 }
            //             });

            //         return true;
            //     }

            std::shared_ptr<FE> fe_ptr() { return fe_; }
            FE &fe() { return *fe_; }
            const FE &fe() const { return *fe_; }

            //     //////////////////////////////////////////////////////////////////////////////////////////////////

            bool update_input(const Vector &x) {
                // auto space = this->discretization()->space();
                // utopia::Field<FunctionSpace> in("x", space, make_ref(const_cast<Vector &>(x)));

                // // FIXME does not work for mixed FE
                // in.set_tensor_size(space->n_var());

                // if (!current_solution_) {
                //     assert(fe_);
                //     current_solution_ = std::make_shared<Field<FE>>(fe_);
                // }

                // convert_field(in, *current_solution_);
                // return true;

                assert(false);
                Utopia::Abort();
                return false;
            }

            void matrix_assembly_begin(Matrix &, AssemblyMode) {
                assert(false);
                Utopia::Abort();
            }

            void matrix_assembly_end(Matrix &matrix, AssemblyMode mode) {
                // this->discretization()->local_to_global({this->matrix_data()}, mode, matrix);
                assert(false);
                Utopia::Abort();
            }

            void vector_assembly_begin(Vector &, AssemblyMode) {}

            void vector_assembly_end(Vector &vector, AssemblyMode mode) {
                // this->discretization()->local_to_global(this->vector_data(), mode, vector);
                assert(false);
                Utopia::Abort();
            }

            void scalar_assembly_begin(Scalar &scalar, AssemblyMode mode) {
                assert(false);
                Utopia::Abort();
            }

            void scalar_assembly_end(Scalar &scalar, AssemblyMode mode) {
                assert(false);
                Utopia::Abort();
            }

        private:
            std::shared_ptr<FE> fe_;
            std::shared_ptr<kokkos::Field<FE>> current_solution_;
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_FEASSEMBLER_NEW_HPP
