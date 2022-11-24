#include "utopia_mars_FEAssembler_new.hpp"

//     bool update_input(const Vector &x) {
//         auto space = this->discretization()->space();
//         utopia::Field<FunctionSpace> in("x", space, make_ref(const_cast<Vector &>(x)));

//         // FIXME does not work for mixed FE
//         in.set_tensor_size(space->n_var());

//         if (!current_solution_) {
//             assert(fe_);
//             current_solution_ = std::make_shared<Field<FE>>(fe_);
//         }

//         convert_field(in, *current_solution_);
//         return true;
//     }

//     void matrix_assembly_begin(Matrix &, AssemblyMode) {}

//     void matrix_assembly_end(Matrix &matrix, AssemblyMode mode) {
//         this->discretization()->local_to_global({this->matrix_data()}, mode, matrix);
//     }

//     void vector_assembly_begin(Vector &, AssemblyMode) {}

//     void vector_assembly_end(Vector &vector, AssemblyMode mode) {
//         // this->discretization()->local_to_global(this->vector_data(), mode, vector);
//     }

//     void scalar_assembly_begin(Scalar &scalar, AssemblyMode mode) {
//         assert(false);
//         Utopia::Abort();
//     }

//     void scalar_assembly_end(Scalar &scalar, AssemblyMode mode) {
//         assert(false);
//         Utopia::Abort();
//     }