// #ifndef UTOPIA_INTREPID2_L2_NORM_HPP
// #define UTOPIA_INTREPID2_L2_NORM_HPP

// #include "utopia_Input.hpp"
// #include "utopia_Traits.hpp"

// #include "utopia_intrepid2_Mass.hpp"

// namespace utopia {

//     namespace intrepid2 {
//         template <class FunctionSpace>
//         class L2Norm : public Configurable {
//         public:
//             using Scalar_t = typename Traits<FunctionSpace>::Scalar;
//             using Mass_t = intrepid2::Mass<Scalar_t>;
//             using MassAssembler_t = intrepid2::Assemble<Mass_t>;
//             using Intrepid2Field_t = intrepid2::Field<Scalar_t>;

//             void read(Input &) override {}

//             Scalar_t compute(const Field<FunctionSpace> &field) {
//                 if (!fe_) {
//                     fe_ = std::make_shared<intrepid2::FE<Scalar_t>>();
//                     create_fe(*field.space(), *fe_);
//                 }

//                 if (!mass_) {
//                     Mass<Mass_t> op;
//                     op.n_components = field.tensor_size();
//                     mass_ = std::make_shared<MassAssembler_t>(fe_, op);
//                     mass_->ensure_vector_accumulator();
//                     mass_->set_mode(OVERWRITE_MODE);
//                 }

//                 if (!interpid2_field_) {
//                     interpid2_field_ = std::make_shared<Intrepid2Field_t>();
//                 }

//                 convert_field(field, *interpid2_field_);

//                 mass_->update(interpid2_field_);
//                 mass_->assemble_vector();

//                 auto local_result = interpid2_field_->dot(mass_->vector_accumulator()->data());

//                 auto global_result = field.data().comm().sum(local_result);
//                 return std::sqrt(global_result);
//             }

//         private:
//             std::shared_ptr<intrepid2::FE<Scalar_t>> fe_;
//             std::shared_ptr<MassAssembler_t> mass_;
//             std::shared_ptr<Intrepid2Field_t> interpid2_field_;
//         };

//         template <class FunctionSpace>
//         auto l2_norm(const Field<FunctionSpace> &field) -> typename Traits<FunctionSpace>::Scalar {
//             return L2Norm<FunctionSpace>().compute(field);
//         }

//     }  // namespace intrepid2
// }  // namespace utopia

// #endif  // UTOPIA_INTREPID2_L2_NORM_HPP
