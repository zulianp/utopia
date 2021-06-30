// #ifndef UTOPIA_IMPLICIT_OBSTACLE_FE_FUNCTION_HPP
// #define UTOPIA_IMPLICIT_OBSTACLE_FE_FUNCTION_HPP

// #include "utopia_BoxConstrainedFEFunction.hpp"
// #include "utopia_fe_Core.hpp"

// namespace utopia {

//     template <class FunctionSpace>
//     class ImplicitObstacleFEFunction final : public BoxConstrainedFEFunction<FunctionSpace> {
//     public:
//         using Super = utopia::BoxConstrainedFEFunction<FunctionSpace>;
//         using Vector_t = typename Traits<FunctionSpace>::Vector;
//         using Matrix_t = typename Traits<FunctionSpace>::Matrix;

//         using Size_t = typename Traits<FunctionSpace>::SizeType;
//         using Scalar_t = typename Traits<FunctionSpace>::Scalar;
//         using Communicator_t = typename Traits<FunctionSpace>::Communicator;
//         using Mesh_t = typename Traits<FunctionSpace>::Mesh;

//         // Use specialized components for function space
//         using ImplicitObstacle_t = utopia::ImplicitObstacle<FunctionSpace>;

//         bool has_nonlinear_constraints() const override { return !linear_obstacle_; }

//         bool update_constraints(const Vector_t &x) {
//             utopia::out() << "update_constraints\n";

//             this->space()->displace(x);
//             bool ok = obstacle_->assemble(*this->space());

//             if (debug_) {
//                 static int iter_debug = 0;
//                 ouput_debug_data(iter_debug++, x);
//             }

//             this->space()->displace(-x);
//             return ok;
//         }

//         bool has_transformation() const override { return true; }

//         void transform(const Vector_t &x, Vector_t &x_constrained) override { obstacle_->transform(x, x_constrained);
//         }

//         void inverse_transform(const Vector_t &x_constrained, Vector_t &x) override {
//             obstacle_->inverse_transform(x_constrained, x);
//         }

//         void transform(const Matrix_t &H, Matrix_t &H_constrained) override { obstacle_->transform(H, H_constrained);
//         }

//         bool constraints_gradient(const Vector_t &x, BoxConstraints<Vector_t> &box) override {
//             update_constraints(x);

//             if (!box.upper_bound()) {
//                 box.upper_bound() = std::make_shared<Vector_t>();
//             }

//             (*box.upper_bound()) = obstacle_->gap();
//             return true;
//         }

//         ImplicitObstacleFEFunction(const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &unconstrained)
//             : Super(unconstrained) {}

//         void read(Input &in) override {
//             Super::read(in);

//             in.get("linear_obstacle", linear_obstacle_);
//             in.get("debug", debug_);

//             if (!obstacle_) {
//                 obstacle_ = std::make_shared<ImplicitObstacle_t>();
//                 in.require("obstacle", obstacle_);
//             }
//         }

//     private:
//         std::shared_ptr<ImplicitObstacle_t> obstacle_;
//         FETransfer<FunctionSpace> transfer_;
//         bool linear_obstacle_{false};
//         bool debug_{false};
//     };

// }  // namespace utopia

// #endif  // UTOPIA_IMPLICIT_OBSTACLE_FE_FUNCTION_HPP
