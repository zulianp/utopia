// #ifndef UTOPIA_MONODOMAIN_FE_HPP
// #define UTOPIA_MONODOMAIN_FE_HPP


// #include "utopia_Function.hpp"
// #include "utopia_LaplacianView.hpp"
// #include "utopia_MassMatrixView.hpp"

// #include "utopia_QuadratureView.hpp"
// #include "utopia_MPITimeStatistics.hpp"
// #include "utopia_NodalInterpolateView.hpp"
// #include "utopia_Algorithms.hpp"
// #include "utopia_DeviceReduce.hpp"
// #include "utopia_DeviceTensorReduce.hpp"

// namespace utopia {

//     template<class FunctionSpace>
//     class MonodomainFE final : public Function<
//                             typename FunctionSpace::Matrix,
//                             typename FunctionSpace::Vector> {
//     public:
//         using Matrix     = typename FunctionSpace::Matrix;
//         using Vector     = typename FunctionSpace::Vector;
//         using Scalar     = typename Traits<Vector>::Scalar;
//         using SizeType   = typename Traits<Vector>::SizeType;
//         using Elem       = typename FunctionSpace::Elem;
//         using Quadrature = utopia::Quadrature<Elem, 2>;
//         using Laplacian  = utopia::Laplacian<FunctionSpace, Quadrature>;
//         using ScaledMassMatrix = utopia::ScaledMassMatrix<FunctionSpace, Quadrature>;

//         using Device     = typename FunctionSpace::Device;

//         static const int Dim    = Elem::Dim;
//         static const int NNodes = Elem::NNodes;
//         using ElementMatrix     = utopia::StaticMatrix<Scalar, NNodes, NNodes>;
//         using ElementVector     = utopia::StaticVector<Scalar, NNodes>;
//         using Point             = typename FunctionSpace::Point;


//         class FitzHughNagumoReactionTerm {
//         public:
//             UTOPIA_INLINE_FUNCTION FitzHughNagumoReactionTerm() {}

//             UTOPIA_INLINE_FUNCTION Scalar eval(const Scalar &u) const
//             {
//                 return alpha * (u - u_min) * (u - u_max) * (u - u_unst);
//             }

//             UTOPIA_INLINE_FUNCTION Scalar eval_derivative(const Scalar &u) const
//             {
//                 return alpha * ((u - u_min * u - u_max) + (u - u_min * u - u_unst) + (u - u_unst * u - u_max));
//             }

//             Scalar u_min = -85., u_max = 30., u_unst = -57.6, u_source = -0.5;
//             Scalar alpha = -1.4e-3;
//         };

//         inline bool value(const Vector &x, Scalar &value) const override
//         {
//             // auto space_view = space_->view_device();

//             // auto x_view = space_->assembly_view_device(x);

//             // auto l_view = laplacian_.view_device();
//             // auto m_view = scaled_mass_matrix_.view_device();

//             // auto reaction_term = UTOPIA_LAMBDA(const Scalar &u) -> Scalar {
//             //     return -lambda_ * device::exp(u);
//             // };

//             // value = 0.0;

//             // Device::parallel_reduce(
//             //     space_->element_range(),
//             //     UTOPIA_LAMBDA(const SizeType &i) -> Scalar
//             // {
//             //     Elem e;
//             //     space_view.elem(i, e);

//             //     ElementVector coeff;

//             //     space_view.coefficients(e, x_view, coeff);

//             //     ElementMatrix el_mat;
//             //     el_mat.set(0.0);

//             //     l_view.add(i, e, el_mat);

//             //     ElementVector el_vec;
//             //     el_vec = el_mat * coeff;

//             //     Scalar v = dot(el_vec, coeff);

//             //     el_mat.set(0.0);
//             //     m_view.add(i, e, reaction_term, el_mat);

//             //     v += sum(el_mat);
//             //     return v;
//             // }, value);

//             // value = x.comm().sum(value);
//             return false;
//         }

//         inline bool update(const Vector &point) override {
//             scaled_mass_matrix_.update(point);
//             return true;
//         }

//         inline bool gradient(const Vector &x, Vector &g) const override
//         {
//             if(empty(g)) {
//                 space_->create_vector(g);
//             } else {
//                 g.set(0.0);
//             }

//             {
//                 auto space_view = space_->view_device();

//                 auto x_view = space_->assembly_view_device(x);
//                 auto g_view = space_->assembly_view_device(g);

//                 auto l_view = laplacian_.view_device();
//                 auto m_view = scaled_mass_matrix_.view_device();

//                 auto reaction_term = UTOPIA_LAMBDA(const Scalar &u) -> Scalar {
//                     return -fhn_.eval(u);
//                 };

//                 Device::parallel_for(
//                     space_->element_range(),
//                     UTOPIA_LAMBDA(const SizeType &i)
//                 {
//                     Elem e;
//                     space_view.elem(i, e);

//                     ElementVector coeff;

//                     space_view.coefficients(e, x_view, coeff);

//                     ElementMatrix el_mat;
//                     el_mat.set(0.0);

//                     l_view.add(i, e, el_mat);

//                     ElementVector el_vec;
//                     el_vec = el_mat * coeff;

//                     el_mat.set(0.0);
//                     m_view.add(i, e, reaction_term, el_mat);

//                     el_vec += row_sum(el_mat);

//                     space_view.add_vector(e, el_vec, g_view);
//                 });
//             }

//             space_->apply_constraints(g);
//             return true;
//         }

//         inline bool hessian(const Vector &, Matrix &H) const override
//         {
//             if(empty(H)) {
//                 space_->create_matrix(H);
//             } else {
//                 H *= 0.0;
//             }

//             {
//                 auto space_view = space_->view_device();

//                 auto H_view = space_->assembly_view_device(H);
//                 auto l_view = laplacian_.view_device();
//                 auto m_view = scaled_mass_matrix_.view_device();

//                 auto f = UTOPIA_LAMBDA(const Scalar &u) -> Scalar {
//                     return -fhn_.eval_derivative(u);
//                 };

//                 Device::parallel_for(
//                     space_->element_range(),
//                     UTOPIA_LAMBDA(const SizeType &i)
//                 {
//                     Elem e;
//                     space_view.elem(i, e);

//                     ElementMatrix el_mat;
//                     el_mat.set(0.0);

//                     l_view.add(i, e, el_mat);
//                     m_view.add(i, e, f, el_mat);

//                     space_view.add_matrix(e, el_mat, H_view);
//                 });
//             }

//             space_->apply_constraints(H);
//             return true;
//         }

//         inline bool has_preconditioner() const override
//         {
//             return false;
//         }

//         inline bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const
//         {
//             space_->create_matrix(H);
//             return true;
//         }

//         MonodomainFE(FunctionSpace &space)
//         : space_(utopia::make_ref(space)),
//           quadrature_(),
//           laplacian_(space, quadrature_),
//           scaled_mass_matrix_(space, quadrature_)
//         {}

//         MonodomainFE(const std::shared_ptr<FunctionSpace> &space)
//         : space_(space),
//           quadrature_(),
//           laplacian_(*space, quadrature_),
//           scaled_mass_matrix_(*space, quadrature_)
//         {}

//     private:
//         std::shared_ptr<FunctionSpace> space_;
//         Quadrature quadrature_;
//         Laplacian laplacian_;
//         ScaledMassMatrix scaled_mass_matrix_;

//         Scalar chi_   = 1400.;
//         Scalar diffusivity_ = 1.33;
//         Scalar dt_ = 0.1;

//         FitzHughNagumoReactionTerm fhn_;
//     };
// }


// #endif //UTOPIA_MONODOMAIN_FE_HPP
