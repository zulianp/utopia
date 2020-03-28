// #ifndef UTOPIA_L2_NORM_HPP
// #define UTOPIA_L2_NORM_HPP

// #include "utopia_AssemblyView.hpp"

// namespace utopia {
//     template<class FunctionSpace, class Quadrature>
//     class L2Norm {
//     public:
//         using Scalar        = typename FunctionSpace::Scalar;
//         using SizeType      = typename FunctionSpace::SizeType;
//         using Device        = typename FunctionSpace::Device;
//         using Elem          = typename FunctionSpace::ViewDevice::Elem;
//         using Differential  = utopia::Differential<FunctionSpace, Quadrature>;
//         using PhysicalPoint = utopia::PhysicalPoint<FunctionSpace, Quadrature>;


//         L2Norm(const FunctionSpace &space, const Quadrature &q)
//         : space_(space), q_(q), differential_(space, q), point_(space, q)
//         {}

//         template<class F>
//         Scalar apply(F f) const
//         {

//             Scalar norm_f = 0.0;

//             {
//                 auto dx_view = differential_.view_device();
//                 auto p_view  = point_.view_device();
//                 auto space_view = space_.view_device();

//                 Device::parallel_reduce(
//                     space_->element_range(),
//                     UTOPIA_LAMBDA(const SizeType &i) -> Scalar
//                 {
//                     Elem e;
//                     space_view.elem(i, e);
//                     Scalar el_norm_f = 0.0;

//                     auto dx = dx_view.make(e);
//                     auto p  = p_view.make(e);

//                     for(SizeType qp = 0; qp < Quadrature::NPoints; ++qp) {
//                         el_norm_f += f(p(qp)) * dx(qp);
//                     }

//                     return el_norm_f;
//                 }, norm_f);
//             }

//             return std::sqrt(space_.comm().sum(norm_f));
//         }

//     private:
//         const FunctionSpace &space_;
//         const Quadrature &q_;
//         Differential differential_;
//         PhysicalPoint point_;
//     };

//     template<class F, class FunctionSpace, class Quadrature>
//     typename FunctionSpace::Scalar l2_norm(
//         F f,
//         const FunctionSpace &space,
//         const Quadrature &q)
//     {
//         return L2Norm<FunctionSpace, Quadrature>(space).apply(f, q);
//     }
// }

// #endif //UTOPIA_L2_NORM_HPP
