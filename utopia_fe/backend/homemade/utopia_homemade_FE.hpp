// #ifndef UTOPIA_HOMEMADE_FE_HPP
// #define UTOPIA_HOMEMADE_FE_HPP

// #include "utopia_Base.hpp"
// #include "utopia_homemade_FEForwardDeclarations.hpp"

// #include "utopia_homemade_FETypes.hpp"
// #include "utopia_intersector.hpp"

// namespace utopia {

//     class FE {
//     public:
//         class Impl;

//         FE();
//         ~FE();
//         void init(const int current_element, Mesh &mesh, int quadrature_order);
//         int n_shape_functions() const;

//         HMDerivative grad;
//         HMFun fun;
//         HMDx dx;

//     private:
//         std::unique_ptr<Impl> impl_ptr;
//     };

// }

// #endif //UTOPIA_HOMEMADE_FE_HPP