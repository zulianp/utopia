#ifndef UTOPIA_KOKKOS_BasisFunction_hpp
#define UTOPIA_KOKKOS_BasisFunction_hpp

#include "utopia_Base.hpp"

namespace utopia {

  template<typename T, int N_>
  class Array {
  public:
    constexpr static const int N = N_; 
    UTOPIA_INLINE_FUNCTION constexpr Array()
    {}

    UTOPIA_INLINE_FUNCTION constexpr int size() const
    {
        return N;
    }

    T values[N];
  };

  template<typename T>
  class Quad4 {

    using Scalar = T;
  
  public:

    UTOPIA_INLINE_FUNCTION constexpr Quad4()
    {}


    UTOPIA_INLINE_FUNCTION constexpr int size() const
    {
        return 4;
    }

   /*evaluation basis function*/
    UTOPIA_INLINE_FUNCTION static void eval_phi(const Scalar *p, Array<Scalar, 4> &fn) {

        const Scalar u = 1.0 - p[0];
        const Scalar v = 1.0 - p[1];

        const Scalar up = p[0];
        const Scalar vp = p[1];


        fn.values[0] = u * v;        
        fn.values[1] = up * v;        
        fn.values[2] = vp * up;       
        fn.values[3] = u * vp;
    }


    /*evaluation gradient of basis function*/
    UTOPIA_INLINE_FUNCTION static void eval_grad(const Scalar *p, Array<Scalar, 8> &grad_fn) {


        const Scalar u = 1.0 - p[0];
        const Scalar v = 1.0 - p[1];


        const Scalar up = p[0];
        const Scalar vp = p[1];


      //gfn 0
        grad_fn.values[0] = -1.0 *  v;
        grad_fn.values[1] = -1.0 *  u;
      
      //gfn 1
        grad_fn.values[2] = v;
        grad_fn.values[3] = -1.0 * vp;
      
      //gfn 2
        grad_fn.values[4] = vp;
        grad_fn.values[5] = up;
     
      //gfn 3
        grad_fn.values[6] =  -1.0 *  vp;
        grad_fn.values[7] =  u;
        
    }

     /*evaluation jacobian transformation*/
    UTOPIA_INLINE_FUNCTION static void eval_J(const Scalar &h, Array<Scalar, 4> &J) {



      //gfn 0
        J.values[0] = h;
        J.values[1] = 0.0;
        J.values[2] = 0.0;
        J.values[3] = h;

        
    }


         /*evaluation jacobian transformation*/
    UTOPIA_INLINE_FUNCTION static void eval_inv_J(const Scalar &h, Array<Scalar, 4> &J) {



      //gfn 0
        J.values[0] = 1./h;
        J.values[1] = 0.0;
        J.values[2] = 0.0;
        J.values[3] = 1./h;

        
    }


         /*evaluation jacobian transformation*/
    UTOPIA_INLINE_FUNCTION static void eval_det_J(const Scalar &h, Scalar &J) {

        J= h*h;

        
    }

   UTOPIA_INLINE_FUNCTION Scalar static measure ()
   {
     return 1.0;
   }

private:

  };


  template<typename T>
  class Hex8 {

    using Scalar = T;
  
  public:

    UTOPIA_INLINE_FUNCTION constexpr Hex8()
    {}


    UTOPIA_INLINE_FUNCTION constexpr int size() const
    {
        return 8;
    }

   /*evaluation basis function*/
    UTOPIA_INLINE_FUNCTION static void eval_phi(const Scalar *p, Array<Scalar, 8> &fn) {

        
        const Scalar u = 1.0 - p[0];
        const Scalar v = 1.0 - p[1];
        const Scalar w = 1.0 - p[2];

        const Scalar up = p[0];
        const Scalar vp = p[1];
        const Scalar wp = p[2];

        fn.values[0] = u * v * w;
        
        fn.values[1] = up * v * w;
        
        fn.values[2] = vp * up * w;
        
        fn.values[3] = u * vp * w;
        
        fn.values[4] = u * v * wp;
        
        fn.values[5] = up * v * wp;
        
        fn.values[6] = up * vp * wp;
        
        fn.values[7] = u * vp * wp;
    }


    /*evaluation gradient of basis function*/
    UTOPIA_INLINE_FUNCTION static void eval_grad(const Scalar *p, Array<Scalar, 23> &grad_fn) {


        const Scalar u = 1.0 - p[0];
        const Scalar v = 1.0 - p[1];
        const Scalar w = 1.0 - p[2];

        const Scalar up = p[0];
        const Scalar vp = p[1];
        const Scalar wp = p[2];

      //gfn 0
        grad_fn.values[0] = -1.0 *  v  *  w;
        grad_fn.values[1] = -1.0 *  u  *  w;
        grad_fn.values[2] = -1.0 *  u  *  v;
      
      //gfn 1
        grad_fn.values[3] =  1.0 *  v   *  w;
        grad_fn.values[4] = -1.0 *  up  *  w;
        grad_fn.values[5] = -1.0 *  up  *  v;
      
      //gfn 2
        grad_fn.values[6] =  1.0 *  vp  *  w;
        grad_fn.values[7] =  1.0 *  up  *  w;
        grad_fn.values[8] = -1.0 *  up  *  vp;
      
      //gfn 3
        grad_fn.values[9]  = -1.0 *  vp *  w;
        grad_fn.values[10] =  1.0 *  u  *  w;
        grad_fn.values[11] = -1.0 *  u  *  vp;
      
      //gfn 4
        grad_fn.values[12] = -1.0 *  v  *  wp;
        grad_fn.values[13] = -1.0 *  u  *  wp;
        grad_fn.values[14] =  1.0 *  u  *  v;
      
      //gfn 5
        grad_fn.values[15] =  1.0 *  v   *  wp;
        grad_fn.values[16] = -1.0 *  up  *  wp;
        grad_fn.values[17] =  1.0 *  up  *  v;
      
      //gfn 6
        grad_fn.values[18] =  1.0 *  vp  *  wp;
        grad_fn.values[19] =  1.0 *  up  *  wp;
        grad_fn.values[20] =  1.0 *  up  *  vp;
      
      //gfn 7
        grad_fn.values[21] = -1.0 *  vp  *  wp;
        grad_fn.values[22] =  1.0 *  u   *  wp;
        grad_fn.values[23] =  1.0 *  u   *  vp;
        
    }


         /*evaluation jacobian transformation*/
    UTOPIA_INLINE_FUNCTION static void eval_J(const Scalar &h, Array<Scalar, 9> &J) {



      //gfn 0
        J.values[0] = h;
        J.values[1] = 0.0;
        J.values[2] = 0.0;
        J.values[3] = 0.0;
        J.values[4] = h;
        J.values[5] = 0.0;
        J.values[6] = 0.0;
        J.values[7] = 0.0;
        J.values[8] = h;


        
    }


         /*evaluation jacobian transformation*/
    UTOPIA_INLINE_FUNCTION static void eval_inv_J(const Scalar &h, Array<Scalar, 9> &J) {



        J.values[0] = 1./h;
        J.values[1] = 0.0;
        J.values[2] = 0.0;
        J.values[3] = 0.0;
        J.values[4] = 1./h;
        J.values[5] = 0.0;
        J.values[6] = 0.0;
        J.values[7] = 0.0;
        J.values[8] = 1./h;

        
    }


         /*evaluation jacobian transformation*/
    UTOPIA_INLINE_FUNCTION static void eval_det_J(const Scalar &h, Scalar &J) {



      //gfn 0
        J= h*h*h;
        
    }

   UTOPIA_INLINE_FUNCTION Scalar static measure ()
   {
     return 1.0;
   }

   private:

  };


  // template<class Q, class FE>
  // class AssembleMassMatrix {
  // public:
  //       Simplex<4, 3 /*KokkosImplementation*/> elem;
  //       ViewVectorType<bool> active;
  //       Kokkos::View<Real> detJ;
  //       ViewMatrixType<Integer> element_matrix;

  //       AssembleMassMatrix(
  //           const Simplex<4, 3 /*KokkosImplementation*/> &el,
  //           ViewVectorType<bool> ac,
  //           Kokkos::View<Real> J,
  //           ViewMatrixType<Integer> element_matrix) 
  //       : elem(el), active(ac), detJ(J), element_matrix(element_matrix) {}

  //       AssembleMassMatrix() {}




  //         UTOPIA_INLINE_FUNCTION void operator()(int index) const
  //         {
             

  //               constexpr const Q q;

  //               constexpr const FE phi;

  //               //element_matrix(index) = compute_ele_mat_matrix(J, el);
              
  //         }
 
  //   };




  

}
#endif


