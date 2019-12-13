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
    using Point= Array<T, 2>;
  
  public:

    UTOPIA_INLINE_FUNCTION constexpr Quad4()
    {}


    UTOPIA_INLINE_FUNCTION constexpr int size() const
    {
        return 4;
    }


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
    UTOPIA_INLINE_FUNCTION void eval_J(Array<Scalar, 4> &J) {


          for(int i = 0; i < 2; ++i) {
               const int i_offset = i * 2;
                J[i_offset]     = node(1)[i] - node(0)[i];
                J[i_offset + 1] = node(3)[i] - node(0)[i];
           }

        
    }

         /*evaluation jacobian transformation*/
    UTOPIA_INLINE_FUNCTION void eval_det_J(Scalar& detJ) {

        Array<Scalar, 4> J;

          for(int i = 0; i < 2; ++i) {
               const int i_offset = i * 2;
                J[i_offset]     = node(1)[i] - node(0)[i];
                J[i_offset + 1] = node(3)[i] - node(0)[i];
           }

        detJ = J[0]*J[3]-J[2]*J[1];
    }

    UTOPIA_INLINE_FUNCTION void eval_inv_J(Array<Scalar, 4> invJ) {
        
        Array<Scalar, 4> J;

          for(int i = 0; i < 2; ++i) {
               const int i_offset = i * 2;
                J[i_offset]     = node(1)[i] - node(0)[i];
                J[i_offset + 1] = node(3)[i] - node(0)[i];
           }

          Scalar detJ = J[0]*J[3]-J[2]*J[1];

          Scalar inv_detJ = 1.0/detJ;

            invJ[0] =  J[3]*inv_detJ;
            invJ[1] = -J[1]*inv_detJ;
            invJ[2] = -J[2]*inv_detJ;
            invJ[3] =  J[1]*inv_detJ;

    }

    Point &node(const int i) override
    {
        return nodes_[i];
    }

    const Point &node(const int i) const override
    {
        return nodes_[i];
    }



    UTOPIA_INLINE_FUNCTION static void eval_uniform_J(const Scalar &h, Array<Scalar, 4> &J) {



      //gfn 0
        J.values[0] = h;
        J.values[1] = 0.0;
        J.values[2] = 0.0;
        J.values[3] = h;

        
    }

    UTOPIA_INLINE_FUNCTION static void eval_uniform_inv_J(const Scalar &h, Array<Scalar, 4> &J) {


        J.values[0] = 1./h;
        J.values[1] = 0.0;
        J.values[2] = 0.0;
        J.values[3] = 1./h;

        
    }


         /*evaluation jacobian transformation*/
    UTOPIA_INLINE_FUNCTION static void eval_uniform_det_J(const Scalar &h, Scalar &J) {

        J= h*h;

        
    }

    UTOPIA_INLINE_FUNCTION Scalar static measure ()
     {
       return 1.0;
     }

    private:

      Array<Point, 4> nodes_;

  };


  template<typename T>
  class Hex8 {

    using Scalar = T;

    using Point= Array<T, 3>;
  
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
    UTOPIA_INLINE_FUNCTION static void eval_uniform_J(const Scalar &h, Array<Scalar, 9> &J) {



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
    UTOPIA_INLINE_FUNCTION static void eval_uniform_inv_J(const Scalar &h, Array<Scalar, 9> &J) {



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
    UTOPIA_INLINE_FUNCTION static void eval_uniform_det_J(const Scalar &h, Scalar &J) {



      //gfn 0
        J= h*h*h;
        
    }


    UTOPIA_INLINE_FUNCTION void eval_J(Array<Scalar, 9> &J) {


          for(int i = 0; i < 3; ++i) {
                const int i_offset = i * 3;

                for(int j = 1; j < 4; ++j) {
                    J[i_offset + j-1] = node(j)[i] - node(0)[i];
                }
            }

        
    }

         /*evaluation jacobian transformation*/
    UTOPIA_INLINE_FUNCTION void eval_det_J(Scalar& detJ) {

        Array<Scalar, 9> J;

          for(int i = 0; i < 3; ++i) {

            const int i_offset = i * 3;

            for(int j = 1; j < 4; ++j) {
                
                J[i_offset + j-1] = node(j)[i] - node(0)[i];
              
              }
          }

        Scalar temp0 = J[8] * J[4] - J[7] * J[5];
        Scalar temp1 = J[8] * J[1] - J[7] * J[2];
        Scalar temp2 = J[5] * J[1] - J[4] * J[2];

        detJ = J[0] * temp0 - J[1] * temp1 + J[2]*temp2;

    }

    UTOPIA_INLINE_FUNCTION void eval_inv_J(Array<Scalar, 9> invJ) {
        
        Array<Scalar, 9> J;

          for(int i = 0; i < 3; ++i) {

            const int i_offset = i * 3;

            for(int j = 1; j < 4; ++j) {
                
                J[i_offset + j-1] = node(j)[i] - node(0)[i];
              
              }
          }

        Scalar temp0 = J[8]*J[4] - J[7]*J[5];
        Scalar temp1 = J[8]*J[1] - J[7]*J[2];
        Scalar temp2 = J[5]*J[1] - J[4]*J[2];

        Scalar detJ = J[0] * temp0 - J[1] * temp1 + J[2]*temp2;

        Scalar inv_detJ = 1.0/detJ;

              
        invJ[0] =   1.0 * temp0 * inv_detJ;
        invJ[1] = - 1.0 * temp1 * inv_detJ;
        invJ[2] =   1.0 * temp2 * inv_detJ;

        invJ[3] = -1.0 * ( J[4] * J[3] - J[6] * J[5] ) * inv_detJ;
        invJ[4] =  ( J[8] * J[0] - J[6] * J[2] ) * inv_detJ;
        invJ[5] = -1.0 * ( J[5] * J[0] - J[3] * J[2] ) * inv_detJ;

        invJ[6] =  ( J[7] * J[3] -  J[6] * J[4] ) * inv_detJ;
        invJ[7] = -1.0 * ( J[7] * J[0] -  J[6] * J[1] ) * inv_detJ;
        invJ[8] =  ( J[3] * J[0] -  J[3] * J[1] ) *inv_detJ;

    }

    Point &node(const int i) override
    {
        return nodes_[i];
    }

    const Point &node(const int i) const override
    {
        return nodes_[i];
    }



   UTOPIA_INLINE_FUNCTION Scalar static measure ()
   {
     return 1.0;
   }

   private:

      Array<Point, 8> nodes_;



  };






  

}
#endif


