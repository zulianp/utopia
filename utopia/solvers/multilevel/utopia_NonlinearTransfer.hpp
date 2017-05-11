// /*
// * @Author: alenakopanicakova
// * @Date:   2017-05-09
// * @Last Modified by:   Alena Kopanicakova
// * @Last Modified time: 2017-05-09
// */

// #ifndef UTOPIA_ML_NONLINEAR_TRANSFER_HPP
// #define UTOPIA_ML_NONLINEAR_TRANSFER_HPP

//      namespace utopia 
//      {
//         /**
//          * @brief      The class for transfer operators. 
//          *
//          * @tparam     Matrix  
//          * @tparam     Vector  
//          */
//         template<class Matrix, class Vector>
//         class NonLinearTransfer : public Transfer<Matrix, Vector>
//         {
//             typedef UTOPIA_SCALAR(Vector)    Scalar;
//             typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


//         public:

//             NonLinearTransfer(const Matrix & I, const Matrix & P, const Vector & x):
//                                     Transfer<Matrix, Vector>(I),
//                                     _P(P), 
//                                     _x_init(x)
//             {

//             }

//             virtual ~NonLinearTransfer(){} 
            
//             /*=====================================================
//                                 initialization
//             =====================================================*/
            

//             /**
//              * @brief      Initialization of projection down operator.
//              *
//              * @param[in]  P_in  The projection operator. 
//              *
//              */
//             virtual bool P_init(const Matrix &P_in)
//             {
//                 _P = P_in; 
//                 return true; 
//             }


//             /**
//              * @brief      Initialization of projection down operator.
//              *
//              * @param[in]  P_in  The projection operator. 
//              *
//              */
//             virtual bool boundary_values_init(const Vector &x_in)
//             {
//                 _x_init = x_in; 
//                 return true; 
//             }



//             /*=====================================================
//                                     actions
//             =====================================================*/

//             /**
//              * @brief      Projection of vector 
//              *            \f$  x_{new} = P * x  \f$
//              * @param[in]  x     
//              * @param      x_new 
//              *
//              */
//             virtual bool project_down(const Vector &x, Vector &x_new)
//             {
//                 if(local_size(_P).get(1)==local_size(x).get(0))
//                 {
//                     std::cout<<"YES projections down ... \n"; 
//                     x_new = _P * x; 
//                     return true; 
//                 }
//                 else
//                     return(Transfer<Matrix, Vector>::restrict(x, x_new)); 
//             }


//             virtual bool get_boundary_values(Vector & x)
//             {   
//                 x = _x_init; 
//                 return true; 
//             }




//         protected:        
//             Matrix  _P;  
//             Vector _x_init;   // vector carring on values on boundary



//     };

// }

// #endif //UTOPIA_ML_NONLINEAR_TRANSFER_HPP

