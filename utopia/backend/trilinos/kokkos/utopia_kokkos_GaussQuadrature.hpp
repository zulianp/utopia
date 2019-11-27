#ifndef UTOPIA_KOKKOS_GaussQuadrature_hpp
#define UTOPIA_KOKKOS_GaussQuadrature_hpp


#include "utopia_Base.hpp"


namespace utopia {

  template<typename T, int Dim, int Order>
  class QGauss {
      using Scalar = T;
  };

  template<typename Scalar>
  class QGauss<Scalar, 2, 1> {
    public:
        static constexpr const Scalar a   = 0.211324865405187;
        static constexpr const Scalar b   = 0.788675134594813;
        static constexpr const Scalar wts = 0.25;
        static constexpr const int N = 4;
    

      UTOPIA_INLINE_FUNCTION constexpr QGauss()
        : points{
            a,a,
            b,a,
            a,b,
            b,b
        },

          weights{
            wts, wts, wts, wts
        }

    {}

        Scalar points[N * 2];
        Scalar weights[N];

  };


  template<typename Scalar>
  class QGauss<Scalar, 2, 2> {
    public:
        static constexpr const Scalar a   = 0.5;
        static constexpr const Scalar b   = 0.98304589153964795245728880523899;
        static constexpr const Scalar c   = 0.72780186391809642112479237299488;
        static constexpr const Scalar d   = 0.13418502421343273531598225407969;
        static constexpr const Scalar e   = 0.92595732665230024565091782018333;
        static constexpr const Scalar g   = 0.074042673347699754349082179816666; 
        static constexpr const Scalar h   = 0.18454360551162298687829339850317;
        static constexpr const Scalar i   = 0.81454360551162298687829339850317; 
        static constexpr const Scalar wts_0   = 0.28571428571428571428571428571428;
        static constexpr const Scalar wts_1   = 0.10989010989010989010989010989011;
        static constexpr const Scalar wts_2   = 0.14151805175188302631601261486295;
        static constexpr const Scalar wts_3   = 0.16067975044591917148618518733485;
        static constexpr const int N = 6;
    

      UTOPIA_INLINE_FUNCTION constexpr QGauss()
        : points{
            a,a,
            b,a,
            c,g,           
            c,e,            
            d,h,
            d,i
        },

          weights{
            wts_0, wts_1, wts_2, wts_2, wts_3, wts_3
        }

    {}

        Scalar points[N * 2];
        Scalar weights[N];

  };

  template<typename Scalar>
    class QGauss<Scalar, 3, 1> {
      public:
          static constexpr const Scalar a   = 0.0;
          static constexpr const Scalar b   = 0.5;
          static constexpr const Scalar c   = 1.0;
          static constexpr const Scalar wts_0   =  0.16666666666666666666666666666667;
          static constexpr const int N = 6;
      

        UTOPIA_INLINE_FUNCTION constexpr QGauss()
          : points{
              a,b,b,
              b,a,b,
              b,b,a,           
              b,b,c,           
              b,c,b,
              c,b,b,
          },

            weights{
              wts_0, wts_0, wts_0, wts_0, wts_0, wts_0
          }

      {}

          Scalar points[N * 3];
          Scalar weights[N];

    };

  template<typename Scalar>
    class QGauss<Scalar, 3, 4> {
      public:
          static constexpr const Scalar a   = 0.112701665379258;
          static constexpr const Scalar b   = 0.5;
          static constexpr const Scalar c   = 0.887298334620742;
          static constexpr const Scalar wts_0   =   0.021433470507545;
          static constexpr const Scalar wts_1   =   0.034293552812071;
          static constexpr const Scalar wts_2   =   0.054869684499314;
          static constexpr const Scalar wts_3   =   0.087791495198903;
          static constexpr const int N = 27;
      

        UTOPIA_INLINE_FUNCTION constexpr QGauss()
          : points{
              a,a,a,
              b,a,a,
              c,a,a,          
              a,b,a,
              b,b,a,
              c,b,a,
              a,c,a,
              b,c,a,
              c,c,a,
              a,a,b,
              b,a,b,
              c,a,b,
              a,b,b,
              b,b,b,
              c,b,b,
              a,c,b,
              b,c,b,
              c,c,b,
              a,a,c,
              b,a,c,
              c,a,c,
              a,b,c,
              b,b,c,
              c,b,c,
              a,c,c,
              b,c,c,
              c,c,c

          },

            weights{
              wts_0, 
              wts_1,
              wts_0, 
              wts_1, 
              wts_2,
              wts_1,
              wts_0,
              wts_1,
              wts_0,
              wts_1,
              wts_2,
              wts_1,
              wts_2,
              wts_3,
              wts_2,
              wts_1,
              wts_2,
              wts_1,
              wts_0,
              wts_1,
              wts_0,
              wts_1,
              wts_2,
              wts_1,
              wts_0,
              wts_2,
              wts_1

          }

      {}

          Scalar points[N * 3];
          Scalar weights[N];

    };

}

#endif

