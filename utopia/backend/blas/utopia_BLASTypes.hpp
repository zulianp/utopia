
#ifndef UTOPIA_BLAS_TYPES_HPP
#define UTOPIA_BLAS_TYPES_HPP

#include "utopia_BLASMatrix.hpp"
#include "utopia_BLASTraits.hpp"
#include "utopia_CRSMatrix.hpp"

#include "utopia_Wrapper.hpp"

namespace utopia {
    typedef utopia::Wrapper<utopia::Matrix<double>, 2>    Matrixd;
    typedef utopia::Wrapper<std::vector<double>, 1>       Vectord;
    typedef utopia::Wrapper<utopia::CRSMatrix<double>, 2> CRSMatrixd;
    typedef utopia::Wrapper<utopia::CCSMatrix<double>, 2> CCSMatrixd;

    template<>
    class Write<CRSMatrixd> {
    public:
        Write(CRSMatrixd &v)
        : _v(v) {

            _v.implementation().assemblyBegin();

        }

        ~Write() {
            _v.implementation().assemblyEnd();
        }

    private:
        CRSMatrixd &_v;
    };


    template<>
    class Write<CCSMatrixd> {
    public:
        Write(CCSMatrixd &v)
        : _v(v) {

            _v.implementation().assemblyBegin();

        }

        ~Write() {
            _v.implementation().assemblyEnd();
        }

    private:
        CCSMatrixd &_v;
    };


        ///////////////////////////////////////////////////////////////////////////////

    inline Matrix<double> &raw_type(Wrapper<Matrix<double>, 2> &utopiaType)
    {
        return utopiaType.implementation();
    }

    inline std::vector<double> &raw_type(Wrapper<std::vector<double>, 1> &utopiaType)
    {
        return utopiaType.implementation();
    }

        ///////////////////////////////////////////////////////////////////////////////

    inline const Matrix<double> &raw_type(const Wrapper<Matrix<double>, 2> &utopiaType)
    {
        return utopiaType.implementation();
    }

    inline const std::vector<double> &raw_type(const Wrapper<std::vector<double>, 1> &utopiaType)
    {
        return utopiaType.implementation();
    }

    // #define UTOPIA_MAKE_EVAL(Tensor_)                      \
    // template<> class Eval<Tensor_,                         \
    //                       utopia::Traits<Tensor_>,         \
    //                       utopia::Traits<Tensor_>::Backend \
    //                       > : public Eval<                 \
    //                                     utopia::Wrapper<Tensor_::Implementation, Tensor_::Order>, \
    //                                     utopia::Traits< Wrapper<Tensor_::Implementation, Tensor_::Order> >,                  \
    //                                     utopia::Traits<Tensor_>::Backend> {};



        // UTOPIA_MAKE_EVAL(Vectord);
}

#endif //UTOPIA_BLAS_TYPES_HPP
