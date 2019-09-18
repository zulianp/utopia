#ifndef UTOPIA_TRILINOS_EACH_HPP
#define UTOPIA_TRILINOS_EACH_HPP

#include "utopia_Traits.hpp"
#include "utopia_Each.hpp"
#include "utopia_trilinos_ForwardDeclaration.hpp"
#include "utopia_trilinos_Traits.hpp"

namespace utopia {

    class TpetraMatrixEach {
    public:
        using SizeType = typename Traits<TpetraMatrix>::SizeType;
        using Scalar   = typename Traits<TpetraMatrix>::Scalar;

        template<class Fun>
        static void apply(TpetraMatrix &mat, Fun fun);

        template<class Fun>
        static void apply_transform(TpetraMatrix &mat, Fun fun);

        static void apply_read(const TpetraMatrix &mat, std::function<void(const Scalar &)> &fun);

        template<class Fun>
        static void apply_read(const TpetraMatrix &mat, Fun fun);
    };

    class TpetraVectorEach {
    public:
        using SizeType = typename Traits<TpetraVector>::SizeType;
        using Scalar   = typename Traits<TpetraVector>::Scalar;

        static void apply_read(const TpetraVector &v, std::function<void(const Scalar &)> &fun);

        template<class Fun>
        static void apply_read(const TpetraVector &v, Fun fun);

        template<class Fun>
        static void apply_write(TpetraVector &v, Fun fun);

        template<class Fun>
        static void apply_transform(const TpetraVector &in, TpetraVector &out, Fun fun);
    };

    template<>
    class Each<TpetraMatrix, 2, FillType::SPARSE> : public TpetraMatrixEach {};

    template<int FILL_TYPE>
    class Each<TpetraVector, 1, FILL_TYPE> : public TpetraVectorEach {};

    template<class Fun>
    inline void each_transform(TpetraVector &in_out, Fun fun)
    {
        each_transform(in_out, in_out, fun);
    }

}

#endif //UTOPIA_TRILINOS_EACH_HPP
