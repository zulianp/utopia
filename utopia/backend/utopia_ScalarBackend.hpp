#ifndef UTOPIA_SCALAR_BACKEND_HPP
#define UTOPIA_SCALAR_BACKEND_HPP

namespace utopia {

    template<typename Scalar>
    class ScalarBackend {
    public:

        static inline void assign(Scalar &result, const Scalar &value)
        {
            result = value;
        }

        static inline void apply_binary(Scalar &result, const Scalar left, const Plus &op, const Scalar right)
        {
            aux_apply_binary(left, right, op, result);
        }

        static inline void apply_binary(Scalar &result, const Scalar left, const Minus &op, const Scalar right)
        {
            aux_apply_binary(left, right, op, result);
        }

        static inline void apply_binary(Scalar &result, const Scalar left, const Multiplies &op, const Scalar right)
        {
            aux_apply_binary(left, right, op, result);
        }

        static inline void apply_binary(Scalar &result, const Scalar left, const EMultiplies &op, const Scalar right)
        {
            aux_apply_binary(left, right, op, result);
        }

        static inline void apply_binary(Scalar &result, const Scalar left, const Divides &op, const Scalar right)
        {
            aux_apply_binary(left, right, op, result);
        }

        static inline void apply_binary(Scalar &result, const Scalar left, const AbsPlus &op, const Scalar right)
        {
            aux_apply_binary(left, right, op, result);
        }

        static inline void apply_binary(bool &result, const Scalar left, const Scalar right, const ApproxEqual &op)
        {
            aux_apply_binary(left, right, op, result);
        }

        static inline void axpy(Scalar &y, const Scalar alpha, const Scalar x)
        {
            y += alpha * x;
        }

    protected:
        ScalarBackend() {}


        template<typename Left, typename Right, typename Operation, typename Result>
        static inline void aux_apply_binary(const Left &left, const Right &right, const Operation &op, Result &result)
        {
            result = op.template apply<Right>(left, right);
        }
    };

}

#endif //UTOPIA_SCALAR_BACKEND_HPP
