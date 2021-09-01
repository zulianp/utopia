#ifndef UTOPIA_VELOCITY_VARIANT_HPP
#define UTOPIA_VELOCITY_VARIANT_HPP

namespace utopia {
    template <class FunctionSpace>
    class VelocityVariant {
    public:
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        virtual ~VelocityVariant() = default;
    };

    // template <class FunctionSpace>
    // VelocityVariant<FunctionSpace> &velocity_variant(TimeDependentFunction<FunctionSpace> &function)
    // {

    // }

}  // namespace utopia

#endif  // UTOPIA_VELOCITY_VARIANT_HPP
