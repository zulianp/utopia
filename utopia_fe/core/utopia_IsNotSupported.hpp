#ifndef UTOPIA_IS_NOT_SUPPORTED_HPP
#define UTOPIA_IS_NOT_SUPPORTED_HPP

namespace utopia {

    template <class Type>
    class IsNotSupported {
    public:
        static constexpr bool value{false};
    };

}  // namespace utopia

#endif