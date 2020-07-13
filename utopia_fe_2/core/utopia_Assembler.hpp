#ifndef UTOPIA_ASSEMBLER_HPP
#define UTOPIA_ASSEMBLER_HPP

namespace utopia {

    class IAssembler {
    public:
        virtual ~IAssembler() = default;
    };

    template <typename...>
    class Assembler {};

}  // namespace utopia

#endif  // UTOPIA_ASSEMBLER_HPP
