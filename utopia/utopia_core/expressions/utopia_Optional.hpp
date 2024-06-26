#ifndef UTOPIA_OPTIONAL_HPP
#define UTOPIA_OPTIONAL_HPP

#include <utility>

namespace utopia {

    template <class Args, int Begin, int End>
    class OptionalIterator {
    public:
        template <class Fun>
        void visit(Fun &fun) {
            fun.parse_arg(args.template get<Begin>());
            OptionalIterator<Args, Begin + 1, End> next(args);
            next.visit(fun);
        }

        OptionalIterator(const Args &args) : args(args) {}
        const Args &args;
    };

    template <class Args, int Begin>
    class OptionalIterator<Args, Begin, Begin> {
    public:
        template <class Fun>
        void visit(Fun &) {}
        OptionalIterator(const Args &) {}
    };

    template <class... Args>
    class Optional {
    public:
        static const int n_args = std::tuple_size<std::tuple<Args...>>::value;

        Optional(const Args &... args) : opts(args...) {}

        template <int Index>
        inline auto get() const -> const typename std::tuple_element<Index, std::tuple<Args...>>::type {
            return std::get<Index>(opts);
        }

        template <class Fun>
        void each(Fun &fun) const {
            OptionalIterator<const Optional, 0, n_args> iter(*this);
            iter.visit(fun);
        }

        std::tuple<Args...> opts;
    };

    template <class... Args>
    Optional<Args...> options(const Args &... args) {
        return Optional<Args...>(args...);
    }
}  // namespace utopia

#endif  // UTOPIA_OPTIONAL_HPP
