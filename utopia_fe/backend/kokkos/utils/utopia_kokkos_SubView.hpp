#ifndef UTOPIA_KOKKOS_SUBVIEW_HPP
#define UTOPIA_KOKKOS_SUBVIEW_HPP

namespace utopia {
    namespace kokkos {

        class IdentityRange {
        public:
            UTOPIA_INLINE_FUNCTION static constexpr int begin() { return 0; }
            UTOPIA_INLINE_FUNCTION static constexpr int index(const int idx) { return idx; }
        };

        template <int Begin, int End>
        class StaticRange {
        public:
            UTOPIA_INLINE_FUNCTION static constexpr int begin() { return Begin; }
            UTOPIA_INLINE_FUNCTION static constexpr int end() { return End; }
            UTOPIA_INLINE_FUNCTION static constexpr int index(const int idx) { return begin() + idx; }
        };

        template <class OpOrView, class... Ranges>
        class SubView {
        public:
            OpOrView view;
        };

        template <class OpOrView, class Range0, class Range1, class Range2, class Range3>
        class SubView<OpOrView, Range0, Range1, Range2, Range3> {
        public:
            using Scalar = typename Traits<OpOrView>::Scalar;

            UTOPIA_INLINE_FUNCTION constexpr SubView(const OpOrView &view,
                                                     const Range0 &s0 = Range0(),
                                                     const Range1 &s1 = Range1(),
                                                     const Range2 &s2 = Range2(),
                                                     const Range3 &s3 = Range3())
                : view(view), s0(s0), s1(s1), s2(s2), s3(s3) {}

            UTOPIA_INLINE_FUNCTION constexpr Scalar operator()(const int cell,
                                                               const int qp,
                                                               const int sub_i,
                                                               const int sub_j) const {
                return view(s0.index(cell), s1.index(qp), s2.index(sub_i), s3.index(sub_j));
            }

            OpOrView view;
            Range0 s0;
            Range1 s1;
            Range2 s2;
            Range3 s3;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_SUBVIEW_HPP