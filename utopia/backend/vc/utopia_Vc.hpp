#ifndef UTOPIA_VC_HPP
#define UTOPIA_VC_HPP

#include <Vc/Vc>

namespace utopia {
    namespace device {
        template <typename T>
        inline Vc::Vector<T> min(const Vc::Vector<T> &left, const Vc::Vector<T> &right) {
            return std::min(left, right);
        }

        template <typename T>
        inline Vc::Vector<T> max(const Vc::Vector<T> &left, const Vc::Vector<T> &right) {
            return std::max(left, right);
        }

        template <typename T>
        inline Vc::Vector<T> abs(const Vc::Vector<T> &v) {
            return std::abs(v);
        }
    }  // namespace device

}  // namespace utopia

#endif  // UTOPIA_VC_HPP
