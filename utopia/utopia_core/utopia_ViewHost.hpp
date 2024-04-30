#ifndef UTOPIA_VIEW_HOST
#define UTOPIA_VIEW_HOST

#include "utopia_DeviceView.hpp"

// FIXME If we ever support PETSc GPU the following code will create problems (i.e., specialization is required)

namespace utopia {

    template <typename T, int Order = Traits<T>::Order>
    class ViewHost : public DeviceView<T, Order> {
    public:
        using Super = utopia::DeviceView<T, Order>;
        using Super::Super;
    };

    template <typename T, int Order = Traits<T>::Order>
    class LocalViewHost : public LocalViewDevice<T, Order> {
    public:
        using Super = utopia::LocalViewDevice<T, Order>;
        using Super::Super;
    };

    template <class Derived, int Order>
    inline ViewHost<Derived, Order> view_host(Tensor<Derived, Order> &t) {
        return t.derived();
    }

    template <class Derived, int Order>
    inline ViewHost<const Derived, Order> view_host(const Tensor<Derived, Order> &t) {
        return t.derived();
    }

    template <class Derived, int Order>
    inline ViewHost<const Derived, Order> const_view_host(const Tensor<Derived, Order> &t) {
        return t.derived();
    }

    template <class Derived, int Order>
    inline LocalViewHost<Derived, Order> local_view_host(Tensor<Derived, Order> &t) {
        return t.derived();
    }

    template <class Derived, int Order>
    inline LocalViewHost<const Derived, Order> local_view_host(const Tensor<Derived, Order> &t) {
        return t.derived();
    }

    template <class Derived, int Order>
    inline LocalViewHost<const Derived, Order> const_local_view_host(const Tensor<Derived, Order> &t) {
        return t.derived();
    }
}  // namespace utopia

#endif  // UTOPIA_VIEW_HOST