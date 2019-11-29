#ifndef UTOPIA_DEVICEVIEW_HPP
#define UTOPIA_DEVICEVIEW_HPP

#include "utopia_Traits.hpp"
#include "utopia_Readable.hpp"
#include <memory>

namespace utopia {

    template<class T, int Order = Traits<T>::Order>
    class DeviceView {};

    template<class T>
    class DeviceView<T, 1> {
    public:
        using Scalar   = typename Traits<T>::Scalar;
        using SizeType = typename Traits<T>::SizeType;

        inline Scalar get(const SizeType &idx) const
        {
            return tensor_.get(idx);
        }

        inline void set(const SizeType &idx, const Scalar &value) const
        {
            tensor_.set(idx, value);
        }

        inline void add(const SizeType &idx, const Scalar &value) const
        {
            tensor_.add(idx, value);
        }

        DeviceView(T &tensor) : tensor_(tensor), lock_(std::make_shared<ReadAndWrite<T>>(tensor)) {}

    private:
        T &tensor_;
        std::shared_ptr<ReadAndWrite<T>> lock_;
    };

    template<class T>
    class DeviceView<const T, 1> {
    public:
        using Scalar   = typename Traits<T>::Scalar;
        using SizeType = typename Traits<T>::SizeType;

        inline Scalar get(const SizeType &idx) const
        {
            return tensor_.get(idx);
        }

        DeviceView(const T &tensor) : tensor_(tensor), lock_(std::make_shared<Read<T>>(tensor)) {}

    private:
        const T &tensor_;
        std::shared_ptr<Read<T>> lock_;
    };

    ///////////////////////////////////////////////////////////////////////////////////////

    template<class T>
    class DeviceView<T, 2> {
    public:
        using Scalar   = typename Traits<T>::Scalar;
        using SizeType = typename Traits<T>::SizeType;

        //FIXME is not atomic
        inline void atomic_add(const SizeType &i, const SizeType &j, const Scalar &value) const
        {
            tensor_.add(i, j, value);
        }

        DeviceView(T &tensor) : tensor_(tensor), lock_(std::make_shared<Write<T>>(tensor)) {}

    private:
        T &tensor_;
        std::shared_ptr<Write<T>> lock_;
    };

    ///////////////////////////////////////////////////////////////////////////////////////

    template<class Derived, int Order>
    inline DeviceView<Derived, Order> device_view(Tensor<Derived, Order> &t)
    {
        return t.derived();
    }

    template<class Derived, int Order>
    inline DeviceView<const Derived, Order> device_view(const Tensor<Derived, Order> &t)
    {
        return t.derived();
    }

    template<class Derived, int Order>
    inline DeviceView<const Derived, Order> const_device_view(const Tensor<Derived, Order> &t)
    {
        return t.derived();
    }

}

#endif //UTOPIA_DEVICEVIEW_HPP