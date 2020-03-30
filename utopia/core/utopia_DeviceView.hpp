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

        inline void atomic_add(const SizeType &idx, const Scalar &value) const
        {
            //FIXME
            tensor_.c_add(idx, value);
        }

        //FIXME is not atomic
        template<class Index, class Values>
        inline void atomic_add_vector(const Index &I,  const Values &V) const
        {
            tensor_.add_vector(I, V);
        }


        inline void atomic_set(const SizeType &idx, const Scalar &value) const
        {
            //FIXME
            tensor_.c_set(idx, value);
        }

        DeviceView(T &tensor, WriteMode wm = utopia::AUTO) : tensor_(tensor), lock_(std::make_shared<ReadAndWrite<T>>(tensor, wm)) {}

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


    template<class T, int Order>
    class LocalViewDevice {};

    template<class T>
    class LocalViewDevice<const T, 1> {
    public:
        using Scalar   = typename Traits<T>::Scalar;
        using SizeType = typename Traits<T>::SizeType;

        inline Scalar get(const SizeType &idx) const
        {
            return tensor_.l_get(idx);
        }

        LocalViewDevice(const T &tensor) : tensor_(tensor), lock_(std::make_shared<Read<T>>(tensor)) {}

    private:
        const T &tensor_;
        std::shared_ptr<Read<T>> lock_;
    };

    template<class T>
    class LocalViewDevice<T, 1> {
    public:
        using Scalar   = typename Traits<T>::Scalar;
        using SizeType = typename Traits<T>::SizeType;

        inline Scalar get(const SizeType &idx) const
        {
            return tensor_.l_get(idx);
        }

        inline void set(const SizeType &idx, const Scalar &val) const
        {
            return tensor_.l_set(idx, val);
        }

        inline void add(const SizeType &idx, const Scalar &val) const
        {
            return tensor_.l_add(idx, val);
        }

        //FIXME
        inline void atomic_add(const SizeType &idx, const Scalar &val) const
        {
            return tensor_.l_add(idx, val);
        }

        LocalViewDevice(T &tensor) : tensor_(tensor), lock_(std::make_shared<ReadAndWrite<T>>(tensor)) {}

    private:
        T &tensor_;
        std::shared_ptr<ReadAndWrite<T>> lock_;
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
            tensor_.c_add(i, j, value);
        }

        //FIXME is not atomic
        template<class Index, class Values>
        inline void atomic_add_matrix(const Index &I, const Index &J, const Values &V) const
        {
            tensor_.add_matrix(I, J, V);
        }

        DeviceView(T &tensor, WriteMode wm = utopia::AUTO) : tensor_(tensor), lock_(std::make_shared<Write<T>>(tensor, wm)) {}

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


    template<class Derived, int Order>
    inline DeviceView<Derived, Order> view_device(Tensor<Derived, Order> &t)
    {
        return t.derived();
    }

    template<class Derived, int Order>
    inline DeviceView<const Derived, Order> view_device(const Tensor<Derived, Order> &t)
    {
        return t.derived();
    }

    template<class Derived, int Order>
    inline DeviceView<const Derived, Order> const_view_device(const Tensor<Derived, Order> &t)
    {
        return t.derived();
    }

    template<class Derived, int Order>
    inline LocalViewDevice<Derived, Order> local_view_device(Tensor<Derived, Order> &t)
    {
        return t.derived();
    }

    template<class Derived, int Order>
    inline LocalViewDevice<const Derived, Order> const_local_view_device(const Tensor<Derived, Order> &t)
    {
        return t.derived();
    }

}

#endif //UTOPIA_DEVICEVIEW_HPP