#ifndef UTOPIA_RANGE_DEVICE_HPP
#define UTOPIA_RANGE_DEVICE_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Range.hpp"
#include "utopia_Wrapper.hpp"

namespace utopia {

    template<class T>
    class RangeDevice  {
    public:
        using SizeType = typename Traits<T>::SizeType;
        using Device   = typename Traits<T>::Device;

        UTOPIA_INLINE_FUNCTION SizeType begin() const
        {
            return begin_;
        }

        UTOPIA_INLINE_FUNCTION SizeType end() const
        {
            return end_;
        }

        UTOPIA_INLINE_FUNCTION SizeType extent()
        {
            return end() - begin();
        }

        UTOPIA_INLINE_FUNCTION RangeDevice(const SizeType &begin, const SizeType &end)
        : begin_(begin), end_(end)
        {}

    private:
        SizeType begin_, end_;
    };

    template<class T>
    inline RangeDevice<T> range_device(const T &t)
    {
        Range r = range(t);
        return RangeDevice<T>(r.begin(), r.end());
    }

    template<class T>
    inline RangeDevice<T> local_range_device(const T &t)
    {
        Range r = range(t);
        return RangeDevice<T>(0, r.extent());
    }

    template<class T, typename F>
    inline static void parallel_for(const RangeDevice<T> &r, F f)
    {
        using Device = typename RangeDevice<T>::Device;
        Device::template parallel_for(r.begin(), r.end(), f);
    }

    template<class T, typename... Args>
    inline static void parallel_reduce(const RangeDevice<T> &r, Args &&...args)
    {
        using Device = typename RangeDevice<T>::Device;
        Device::template parallel_reduce(r.begin(), r.end(), std::forward<Args>(args)...);
    }
}

#endif //UTOPIA_RANGE_DEVICE_HPP