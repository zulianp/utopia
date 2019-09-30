#ifndef UTOPIA_TRILINOS_DEVICEVIEW_HPP
#define UTOPIA_TRILINOS_DEVICEVIEW_HPP

#include "utopia_Traits.hpp"
#include "utopia_Readable.hpp"
#include "utopia_DeviceView.hpp"

namespace utopia {

    // class DeviceView<TpetraVector, 1> {
    // public:
    //     using Scalar   = typename Traits<T>::Scalar;
    //     using SizeType = typename Traits<T>::SizeType;

    //     inline Scalar get(const SizeType &idx) const
    //     {
    //         return tensor_.get(idx);
    //     }

    //     inline void set(const SizeType &idx, const Scalar &value) const
    //     {
    //         return tensor_.set(idx, value);
    //     }

    //     inline void add(const SizeType &idx, const Scalar &value) const
    //     {
    //         return tensor_.add(idx, value);
    //     }

    //     DeviceView(T &tensor) : tensor_(tensor), lock_(tensor) {}

    // private:
    //     T &tensor_;
    //     ReadAndWrite<T> lock_;
    // };

    template<>
    class DeviceView<const TpetraVector, 1> {
    public:
        using Scalar          = typename Traits<TpetraVector>::Scalar;
        using SizeType        = typename Traits<TpetraVector>::SizeType;
        using ExecutionSpaceT = typename TpetraVector::vector_type::execution_space;
        using DualViewType    = typename TpetraVector::vector_type::dual_view_type;
        using DeviceViewType  = typename DualViewType::t_dev;
        using LocalMapType    = typename TpetraVector::vector_type::map_type::local_map_type;

        inline Scalar get(const SizeType &idx) const
        {
            auto local_idx = map_.getLocalElement(idx);
            return view_(local_idx, 0);
        }

        DeviceView(const TpetraVector &tensor) : 
        tensor_(tensor),
        view_(tensor.raw_type()->template getLocalView<ExecutionSpaceT>()) 
        {}

    private:
        const TpetraVector &tensor_;
        DeviceViewType view_;
        LocalMapType   map_;
    };


    template<>
    class DeviceView<TpetraVector, 1> {
    public:
        using Scalar          = typename Traits<TpetraVector>::Scalar;
        using SizeType        = typename Traits<TpetraVector>::SizeType;
        using ExecutionSpaceT = typename TpetraVector::vector_type::execution_space;
        using DualViewType    = typename TpetraVector::vector_type::dual_view_type;
        using DeviceViewType  = typename DualViewType::t_dev;
        using LocalMapType    = typename TpetraVector::vector_type::map_type::local_map_type;

        inline Scalar get(const SizeType &idx) const
        {
            auto local_idx = map_.getLocalElement(idx);
            return view_(idx, 0);
        }


        inline void set(const SizeType &idx, const Scalar &value) const
        {
            auto local_idx = map_.getLocalElement(idx);
            view_(idx, 0) = value;
        }

        inline void add(const SizeType &idx, const Scalar &value) const
        {
            view_(idx, 0) += value;
        }

        DeviceView(const TpetraVector &tensor) : 
        tensor_(tensor),
        view_(tensor.raw_type()->template getLocalView<ExecutionSpaceT>()) 
        {}

    private:
        const TpetraVector &tensor_;
        DeviceViewType view_;
        LocalMapType   map_;
    };


}

#endif //UTOPIA_TRILINOS_DEVICEVIEW_HPP