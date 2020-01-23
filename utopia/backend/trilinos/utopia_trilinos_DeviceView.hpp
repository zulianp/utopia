#ifndef UTOPIA_TRILINOS_DEVICEVIEW_HPP
#define UTOPIA_TRILINOS_DEVICEVIEW_HPP

#include "utopia_Traits.hpp"
#include "utopia_Readable.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_Algorithms.hpp"

namespace utopia {

    template<>
    class DeviceView<const TpetraVector, 1> {
    public:
        using Scalar          = typename Traits<TpetraVector>::Scalar;
        using SizeType        = typename Traits<TpetraVector>::SizeType;
        using ExecutionSpaceT = typename TpetraVector::ExecutionSpace;
        using DualViewType    = typename TpetraVector::VectorType::dual_view_type;
        using DeviceViewType  = typename DualViewType::t_dev;
        using LocalMapType    = typename TpetraVector::MapType::local_map_type;

        UTOPIA_INLINE_FUNCTION Scalar get(const SizeType &idx) const
        {
            auto local_idx = map_.getLocalElement(idx);
            return view_(local_idx, 0);
        }

        DeviceView(const TpetraVector &tensor) :
        view_(tensor.raw_type()->template getLocalView<ExecutionSpaceT>())
        {
            assert(!tensor.has_ghosts() && "GHOST HANDLING NOT IMPLEMENTED YET");
            map_ = tensor.raw_type()->getMap()->getLocalMap();
        }

    private:
        DeviceViewType view_;
        LocalMapType   map_;
    };

    template<>
    class DeviceView<TpetraVector, 1> {
    public:
        using Scalar          = typename Traits<TpetraVector>::Scalar;
        using SizeType        = typename Traits<TpetraVector>::SizeType;
        using ExecutionSpaceT = typename TpetraVector::ExecutionSpace;
        using DualViewType    = typename TpetraVector::VectorType::dual_view_type;
        using DeviceViewType  = typename DualViewType::t_dev;
        using LocalMapType    = typename TpetraVector::MapType::local_map_type;

        UTOPIA_INLINE_FUNCTION Scalar get(const SizeType &idx) const
        {
            auto local_idx = map_.getLocalElement(idx);
            return view_(local_idx, 0);
        }

        UTOPIA_INLINE_FUNCTION void set(const SizeType &idx, const Scalar &value) const
        {
            auto local_idx = map_.getLocalElement(idx);
            view_(local_idx, 0) = value;
        }

        UTOPIA_INLINE_FUNCTION void add(const SizeType &idx, const Scalar &value) const
        {
            auto local_idx = map_.getLocalElement(idx);
            view_(local_idx, 0) += value;
        }

        /**
         * Atomic with respect to the value, no the map or idx
         */
        UTOPIA_INLINE_FUNCTION void atomic_add(const SizeType &idx, const Scalar &value) const
        {
            auto local_idx = map_.getLocalElement(idx);
            device::atomic_add(&view_(local_idx, 0), value);
        }

        DeviceView(const TpetraVector &tensor) :
        view_(tensor.raw_type()->template getLocalView<ExecutionSpaceT>())
        {
            assert(!tensor.has_ghosts() && "GHOST HANDLING NOT IMPLEMENTED YET");
            map_ = tensor.raw_type()->getMap()->getLocalMap();
        }

    private:
        DeviceViewType view_;
        LocalMapType   map_;
    };

    template<>
    class DeviceView<TpetraMatrix, 2> {
    public:
        using Scalar          = typename Traits<TpetraMatrix>::Scalar;
        using SizeType        = typename Traits<TpetraMatrix>::SizeType;
        using CrsMatType      = typename TpetraMatrix::CrsMatrixType;
        using LocalMatrixType = typename CrsMatType::local_matrix_type;
        using LocalMapType    = typename CrsMatType::map_type::local_map_type;

        /**
         * Atomic with respect to the value, no the map or idx
         */
        UTOPIA_INLINE_FUNCTION void atomic_add(
            const SizeType &i,
            const SizeType &j,
            const Scalar &value) const
        {
            auto local_i = row_map_.getLocalElement(i);
            auto local_j = col_map_.getLocalElement(j);
            auto row = view_.row(local_i);
            //We assume that rows are sorted
            auto offset = KokkosSparse::findRelOffset(&row.colidx(0), row.length, local_j, 0, true);
            device::atomic_add(&row.value(offset), value);
        }

        DeviceView(TpetraMatrix &tensor) :
        view_(tensor.raw_type()->getLocalMatrix()),
        row_map_(tensor.raw_type()->getRowMap()->getLocalMap()),
        col_map_(tensor.raw_type()->getColMap()->getLocalMap())
        {}

    private:
        LocalMatrixType view_;
        LocalMapType row_map_;
        LocalMapType col_map_;
    };

}

#endif //UTOPIA_TRILINOS_DEVICEVIEW_HPP