#ifndef UTOPIA_TRILINOS_VIEW_HOST_HPP
#define UTOPIA_TRILINOS_VIEW_HOST_HPP

#include "utopia_Algorithms.hpp"
#include "utopia_Readable.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ViewHost.hpp"


#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_Vector.hpp"

#include <Trilinos_version.h>

#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
#include <Tpetra_Access.hpp>
#endif

namespace utopia {

    template <>
    class ViewHost<const TpetraVector, 1> {
    public:
        using Scalar = typename Traits<TpetraVector>::Scalar;
        using SizeType = typename Traits<TpetraVector>::SizeType;
        using ExecutionSpaceT = typename TpetraVector::ExecutionSpace;
        using DualViewType = typename TpetraVector::VectorType::dual_view_type;
        using HostExecutionSpaceT = typename DualViewType::host_mirror_space;
        using ViewHostType = typename DualViewType::t_host;
        using LocalMapType = typename TpetraVector::MapType::local_map_type;

        UTOPIA_INLINE_FUNCTION Scalar get(const SizeType &idx) const {
            auto local_idx = map_.getLocalElement(idx);
            return view_(local_idx, 0);
        }

        ViewHost(const TpetraVector &tensor)
            :
#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
              view_(tensor.raw_type()->template getLocalView<HostExecutionSpaceT>(Tpetra::Access::ReadWrite))
#else
              view_(tensor.raw_type()->template getLocalView<HostExecutionSpaceT>())
#endif
        {
            assert(!tensor.has_ghosts() && "GHOST HANDLING NOT IMPLEMENTED YET");
            map_ = tensor.raw_type()->getMap()->getLocalMap();
        }

    private:
        ViewHostType view_;
        LocalMapType map_;
    };

    template <>
    class ViewHost<TpetraVector, 1> {
    public:
        using Scalar = typename Traits<TpetraVector>::Scalar;
        using SizeType = typename Traits<TpetraVector>::SizeType;
        using ExecutionSpaceT = typename TpetraVector::ExecutionSpace;
        using DualViewType = typename TpetraVector::VectorType::dual_view_type;
        using HostExecutionSpaceT = typename DualViewType::host_mirror_space;
        using ViewHostType = typename DualViewType::t_host;
        using LocalMapType = typename TpetraVector::MapType::local_map_type;

        UTOPIA_INLINE_FUNCTION Scalar get(const SizeType &idx) const {
            auto local_idx = map_.getLocalElement(idx);
            return view_(local_idx, 0);
        }

        UTOPIA_INLINE_FUNCTION void set(const SizeType &idx, const Scalar &value) const {
            auto local_idx = map_.getLocalElement(idx);
            view_(local_idx, 0) = value;
        }

        UTOPIA_INLINE_FUNCTION void add(const SizeType &idx, const Scalar &value) const {
            auto local_idx = map_.getLocalElement(idx);
            view_(local_idx, 0) += value;
        }

        /**
         * Atomic with respect to the value, no the map or idx
         */
        UTOPIA_INLINE_FUNCTION void atomic_add(const SizeType &idx, const Scalar &value) const {
            auto local_idx = map_.getLocalElement(idx);
            device::atomic_add(&view_(local_idx, 0), value);
        }

        // UTOPIA_INLINE_FUNCTION void atomic_set(const SizeType &idx, const Scalar &value) const {
        //     auto local_idx = map_.getLocalElement(idx);
        //     device::atomic_set(&view_(local_idx, 0), value);
        // }

        ViewHost(const TpetraVector &tensor)
            :
#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
              view_(tensor.raw_type()->template getLocalView<HostExecutionSpaceT>(Tpetra::Access::ReadWrite))
#else
              view_(tensor.raw_type()->template getLocalView<HostExecutionSpaceT>())
#endif
        {
            assert(!tensor.has_ghosts() && "GHOST HANDLING NOT IMPLEMENTED YET");
            map_ = tensor.raw_type()->getMap()->getLocalMap();
        }

    private:
        ViewHostType view_;
        LocalMapType map_;
    };

    ////////////////////////////////////////////////////////////////////////

    template <>
    class LocalViewHost<const TpetraVector, 1> {
    public:
        using Scalar = typename Traits<TpetraVector>::Scalar;
        using SizeType = typename Traits<TpetraVector>::SizeType;
        using ExecutionSpaceT = typename TpetraVector::ExecutionSpace;
        
        using DualViewType = typename TpetraVector::VectorType::dual_view_type;
        using HostExecutionSpaceT = typename DualViewType::host_mirror_space;
        using ViewHostType = typename DualViewType::t_host;
        using LocalMapType = typename TpetraVector::MapType::local_map_type;

        UTOPIA_INLINE_FUNCTION Scalar get(const SizeType &idx) const { return view_(idx, 0); }

        LocalViewHost(const TpetraVector &tensor)
            :
#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
              view_(tensor.raw_type()->template getLocalView<HostExecutionSpaceT>(Tpetra::Access::ReadWrite))
#else
              view_(tensor.raw_type()->template getLocalView<HostExecutionSpaceT>())
#endif
        {
        }

        inline const ViewHostType &raw_type() { return view_; }

    private:
        const ViewHostType view_;
    };

    template <>
    class LocalViewHost<TpetraVector, 1> {
    public:
        using Scalar = typename Traits<TpetraVector>::Scalar;
        using SizeType = typename Traits<TpetraVector>::SizeType;
        using ExecutionSpaceT = typename TpetraVector::ExecutionSpace;
        using DualViewType = typename TpetraVector::VectorType::dual_view_type;
        using HostExecutionSpaceT = typename DualViewType::host_mirror_space;
        using ViewHostType = typename DualViewType::t_host;
        using LocalMapType = typename TpetraVector::MapType::local_map_type;

        UTOPIA_INLINE_FUNCTION Scalar get(const SizeType &idx) const { return view_(idx, 0); }

        UTOPIA_INLINE_FUNCTION void set(const SizeType &idx, const Scalar &val) const { view_(idx, 0) = val; }

        LocalViewHost(TpetraVector &tensor)
            :
#if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
              view_(tensor.raw_type()->template getLocalView<HostExecutionSpaceT>(Tpetra::Access::ReadWrite))
#else
              view_(tensor.raw_type()->template getLocalView<HostExecutionSpaceT>())
#endif
        {
        }

        inline ViewHostType &raw_type() { return view_; }

    private:
        ViewHostType view_;
    };

    ///////////////////////////////////////////////////////////////////////

//     template <>
//     class ViewHost<TpetraMatrix, 2> {
//     public:
//         using Scalar = typename Traits<TpetraMatrix>::Scalar;
//         using SizeType = typename Traits<TpetraMatrix>::SizeType;
//         using CrsMatType = typename TpetraMatrix::CrsMatrixType;
// #if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
//         using LocalMatrixType = typename CrsMatType::local_matrix_host_type;
// #else
//         using LocalMatrixType = typename CrsMatType::local_matrix_type;
// #endif
//         using LocalMapType = typename CrsMatType::map_type::local_map_type;

//         /**
//          * Atomic with respect to the value, no the map or idx
//          */
//         UTOPIA_INLINE_FUNCTION void atomic_add(const SizeType &i, const SizeType &j, const Scalar &value) const {
//             auto local_i = row_map_.getLocalElement(i);
//             auto local_j = col_map_.getLocalElement(j);
//             auto row = view_.row(local_i);
//             // We assume that rows are sorted
//             auto offset = KokkosSparse::findRelOffset(&row.colidx(0), row.length, local_j, 0, true);
//             device::atomic_add(&row.value(offset), value);
//         }

//         ViewHost(TpetraMatrix &tensor)
//             :
// #if (TRILINOS_MAJOR_MINOR_VERSION >= 130100 && UTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE)
//               view_(tensor.raw_type()->getLocalMatrixHost()),
// #else
//               view_(tensor.raw_type()->getLocalMatrix()),
// #endif
//               row_map_(tensor.raw_type()->getRowMap()->getLocalMap()),
//               col_map_(tensor.raw_type()->getColMap()->getLocalMap()) {
//         }

//     private:
//         LocalMatrixType view_;
//         LocalMapType row_map_;
//         LocalMapType col_map_;
//     };

}  // namespace utopia

#endif  // UTOPIA_TRILINOS_VIEW_HOST_HPP
