#ifndef UTOPIA_TRILINOS_EACH_IMPL_HPP
#define UTOPIA_TRILINOS_EACH_IMPL_HPP

#include "utopia_For.hpp"
#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_Vector.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_trilinos_Each.hpp"

#include <Trilinos_version.h>

#if TRILINOS_MAJOR_MINOR_VERSION >= 130100
#include <Tpetra_Access.hpp>
#endif

namespace utopia {

    template <class Fun>
    void TpetraMatrixEach::apply(TpetraMatrix &mat, Fun fun) {
        auto r = row_range(mat);

        if (r.empty()) {
            return;
        }

        auto impl = raw_type(mat);

#if TRILINOS_MAJOR_MINOR_VERSION >= 130100
        auto local_mat = impl->getLocalMatrixDevice();
#else
        auto local_mat = impl->getLocalMatrix();
#endif

        auto n = local_mat.numRows();

        for (decltype(n) i = 0; i < n; ++i) {
            auto row = local_mat.row(i);
            auto n_values = row.length;

            for (decltype(n_values) k = 0; k < n_values; ++k) {
                auto &val = row.value(k);
                val = fun(val);
            }
        }
    }

    template <class Fun>
    void TpetraMatrixEach::apply_transform(TpetraMatrix &mat, Fun fun) {
        auto rr = row_range(mat);
        if (rr.empty()) return;

        auto impl = raw_type(mat);
        auto col_map = impl->getColMap()->getLocalMap();
        auto row_map = impl->getRowMap()->getLocalMap();

#if TRILINOS_MAJOR_MINOR_VERSION >= 130100
        auto local_mat = impl->getLocalMatrixDevice();
#else
        auto local_mat = impl->getLocalMatrix();
#endif

        auto n = local_mat.numRows();

        for (decltype(n) i = 0; i < n; ++i) {
            auto row = local_mat.row(i);
            auto n_values = row.length;

            for (decltype(n_values) k = 0; k < n_values; ++k) {
                auto &val = row.value(k);
                const auto global_row = row_map.getGlobalElement(i);
                const auto global_col = col_map.getGlobalElement(row.colidx(k));

                val = fun(global_row, global_col, val);
            }
        }
    }

    template <class Fun>
    void TpetraMatrixEach::apply_read(const TpetraMatrix &mat, Fun fun) {
        auto rr = row_range(mat);
        if (rr.empty()) return;

        auto impl = raw_type(mat);
        auto col_map = impl->getColMap()->getLocalMap();
        auto row_map = impl->getRowMap()->getLocalMap();

#if TRILINOS_MAJOR_MINOR_VERSION >= 130100
        auto local_mat = impl->getLocalMatrixDevice();
#else
        auto local_mat = impl->getLocalMatrix();
#endif

        auto n = local_mat.numRows();

        for (decltype(n) i = 0; i < n; ++i) {
            auto row = local_mat.row(i);
            auto n_values = row.length;

            for (decltype(n_values) k = 0; k < n_values; ++k) {
                const auto &val = row.value(k);
                const auto global_row = row_map.getGlobalElement(i);
                const auto global_col = col_map.getGlobalElement(row.colidx(k));

                fun(global_row, global_col, val);
            }
        }
    }

    template <class Fun>
    void TpetraVectorEach::apply_read(const TpetraVector &v, Fun fun) {
        using ExecutionSpaceT = TpetraVector::ExecutionSpace;

        auto impl = raw_type(v);

        auto map = impl->getMap()->getLocalMap();

#if TRILINOS_MAJOR_MINOR_VERSION >= 130100

        auto view = impl->template getLocalView<ExecutionSpaceT>(Tpetra::Access::ReadOnly);
#else
        auto view = impl->getLocalView<ExecutionSpaceT>();
#endif

        const auto r = range(v);

        For<>::apply(
            0, r.extent(), [&map, &view, &fun](const std::size_t i) { fun(map.getGlobalElement(i), view(i, 0)); });
    }

    template <class Fun>
    void TpetraVectorEach::apply_write(TpetraVector &v, Fun fun) {
        using ExecutionSpaceT = TpetraVector::ExecutionSpace;

        auto impl = raw_type(v);
        auto map = impl->getMap()->getLocalMap();

#if TRILINOS_MAJOR_MINOR_VERSION >= 130100
        auto view = impl->template getLocalView<ExecutionSpaceT>(Tpetra::Access::OverwriteAll);
#else
        auto view = impl->getLocalView<ExecutionSpaceT>();
#endif

        const auto r = range(v);

        For<>::apply(
            0, r.extent(), [&map, &view, &fun](const std::size_t i) { view(i, 0) = fun(map.getGlobalElement(i)); });
    }

    template <class Fun>
    void TpetraVectorEach::apply_transform(const TpetraVector &in, TpetraVector &out, Fun fun) {
        const auto r = range(in);

        if (out.empty() || in.size() != out.size()) {
            // make copy
            out = in;
        }

        assert(r == range(out));

        if (in.is_alias(out)) {
            auto impl = raw_type(out);
#if TRILINOS_MAJOR_MINOR_VERSION >= 130100
            auto view = impl->getLocalViewHost(Tpetra::Access::ReadWrite);
#else
            auto view = impl->getLocalViewHost();
#endif
            auto map = impl->getMap()->getLocalMap();

            For<>::apply(0, r.extent(), [&map, &view, &fun](const std::size_t i) {
                auto &val = view(i, 0);
                val = fun(map.getGlobalElement(i), val);
            });

        } else {
            auto impl_in = raw_type(in);

            auto map_in = impl_in->getMap()->getLocalMap();

#if TRILINOS_MAJOR_MINOR_VERSION >= 130100
            auto view_in = impl_in->getLocalViewHost(Tpetra::Access::ReadWrite);
#else
            auto view_in = impl_in->getLocalViewHost();
#endif

            auto impl_out = raw_type(out);
            auto map_out = impl_out->getMap()->getLocalMap();

#if TRILINOS_MAJOR_MINOR_VERSION >= 130100
            auto view_out = impl_out->getLocalViewHost(Tpetra::Access::ReadWrite);
#else
            auto view_out = impl_out->getLocalViewHost();
#endif

            For<>::apply(0, r.extent(), [&view_in, &map_in, &view_out, map_out, &fun](const std::size_t i) {
                ;
                assert(map_in.getGlobalElement(i) == map_out.getGlobalElement(i));

                const auto &val = view_in(i, 0);
                view_out(i, 0) = fun(map_in.getGlobalElement(i), val);
            });
        }
    }

}  // namespace utopia

#endif  // UTOPIA_TRILINOS_EACH_IMPL_HPP
