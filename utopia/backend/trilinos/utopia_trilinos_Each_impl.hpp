#ifndef UTOPIA_TRILINOS_EACH_IMPL_HPP
#define UTOPIA_TRILINOS_EACH_IMPL_HPP

#include "utopia_trilinos_Each.hpp"
#include "utopia_For.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Tpetra_Vector.hpp"
#include "utopia_Tpetra_Matrix.hpp"

namespace utopia {

    template<class Fun>
    void TpetraMatrixEach::apply(TpetraMatrix &mat, Fun fun)
    {
        auto r = row_range(mat);

        if(r.empty()) {
            return;
        }

        auto impl = raw_type(mat);
        auto local_mat = impl->getLocalMatrix();
        auto n = local_mat.numRows();

        for(decltype(n) i = 0; i < n; ++i) {
            auto row = local_mat.row(i);
            auto n_values = row.length;

            for(decltype(n_values) k = 0; k < n_values; ++k) {
                auto &val = row.value(k);
                val = fun(val);
            }
        }
    }

    template<class Fun>
    void TpetraMatrixEach::apply_transform(TpetraMatrix &mat, Fun fun)
    {
        auto rr = row_range(mat);
        if(rr.empty()) return;

        auto impl = raw_type(mat);
        auto col_map   = impl->getColMap()->getLocalMap();
        auto row_map   = impl->getRowMap()->getLocalMap();
        auto local_mat = impl->getLocalMatrix();

        auto n = local_mat.numRows();

        for(decltype(n) i = 0; i < n; ++i) {
            auto row = local_mat.row(i);
            auto n_values = row.length;

            for(decltype(n_values) k = 0; k < n_values; ++k) {
                auto &val = row.value(k);
                const auto global_row = row_map.getGlobalElement(i);
                const auto global_col = col_map.getGlobalElement(row.colidx(k));

                val = fun(global_row, global_col, val);
            }
        }
    }

    template<class Fun>
    void TpetraMatrixEach::apply_read(const TpetraMatrix &mat, Fun fun)
    {
        auto rr = row_range(mat);
        if(rr.empty()) return;

        auto impl = raw_type(mat);
        auto col_map   = impl->getColMap()->getLocalMap();
        auto row_map   = impl->getRowMap()->getLocalMap();
        auto local_mat = impl->getLocalMatrix();

        auto n = local_mat.numRows();

        for(decltype(n) i = 0; i < n; ++i) {
            auto row = local_mat.row(i);
            auto n_values = row.length;

            for(decltype(n_values) k = 0; k < n_values; ++k) {
                const auto &val = row.value(k);
                const auto global_row = row_map.getGlobalElement(i);
                const auto global_col = col_map.getGlobalElement(row.colidx(k));

                fun(global_row, global_col, val);
            }
        }
    }

    template<class Fun>
    void TpetraVectorEach::apply_read(const TpetraVector &v, Fun fun)
    {
        auto impl = raw_type(v);
        auto view = impl->getLocalView<Kokkos::HostSpace>();
        auto map  = impl->getMap()->getLocalMap();

        const auto r = range(v);

        For<>::apply(
            0,
            r.extent(),
            [&map, &view, &fun](const std::size_t i) {
                fun(map.getGlobalElement(i), view(i, 0));
            }
        );
    }

    template<class Fun>
    void TpetraVectorEach::apply_write(TpetraVector &v, Fun fun)
    {
        auto impl = raw_type(v);
        auto view = impl->getLocalView<Kokkos::HostSpace>();
        auto map  = impl->getMap()->getLocalMap();

        const auto r = range(v);

        For<>::apply(
            0,
            r.extent(),
            [&map, &view, &fun](const std::size_t i) {
                view(i, 0) = fun(map.getGlobalElement(i));
            }
        );
    }

    template<class Fun>
    void TpetraVectorEach::apply_transform(const TpetraVector &in, TpetraVector &out, Fun fun)
    {
        const auto r = range(out);


        if(&in == &out) {
            auto impl = raw_type(out);
            auto view = impl->getLocalView<Kokkos::HostSpace>();
            auto map  = impl->getMap()->getLocalMap();

            For<>::apply(
                0,
                r.extent(),
                [&map, &view, &fun](const std::size_t i) {
                    auto &val = view(i, 0);
                    val = fun(map.getGlobalElement(i), val);
                }
            );

        } else {

            auto impl_in = raw_type(in);
            auto view_in = impl_in->getLocalView<Kokkos::HostSpace>();
            auto map_in  = impl_in->getMap()->getLocalMap();


            auto impl_out = raw_type(out);
            auto view_out = impl_out->getLocalView<Kokkos::HostSpace>();
            auto map_out  = impl_out->getMap()->getLocalMap();

            For<>::apply(
                0,
                r.extent(),
                [&view_in, &map_in, &view_out, map_out, &fun](const std::size_t i) {;
                    assert(map_in.getGlobalElement(i) == map_out.getGlobalElement(i));

                    const auto &val = view_in(i, 0);
                    view_out(i, 0) = fun(map_in.getGlobalElement(i), val);
                }
            );
        }
    }

}

#endif //UTOPIA_TRILINOS_EACH_IMPL_HPP
