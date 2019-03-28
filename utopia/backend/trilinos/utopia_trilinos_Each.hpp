#ifndef UTOPIA_TRILINOS_EACH_HPP
#define UTOPIA_TRILINOS_EACH_HPP

#include "utopia_trilinos_Types.hpp"
#include "utopia_trilinos_Traits.hpp"
#include "utopia_For.hpp"

namespace utopia {

    template<>
    class Each<TSMatrixd, 2, FillType::SPARSE>  {
    public:

        using SizeType = typename Traits<TSMatrixd>::SizeType;
        using Scalar   = typename Traits<TSMatrixd>::Scalar;

        template<class Fun>
        inline static void apply(TSMatrixd &mat, Fun fun)
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
        inline static void apply_transform(TSMatrixd &mat, Fun fun)
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

        inline static void apply_read(const TSMatrixd &mat, std::function<void(const Scalar &)> &fun)
        {
            //FIXME make performance version
            apply_read(mat, [&fun](const SizeType i, const SizeType j, const Scalar val) {
                UTOPIA_UNUSED(i);
                UTOPIA_UNUSED(j);

                fun(val);
            });
        }

        template<class Fun>
        inline static void apply_read(const TSMatrixd &mat, Fun fun)
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
    };



    template<int FILL_TYPE>
    class Each<TVectord, 1, FILL_TYPE> {
    public:
        using SizeType = typename Traits<TVectord>::SizeType;
        using Scalar   = typename Traits<TVectord>::Scalar;

        inline static void apply_read(const TVectord &v, std::function<void(const Scalar &)> &fun)
        {
            auto impl = raw_type(v);
            auto view = impl->getLocalView<Kokkos::HostSpace>();

            const auto r = range(v);

            For<>::apply(
                0,
                r.extent(),
                [&view, &fun](const std::size_t i) {
                    fun(view(i, 0));
                }
            );
        }

        template<class Fun>
        inline static void apply_read(const TVectord &v, Fun fun)
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
        inline static void apply_write(TVectord &v, Fun fun)
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
        inline static void apply_transform(const TVectord &in, TVectord &out, Fun fun)
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

    };

    template<class Fun>
    inline void each_transform(TVectord &in_out, Fun fun)
    {
        each_transform(in_out, in_out, fun);
    }

}

#endif //UTOPIA_TRILINOS_EACH_HPP
