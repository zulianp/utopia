#ifndef UTOPIA_KOKKOS_PARALLEL_EACH_HPP
#define UTOPIA_KOKKOS_PARALLEL_EACH_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Base.hpp"
#include "utopia_Each.hpp"
#include "utopia_trilinos_Types.hpp"
#include "utopia_ParallelEach.hpp"
#include "utopia_For.hpp"


#include <Kokkos_Core.hpp>

#include <memory>

namespace utopia {

    // template<class Tensor, int Order = Tensor::Order, int FILL_TYPE = Tensor::FILL_TYPE>
    // class ParallelEach {};

    template<int FILL_TYPE>
    class ParallelEach<TpetraVector, 1, FILL_TYPE>{
    public:
        template<class Fun>
        inline static void apply_write(const TpetraVector &v, Fun fun, const std::string &name)
        {
            using ExecutionSpaceT = TpetraVector::vector_type::execution_space;

            auto k_v = raw_type(v)->getLocalView<ExecutionSpaceT>();
            auto offset = range(v).begin();
            Kokkos::parallel_for(
                name,
                k_v.extent(0),
                KOKKOS_LAMBDA(const int i) {
                    k_v(i, 0) = fun(offset + i);
            });
        }

        template<class Fun>
        inline static void apply_read(const TpetraVector &v, Fun fun, const std::string &name)
        {
            using ExecutionSpaceT = TpetraVector::vector_type::execution_space;

            auto k_v = raw_type(v)->getLocalView<ExecutionSpaceT>();
            auto offset = range(v).begin();
            Kokkos::parallel_for(
                name,
                k_v.extent(0),
                KOKKOS_LAMBDA(const int i) {
                    fun(offset + i, k_v(i, 0));
            });
        }


        template<class Fun>
        inline static void apply_transform(const TpetraVector &v, Fun fun, const std::string &name)
        {
            using ExecutionSpaceT = TpetraVector::vector_type::execution_space;

            auto k_v = raw_type(v)->getLocalView<ExecutionSpaceT>();
            auto offset = range(v).begin();
            Kokkos::parallel_for(
                name,
                k_v.extent(0),
                KOKKOS_LAMBDA(const int i) {
                   k_v(i, 0) =  fun(offset + i, k_v(i, 0));
            });
        }
    };


    template<int FILL_TYPE>
    class ParallelEach<TpetraMatrixd, 2, FILL_TYPE>{
    public:
        template<class Fun>
        inline static void apply_write(TpetraMatrixd &mat, Fun fun, const std::string &name)
        {
             typedef Kokkos::TeamPolicy<>               team_policy;
             typedef Kokkos::TeamPolicy<>::member_type  member_type;

             auto r = row_range(mat);

             if(r.empty()) {
                 return;
             }

             auto impl = raw_type(mat);
             auto local_mat = impl->getLocalMatrix();

             auto n = local_mat.numRows();

             auto row_map = impl->getRowMap()->getLocalMap();
             auto col_map = impl->getColMap()->getLocalMap();

             Kokkos::parallel_for(name, team_policy(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member_type &team_member) {
                 const int row_ind = team_member.league_rank();
                 auto row = local_mat.row(row_ind);
                 auto n_values = row.length;

                 Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_values), [&] (const int k) {
                     auto &val    = row.value(k);
                     auto col_ind = row.colidx(k);
                     val = fun(
                        row_map.getGlobalElement(row_ind),
                        col_map.getGlobalElement(col_ind)
                        );
                 });
             });
        }

        template<class Fun>
        inline static void apply_read(const TpetraMatrixd &mat, Fun fun, const std::string &name)
        {
            typedef Kokkos::TeamPolicy<>               team_policy;
            typedef Kokkos::TeamPolicy<>::member_type  member_type;

            auto r = row_range(mat);

            if(r.empty()) {
                return;
            }

            auto impl = raw_type(mat);
            auto local_mat = impl->getLocalMatrix();

            auto n = local_mat.numRows();

            auto row_map = impl->getRowMap()->getLocalMap();
            auto col_map = impl->getColMap()->getLocalMap();

            Kokkos::parallel_for(name, team_policy(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member_type &team_member) {
                const int row_ind = team_member.league_rank();
                auto row = local_mat.row(row_ind);
                auto n_values = row.length;

                Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_values), [&] (const int k) {
                    auto &val    = row.value(k);
                    auto col_ind = row.colidx(k);
                    fun(
                       row_map.getGlobalElement(row_ind),
                       col_map.getGlobalElement(col_ind),
                       val
                       );
                });
            });
        }

        template<class Fun>
        inline static void apply_transform(TpetraMatrixd &mat, Fun fun)
        {
            typedef Kokkos::TeamPolicy<>               team_policy;
            typedef Kokkos::TeamPolicy<>::member_type  member_type;

            auto r = row_range(mat);

            if(r.empty()) {
                return;
            }

            auto impl = raw_type(mat);
            auto local_mat = impl->getLocalMatrix();

            auto n = local_mat.numRows();

            auto row_map = impl->getRowMap()->getLocalMap();
            auto col_map = impl->getColMap()->getLocalMap();

            Kokkos::parallel_for(team_policy(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member_type &team_member) {
                const int row_ind = team_member.league_rank();
                auto row = local_mat.row(row_ind);
                auto n_values = row.length;

                Kokkos::parallel_for(Kokkos::TeamThreadRange(team_member, n_values), [&](const int k) {
                    auto &val    = row.value(k);
                    auto col_ind = row.colidx(k);
                    val = fun(
                       row_map.getGlobalElement(row_ind),
                       col_map.getGlobalElement(col_ind),
                       val
                       );
                });
            });
        }
    };

    template<class Fun>
    inline void parallel_transform(TpetraMatrixd &mat, Fun fun)
    {
        typedef Kokkos::TeamPolicy<>               team_policy;
        typedef Kokkos::TeamPolicy<>::member_type  member_type;

        auto r = row_range(mat);

        if(r.empty()) {
            return;
        }

        auto impl = raw_type(mat);
        auto local_mat = impl->getLocalMatrix();

        auto n = local_mat.numRows();

        Kokkos::parallel_for(team_policy(n, Kokkos::AUTO), KOKKOS_LAMBDA(const member_type &team_member) {
            const int j = team_member.league_rank();
            auto row      = local_mat.row(j);
            auto n_values = row.length;

            Kokkos::parallel_for( Kokkos::TeamThreadRange(team_member, n_values), [&] (const int i) {
                auto &val = row.value(i);
                val = fun(val);
            });
        });
    }

    template<>
    class ParallelFor<TRILINOS> {
    public:
        template<typename F>
        inline static void apply(const Range &r, F f)
        {
            apply(r.begin(), r.end(), f);
        }

        template<typename F>
        inline static void apply(
            const std::size_t &begin,
            const std::size_t &end,
            F f)
        {
            auto extent = end - begin;
            Kokkos::parallel_for(extent, KOKKOS_LAMBDA(const int i) {
                f(begin + i);
            });
        }

        template<typename F>
        inline static void apply(
            const std::size_t &n,
            F f)
        {
            Kokkos::parallel_for(n, KOKKOS_LAMBDA(const int i) {
                f(i);
            });
        }
    };

}

#endif //UTOPIA_KOKKOS_PARALLEL_EACH_HPP
