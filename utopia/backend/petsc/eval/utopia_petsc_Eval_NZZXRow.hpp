#ifndef UTOPIA_PETSC_EVALN_ZZ_X_ROW_HPP
#define UTOPIA_PETSC_EVALN_ZZ_X_ROW_HPP

namespace utopia {

    // template<class Left, class Index, class Traits>
    // class Eval< Construct<Left, Factory<NNZXRow<Index>, 2>  >, Traits, PETSC> {
    // public:
    //     typedef utopia::Construct<Left, Factory<NNZXRow<Index>, 2>  > Expr;

    //     inline static void apply(const Expr &expr) {
    //         UTOPIA_TRACE_BEGIN(expr);

    //         const auto &r = expr.right();
    //         const auto &d_nnz = r.type().d_nnz;
    //         const auto &o_nnz = r.type().o_nnz;
    //         const auto &gs    = r.size();
    //         auto &A = Eval<Left, Traits>::apply(expr.left());

    //         //FIXME
    //         A.matij_init(
    //             A.comm().get(),
    //             A.type_override(),
    //             d_nnz.size(),
    //             PETSC_DETERMINE,
    //             gs.get(0),
    //             gs.get(1),
    //             d_nnz,
    //             o_nnz
    //             );

    //         UTOPIA_TRACE_END(expr);
    //     }
    // };

    template <class Left, class Index, class Traits>
    class Eval<Assign<Left, Factory<NNZXRow<Index>, 2> >, Traits, PETSC> {
    public:
        typedef utopia::Assign<Left, Factory<NNZXRow<Index>, 2> > Expr;

        inline static void apply(const Expr &expr) {
            UTOPIA_TRACE_BEGIN(expr);

            const auto &r = expr.right();
            const auto &d_nnz = r.type().d_nnz;
            const auto &o_nnz = r.type().o_nnz;
            const auto &gs = r.size();
            auto &A = Eval<Left, Traits>::apply(expr.left());

            // FIXME
            A.matij_init(
                A.comm().get(), A.type_override(), d_nnz.size(), PETSC_DETERMINE, gs.get(0), gs.get(1), d_nnz, o_nnz);

            UTOPIA_TRACE_END(expr);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_PETSC_EVALN_ZZ_X_ROW_HPP
