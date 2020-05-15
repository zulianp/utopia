#ifndef UTOPIA_LIBMESH_HPP
#define UTOPIA_LIBMESH_HPP

#include "utopia_fe_kokkos_fix.hpp"

#include "utopia_Base.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_libmesh_AssemblyContext.hpp"
#include "utopia_libmesh_FEBackend.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_FormEval.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_FunctionalTraits.hpp"
#include "utopia_libmesh_LambdaFunction.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"

#include "utopia_FEFunction.hpp"

#include "libmesh/parallel_mesh.h"

#ifdef WITH_TINY_EXPR
#include "utopia_SymbolicFunction.hpp"
namespace utopia {

    inline ContextFunction<std::vector<double>,
                           std::function<std::vector<double>(const AssemblyContext<LIBMESH_TAG> &)> >
    symbolic_to_ctx_fun(const std::string &expr) {
        std::function<std::vector<double>(const AssemblyContext<LIBMESH_TAG> &)> f =
            [expr](const AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<double> {
            SymbolicFunction fun(expr);
            const auto &pts = ctx.fe()[0]->get_xyz();

            const auto n = pts.size();
            std::vector<double> ret(n);

            for (std::size_t i = 0; i != n; ++i) {
                ret[i] = fun.eval(pts[i](0), pts[i](1), pts[i](2));
            }

            return ret;
        };

        return ctx_fun<std::vector<double> >(f);
    }
}  // namespace utopia

#endif  // WITH_TINY_EXPR

#endif  // UTOPIA_LIBMESH_HPP
