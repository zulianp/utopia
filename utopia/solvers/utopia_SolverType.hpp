#ifndef UTOPIA_SOLVER_TYPE_HPP
#define UTOPIA_SOLVER_TYPE_HPP

#include "utopia_Base.hpp"

namespace utopia {
    typedef const char * SolverType;
    typedef const char * SolverPackage;

    // class SolverPackage {
    // public:
    // 	constexpr SolverPackage(const char * value)
    // 	: value(value) {}

    // 	constexpr operator const char *() const
    // 	{
    // 		return value;
    // 	}

    // 	const char * value;
    // };

    class Solver {
    public:
        //all solvers
        inline constexpr static SolverType automatic() { return "auto"; }

        //linear solvers
        inline constexpr static SolverType utopia_cg() { return "utopia_cg"; }
        inline constexpr static SolverType cg() 	   { return "cg"; }
        inline constexpr static SolverType bicgstab()  { return "bicgstab"; }
        inline constexpr static SolverType direct()    { return "direct"; }
        inline constexpr static SolverType ksp()       { return "ksp"; }

        inline constexpr static SolverType lu_decomposition() 		{ return "lu_decomposition"; }
        inline constexpr static SolverType cholesky_decomposition() { return "cholesky_decomposition"; }

// #ifdef WITH_UMFPACK
        inline constexpr static SolverPackage umfpack()   { return "umfpack"; }
// #endif //WITH_UMFPACK

// #ifdef PETSC_HAVE_MUMPS
        inline constexpr static SolverPackage mumps()   { return "mumps"; }
// #endif //PETSC_HAVE_MUMPS

// #ifdef WITH_SUPERLU_DIST
        inline constexpr static SolverType superlu_dist() { return "superlu_dist"; }
// #endif //WITH_SUPERLU_DIST

        inline constexpr static SolverPackage petsc() { return "petsc"; }

        //nonlinear solvers
        static constexpr SolverType newton() 	   { return "newton"; }
        static constexpr SolverType line_search()  { return "line_search"; }
        static constexpr SolverType trust_region() { return "trust_region"; }

        //trust-region solvers
        inline constexpr static SolverType cauchypoint()    { return "cauchypoint"; }
        inline constexpr static SolverType dogleg()         { return "dogleg"; }
        inline constexpr static SolverType steihaug_toint() { return "steihaug_toint"; }
        inline constexpr static SolverType toint()          { return "toint"; }
        inline constexpr static SolverType nash()           { return "nash"; }
        inline constexpr static SolverType lanczos()        { return "lanczos"; }
        inline constexpr static SolverType cgne()           { return "cgne"; }
        inline constexpr static SolverType stcg()           { return "stcg"; }

        //backtracking
        inline constexpr static SolverType simple_backtracking() { return "simple_backtracking"; }
        inline constexpr static SolverType backtracking() 		 { return "backtracking"; }
    };
}


#endif //UTOPIA_SOLVER_TYPE_HPP
