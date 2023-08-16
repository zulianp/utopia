#include "utopia_Base.hpp"
#include "utopia_ConjugateGradient.hpp"  // FIXME(zulianp)
#include "utopia_InputParameters.hpp"
#include "utopia_polymorphic_LinearSolver.hpp"

#ifdef UTOPIA_ENABLE_PETSC
#include "utopia_petsc.hpp"
#include "utopia_petsc_BDDLinearSolver.hpp"
#include "utopia_petsc_LinearSolverFactory.hpp"

using Matrix_t = utopia::PetscMatrix;
using Solver_t = utopia::KSPSolver<utopia::PetscMatrix, utopia::PetscVector>;
// using Solver_t = utopia::BDDLinearSolver<utopia::PetscMatrix, utopia::PetscVector>;
#else
#ifdef UTOPIA_ENABLE_TRILINOS
#include "utopia_trilinos_Types.hpp"
using Matrix_t = utopia::TpetraMatrixd;
using Solver_t = utopia::ConjugateGradient<utopia::TpetraMatrixd, utopia::TpetraVectord>;
#else

#ifdef UTOPIA_ENABLE_BLAS
#error "Blas backend not enough for this feature!\nPlease use -DUTOPIA_ENABLE_ISOLVER=OFF"
#endif  // UTOPIA_ENABLE_BLAS
#endif  // UTOPIA_ENABLE_TRILINOS
#endif  // UTOPIA_ENABLE_PETSC

using Vector_t = typename utopia::Traits<Matrix_t>::Vector;
using Scalar_t = typename utopia::Traits<Matrix_t>::Scalar;
using Size_t = typename utopia::Traits<Matrix_t>::SizeType;

#define Scalar_t ISOLVER_SCALAR_T;
#define Size_t ISOLVER_IDX_T;

#define UTOPIA_READ_ENV(name, conversion) \
    do {                                  \
        char *var = getenv(#name);        \
        if (var) {                        \
            name = conversion(var);       \
        }                                 \
    } while (0)

extern "C" {
#include "isolver_lsolve.h"

int ISOLVER_EXPORT isolver_lsolve_init(isolver_lsolve_t *info) {
#ifdef UTOPIA_ENABLE_PETSC
    PetscInitializeNoArguments();
#endif  // UTOPIA_ENABLE_PETSC

    auto solver = new Solver_t();
#ifdef UTOPIA_ENABLE_PETSC
    solver->pc_type("hypre");
    solver->ksp_type("gmres");
#endif  // UTOPIA_ENABLE_PETSC

    info->private_data = (void *)solver;

    const char *UTOPIA_LINEAR_SOLVER_CONFIG = 0;
    UTOPIA_READ_ENV(UTOPIA_LINEAR_SOLVER_CONFIG, );

    if (UTOPIA_LINEAR_SOLVER_CONFIG) {
        printf("UTOPIA_LINEAR_SOLVER_CONFIG=%s\n", UTOPIA_LINEAR_SOLVER_CONFIG);
        solver->import(UTOPIA_LINEAR_SOLVER_CONFIG);
    }
    // else {
    //     utopia::param_list(utopia::param("type", "bdd"), utopia::param("num_blocks", 8));
    // }

    return 0;
}

int ISOLVER_EXPORT isolver_lsolve_update_crs(const isolver_lsolve_t *info,
                                             const ptrdiff_t n_local,
                                             const ptrdiff_t n_global,
                                             const isolver_idx_t *const rowptr,
                                             const isolver_idx_t *const colidx,
                                             const isolver_scalar_t *const values) {
    auto mat = std::make_shared<Matrix_t>();

    mat->wrap(info->comm,
              n_local,
              n_local,
              n_global,
              n_global,
              (isolver_idx_t *)rowptr,
              (isolver_idx_t *)colidx,
              (isolver_scalar_t *)values,
              []() {
                  // Resources are freed outside!
              });

    auto solver = (Solver_t *)info->private_data;
    solver->update(mat);
    return 0;
}

int ISOLVER_EXPORT isolver_lsolve_update_with_preconditioner_crs(const isolver_lsolve_t *info,
                                                                 const ptrdiff_t n_local,
                                                                 const ptrdiff_t n_global,
                                                                 const isolver_idx_t *const rowptr,
                                                                 const isolver_idx_t *const colidx,
                                                                 const isolver_scalar_t *const values,
                                                                 const isolver_idx_t *const prec_rowptr,
                                                                 const isolver_idx_t *const prec_colidx,
                                                                 const isolver_scalar_t *const prec_values) {
    auto mat = std::make_shared<Matrix_t>();

    mat->wrap(info->comm,
              n_local,
              n_local,
              n_global,
              n_global,
              (isolver_idx_t *)rowptr,
              (isolver_idx_t *)colidx,
              (isolver_scalar_t *)values,
              []() {
                  // Resources are freed outside!
              });

    auto prec = std::make_shared<Matrix_t>();

    prec->wrap(info->comm,
               n_local,
               n_local,
               n_global,
               n_global,
               (isolver_idx_t *)prec_rowptr,
               (isolver_idx_t *)prec_colidx,
               (isolver_scalar_t *)prec_values,
               []() {
                   // Resources are freed outside!
               });

    auto solver = (Solver_t *)info->private_data;
    solver->update(mat, prec);

    return 0;
}

int ISOLVER_EXPORT isolver_lsolve_apply(const isolver_lsolve_t *info,
                                        const isolver_scalar_t *const rhs,
                                        isolver_scalar_t *const x) {
    auto solver = (Solver_t *)info->private_data;
    auto op = solver->get_operator();

    Vector_t w_rhs, w_x;

    auto ls = op->local_size().get(0);
    auto gs = op->size().get(0);

    w_x.wrap(info->comm, ls, gs, x, []() {
        // Resources are freed outside!
    });

    w_rhs.wrap(info->comm, ls, gs, rhs, []() {
        // Resources are freed outside!
    });

    solver->apply(w_rhs, w_x);

    return 0;
}

int ISOLVER_EXPORT isolver_lsolve_set_max_iterations(const isolver_lsolve_t *info, const int max_it) {
    auto solver = (Solver_t *)info->private_data;

    auto p = utopia::param_list(utopia::param("max_it", max_it));
    // solver->max_it(max_it);
    solver->read(p);
    return 0;
}

int ISOLVER_EXPORT isolver_lsolve_set_atol(const isolver_lsolve_t *info, const isolver_scalar_t atol) {
    auto solver = (Solver_t *)info->private_data;
    // solver->atol(atol);
    auto p = utopia::param_list(utopia::param("atol", atol));
    // solver->max_it(max_it);
    solver->read(p);
    return 0;
}

int ISOLVER_EXPORT isolver_lsolve_set_stol(const isolver_lsolve_t *info, const isolver_scalar_t stol) {
    auto solver = (Solver_t *)info->private_data;
    auto p = utopia::param_list(utopia::param("stol", stol));
    // solver->stol(stol);
    solver->read(p);
    return 0;
}

int ISOLVER_EXPORT isolver_lsolve_set_rtol(const isolver_lsolve_t *info, const isolver_scalar_t rtol) {
    auto solver = (Solver_t *)info->private_data;
    auto p = utopia::param_list(utopia::param("rtol", rtol));
    // solver->rtol(rtol);
    solver->read(p);
    return 0;
}

int ISOLVER_EXPORT isolver_lsolve_set_verbosity(const isolver_lsolve_t *info, const int verbosity_level) {
    auto solver = (Solver_t *)info->private_data;

    auto p = utopia::param_list(utopia::param("verbose", verbosity_level > 0));
    // solver->verbose(verbosity_level);
    solver->read(p);
    return 0;
}

int ISOLVER_EXPORT isolver_lsolve_destroy(isolver_lsolve_t *info) {
    auto solver = (Solver_t *)info->private_data;
    delete solver;
    info->private_data = nullptr;

#ifdef UTOPIA_ENABLE_PETSC
    PetscFinalize();
#endif  // UTOPIA_ENABLE_PETSC
    return 0;
}
}
