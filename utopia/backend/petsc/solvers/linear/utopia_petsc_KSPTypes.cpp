#include "utopia_petsc_KSPTypes.hpp"
#include "utopia_petsc_Base.hpp"
#include "utopia_petsc_Types.hpp"

#include <petscksp.h>
#include <petscpc.h>

namespace utopia {
    // TODO(zulianp): check what solvers are available in petsc based on the compilation and version

    KSPTypes::KSPTypes()
        : ksp_({
              KSPRICHARDSON, KSPCHEBYSHEV, KSPCG, KSPGROPPCG, KSPPIPECG, KSPPIPECGRR,
                  // Not reallu there in 3.9.0
                  // #if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3,9,0)
                  // 		KSPPIPELCG,
                  // #endif
                  KSPCGNE,
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 8, 0)
                  KSPFETIDP, KSPPIPEBCGS,
#endif
        // Seriously??
#if (UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 8, 0) && UTOPIA_PETSC_VERSION_LESS_THAN(3, 12, 0))
                  KSPCGNASH, KSPCGSTCG, KSPCGGLTR,
#else
                  KSPSTCG, KSPGLTR, KSPNASH,
#endif
                  KSPFCG, KSPPIPEFCG, KSPGMRES, KSPPIPEFGMRES, KSPFGMRES, KSPLGMRES, KSPDGMRES, KSPPGMRES, KSPTCQMR,
                  KSPBCGS, KSPIBCGS, KSPFBCGS, KSPFBCGSR, KSPBCGSL, KSPCGS, KSPTFQMR, KSPCR, KSPPIPECR, KSPLSQR,
                  KSPPREONLY, KSPQCG, KSPBICG, KSPMINRES, KSPSYMMLQ, KSPLCD, KSPPYTHON, KSPGCR, KSPPIPEGCR, KSPTSIRM,
                  KSPCGLS
          }),
          pc_({
              PCJACOBI, PCSOR, PCLU, PCSHELL, PCBJACOBI, PCMG, PCEISENSTAT, PCILU, PCICC, PCASM, PCGASM, PCKSP,
                  PCCOMPOSITE, PCREDUNDANT, PCSPAI, PCNN, PCCHOLESKY, PCPBJACOBI, PCMAT, PCHYPRE, PCPARMS, PCFIELDSPLIT,
                  PCTFS, PCML, PCGALERKIN, PCEXOTIC, PCCP, PCBFBT, PCLSC, PCPYTHON, PCPFMG, PCSYSPFMG, PCREDISTRIBUTE,
                  PCSVD, PCGAMG,
#if UTOPIA_PETSC_VERSION_GREATER_EQUAL_THAN(3, 8, 0)
                  PCCHOWILUVIENNACL, PCROWSCALINGVIENNACL, PCSAVIENNACL,
#endif
                  PCBDDC, PCKACZMARZ, PCTELESCOPE, PCNONE
          }),
          package_({
#ifdef PETSC_HAVE_SUPERLU
              MATSOLVERSUPERLU,
#endif
#ifdef PETSC_HAVE_SUPERLU_DIST
              MATSOLVERSUPERLU_DIST,
#endif
              // MATSOLVERSTRUMPACK,
              // MATSOLVERUMFPACK,
              // MATSOLVERCHOLMOD,
              // MATSOLVERKLU,
              // MATSOLVERSPARSEELEMENTAL,
              // MATSOLVERELEMENTAL,
              // MATSOLVERESSL,
              MATSOLVERLUSOL,
#ifdef PETSC_HAVE_MUMPS
              MATSOLVERMUMPS,
#endif
              // MATSOLVERMKL_PARDISO,
              // MATSOLVERMKL_CPARDISO,
              // MATSOLVERPASTIX,
              // MATSOLVERMATLAB,
              MATSOLVERPETSC,
              // MATSOLVERBAS,
              // MATSOLVERCUSPARSE,
              " "}) {
    }
}  // namespace utopia
