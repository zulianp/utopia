
#include "utopia_Main.hpp"

#include "utopia_fe_base.hpp"

#ifdef UTOPIA_WITH_LIBMESH

#include "utopia_ObstacleApp.hpp"

#include "utopia_ui.hpp"

#include "utopia_Agglomerate.hpp"
#include "utopia_BlockAgglomerate.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_ILU.hpp"

#include "utopia_MPITimeStatistics.hpp"
#include "utopia_petsc_AdditiveCorrectionTransfer.hpp"
#include "utopia_petsc_DILUAlgorithm.hpp"

#include "utopia_libmesh_Obstacle.hpp"
#include "utopia_libmesh_OmniAssembler.hpp"

void obs(utopia::Input &in) {
    utopia::ObstacleProblem<utopia::libmesh::FunctionSpace> obs;
    obs.read(in);
    if (obs.valid) {
        obs.run();
    } else {
        utopia::err() << "[Error] invalid set-up\n";
    }
}

UTOPIA_REGISTER_APP(obs);

#endif  // UTOPIA_WITH_LIBMESH
