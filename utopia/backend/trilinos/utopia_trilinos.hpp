#ifndef UTOPIA_TRILINOS_HPP
#define UTOPIA_TRILINOS_HPP

#include "utopia_Base.hpp"

#ifdef WITH_TRILINOS
// #include <Teuchos_DefaultMpiComm.hpp>
// #include <Tpetra_CrsMatrix.hpp>
// #include <Tpetra_Map.hpp>
// #include <Tpetra_MultiVector.hpp>
// #include <Tpetra_Vector.hpp>
// #include <Tpetra_Version.hpp>
// #include <Teuchos_GlobalMPISession.hpp>
// #include <Teuchos_oblackholestream.hpp>
// #include <Teuchos_Array.hpp>
// #include <Teuchos_ScalarTraits.hpp>
// #include <Teuchos_RCP.hpp>
// #include <Kokkos_Core.hpp>

#include "utopia_trilinos_ForwardDeclaration.hpp"
#include "utopia_trilinos_Traits.hpp"
#include "utopia_trilinos_Types.hpp"
#include "utopia_trilinos_RowView.hpp"
#include "utopia_trilinos_solvers.hpp"
#include "utopia_trilinos_Each.hpp"


//FIXME re-introduce later
// #include "utopia_trilinos_Eval_Factory.hpp"
// #include "utopia_trilinos_Eval_RAP.hpp"
// #include "utopia_trilinos_Eval_Distance.hpp"


#endif //WITH_TRILINOS
#endif //UTOPIA_TRILINOS_HPP
