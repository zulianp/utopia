#ifndef UTOPIA_TRILINOS_HPP
#define UTOPIA_TRILINOS_HPP 


#ifdef WITH_TRILINOS
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>

#include "utopia_trilinos_ForwardDeclaration.hpp"
//#include "utopia_trilinos_Vector.hpp"
//#include "utopia_trilinos_Matrix.hpp"
//#include "utopia_trilinos_Backend.hpp"
#include "utopia_trilinos_Traits.hpp"
#include "utopia_trilinos_Types.hpp"


#endif //WITH_TRILINOS
#endif //UTOPIA_TRILINOS_HPP
