#include "utopia_trilinos.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Tpetra_Matrix.hpp"
#include "utopia_Tpetra_Vector.hpp"
#include "utopia_trilinos_RowView.hpp"

namespace utopia{
	template class Wrapper<utopia::TpetraMatrix, 2>;
	template class Wrapper<utopia::TpetraVector, 1>;
	template class RowView<utopia::TSMatrixd>;
}
