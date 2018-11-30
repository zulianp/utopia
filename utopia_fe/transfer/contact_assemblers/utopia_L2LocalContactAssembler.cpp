#include "utopia_L2LocalContactAssembler.hpp"

namespace utopia {

	bool L2LocalContactAssembler::assemble(
		const Elem &master,
		const int master_side,
		FEType master_type,
		const Elem &slave,
		const int slave_side,
		FEType slave_type,
		L2LocalContactAssembler::Result &result
		)
	{
		return false;
	}

	void L2LocalContactAssembler::print_stats(std::ostream &os) const
	{

	}

}
