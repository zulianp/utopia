#include "utopia_MixedBidirectionalLocalAssembler.hpp"
#include "utopia.hpp"

namespace utopia {
	MixedBidirectionalLocalAssembler::	MixedBidirectionalLocalAssembler(
			const std::shared_ptr<LocalAssembler> &master_to_slave,
			const std::shared_ptr<LocalAssembler> &slave_to_master)
	: assemblers_{{master_to_slave, slave_to_master}}
	{}

	MixedBidirectionalLocalAssembler::~MixedBidirectionalLocalAssembler() {}

	bool MixedBidirectionalLocalAssembler::assemble(
		const Elem &master,
		FEType master_type,
		const Elem &slave,
		FEType slave_type,
		Matrix &mat
		) 
	{
		assert(false);
		return false;
	}

	bool MixedBidirectionalLocalAssembler::assemble(
		const Elem &master,
		FEType master_type,
		const Elem &slave,
		FEType slave_type,
		std::vector<Matrix> &mat
		) 
	{
		const std::size_t n = n_forms();
		mat.reserve(n);

		if(!assemblers_[0]->assemble(master, master_type, slave, slave_type, mat)) {
			return false;
		}

		std::vector<Matrix> temp;
		if(!assemblers_[1]->assemble(slave, slave_type, master, master_type, temp)) {
			return false;
		}

		mat.insert(mat.end(), temp.begin(), temp.end());
		return true;
	}

	int MixedBidirectionalLocalAssembler::n_forms() const
	{
		int ret = 0;
		for(const auto &a_ptr : assemblers_) {
			ret += a_ptr->n_forms();
		}

		return ret;
	}

	MixedBidirectionalLocalAssembler::Type MixedBidirectionalLocalAssembler::type(const int index) const
	{
		int range_end = assemblers_[0]->n_forms();

		if(index < range_end) {
			return assemblers_[0]->type(index);
		}

		assert(index < range_end + assemblers_[1]->n_forms());

		if(index == range_end) {
			return SLAVE_X_MASTER;
		}

		if(index == range_end + 1) {
			return MASTER_X_MASTER;
		}

		assert(false);
		m_utopia_error("should never happen");
		return MASTER_X_SLAVE;
	}

	void MixedBidirectionalLocalAssembler::print_stats(std::ostream &os) const
	{
		for(const auto &a_ptr : assemblers_) {
			a_ptr->print_stats(os);
		}
	}
}
