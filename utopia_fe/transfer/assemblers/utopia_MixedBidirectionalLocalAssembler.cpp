#include "utopia_MixedBidirectionalLocalAssembler.hpp"
#include "utopia.hpp"

namespace utopia {
	MixedBidirectionalLocalAssembler::	MixedBidirectionalLocalAssembler(
			const std::shared_ptr<LocalAssembler> &master_to_slave,
			const std::shared_ptr<LocalAssembler> &slave_to_master)
	: assemblers_{{master_to_slave, slave_to_master}}, return_false_only_if_both_fail_(false)
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
		mat.clear();
		mat.reserve(n);

		bool ok[2] = { true, true };
		if(!assemblers_[0]->assemble(master, master_type, slave, slave_type, mat)) {
			ok[0] = false;

			if(!return_false_only_if_both_fail_) {
				return false;
			}
		}

		std::vector<Matrix> temp;
		if(!assemblers_[1]->assemble(slave, slave_type, master, master_type, temp)) {
			ok[1] = false;

			if(!return_false_only_if_both_fail_ || !ok[0]) {
				return false;
			}

		}

		if(!ok[0]) {
			const int n = assemblers_[0]->n_forms();
				
			mat.resize(n);
			for(int i = 0; i < n; ++i) {
				mat[i].resize(0, 0);
				mat[i].get_values().clear();
			}

			mat.insert(mat.end(), temp.begin(), temp.end());
		} else {
			if(ok[1]) {
				mat.insert(mat.end(), temp.begin(), temp.end());
			} else {
				const int n = assemblers_[1]->n_forms();
				mat.resize(mat.size() + n, Matrix());
			}
		}

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
