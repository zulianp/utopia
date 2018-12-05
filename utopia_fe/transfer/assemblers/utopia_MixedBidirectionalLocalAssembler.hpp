#ifndef UTOPIA_COMPOSITE_LOCAL_ASSEMBLER_HPP
#define UTOPIA_COMPOSITE_LOCAL_ASSEMBLER_HPP

#include "utopia_LocalAssembler.hpp"

#include "libmesh/dense_matrix.h"
#include <iostream>
#include <memory>

namespace utopia {

	class MixedBidirectionalLocalAssembler final : public LocalAssembler {
	public:
		using Elem   = LocalAssembler::Elem;
		using FEType = LocalAssembler::FEType;
		using Matrix = LocalAssembler::Matrix;

		MixedBidirectionalLocalAssembler(
			const std::shared_ptr<LocalAssembler> &master_to_slave,
			const std::shared_ptr<LocalAssembler> &slave_to_master
		);

		~MixedBidirectionalLocalAssembler();

		bool assemble(
			const Elem &master,
			FEType master_type,
			const Elem &slave,
			FEType slave_type,
			Matrix &mat
			) override;

		bool assemble(
			const Elem &master,
			FEType master_type,
			const Elem &slave,
			FEType slave_type,
			std::vector<Matrix> &mat
			) override;

		int n_forms() const override;

		Type type(const int index) const override;

		void print_stats(std::ostream &os = std::cout) const override;
	
	private:
		std::vector<std::shared_ptr<LocalAssembler>> assemblers_;
	};
}

#endif //UTOPIA_COMPOSITE_LOCAL_ASSEMBLER_HPP