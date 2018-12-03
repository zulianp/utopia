#ifndef UTOPIA_L2_LOCAL_ASSEMBLER_HPP
#define UTOPIA_L2_LOCAL_ASSEMBLER_HPP

#include "utopia_LocalContactAssembler.hpp"
#include "utopia_ContactQMortarBuilder.hpp"

#include "libmesh/elem.h"
#include "libmesh/fe.h"

#include <cassert>
#include <ostream>
#include <memory>

namespace utopia {

	class L2LocalContactAssembler final : public LocalContactAssembler {
	public:
		using Elem   = libMesh::Elem;
		using FEType = libMesh::FEType;
		using Matrix = libMesh::DenseMatrix<libMesh::Real>;

		L2LocalContactAssembler(
			const int dim,
			const double search_radius,
			const bool use_biorth);
		~L2LocalContactAssembler() {}

		bool assemble(
			const Elem &master,
			const int master_side,
			FEType master_type,
			const Elem &slave,
			const int slave_side,
			FEType slave_type,
			Result &result
			) override;

		void print_stats(std::ostream &os = std::cout) const override;
	private:
		std::unique_ptr<ContactQMortarBuilder> q_builder_;
		std::unique_ptr<libMesh::FEBase> trial_fe_, test_fe_;
		QMortar q_trial_, q_test_;
		double search_radius_;
		bool use_biorth_;

		void init_fe(const Elem &trial,
				FEType trial_type,
				const Elem &test,
				FEType test_type);
		
	};
}

#endif //UTOPIA_L2_LOCAL_ASSEMBLER_HPP
