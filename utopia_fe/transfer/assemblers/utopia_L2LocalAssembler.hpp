#ifndef UTOPIA_LOCAL_L2_ASSEMBLER_HPP
#define UTOPIA_LOCAL_L2_ASSEMBLER_HPP

#include "utopia_LocalAssembler.hpp"
#include "MortarAssemble.hpp"

namespace utopia {
	class QMortarBuilder;

	class L2LocalAssembler final : public LocalAssembler {
	public:
		L2LocalAssembler(const int dim, const bool use_biorth);
		~L2LocalAssembler();

		/**
		* @brief if you are performing volume to surface transfer
		* the method does not provide reliable results if the volume 
		* element has facets aligned with the surface ones
		*/
		bool assemble(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			Matrix &mat
			) override;

		inline const QMortarBuilder &get_q_builder() const
		{
			assert(q_builder);
			return *q_builder;
		}

	private:
		int dim;
		bool use_biorth;
		bool must_compute_biorth;
		QMortar composite_ir;
		QMortar q_trial;
		QMortar q_test;

		Matrix biorth_weights;

		std::shared_ptr<QMortarBuilder> q_builder;
		std::unique_ptr<libMesh::FEBase> trial_fe, test_fe;

		void init_fe(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type);

		void init_biorth(const Elem &test, FEType test_type);

		static void assemble_biorth_weights(
			const libMesh::Elem &el,
			const int dim,
			const libMesh::FEType &var_type,
			const int el_order,
			libMesh::DenseMatrix<libMesh::Real> &weights);
	};
}


#endif //UTOPIA_LOCAL_L2_ASSEMBLER_HPP
