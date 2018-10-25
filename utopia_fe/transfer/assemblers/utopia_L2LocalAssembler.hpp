#ifndef UTOPIA_LOCAL_L2_ASSEMBLER_HPP
#define UTOPIA_LOCAL_L2_ASSEMBLER_HPP

#include "utopia_LocalAssembler.hpp"
#include "MortarAssemble.hpp"

#include <vector>

namespace utopia {
	class QMortarBuilder;

	class L2LocalAssembler final : public LocalAssembler {
	public:
		using Matrix = LocalAssembler::Matrix;

		L2LocalAssembler(
			const int dim,
			const bool use_biorth,
			const bool assemble_mass_mat = false,
			const bool is_shell = false);

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

		bool assemble(
			const Elem &trial,
			FEType trial_type,
			const Elem &test,
			FEType test_type,
			std::vector<Matrix> &mat
			) override;

		// bool volume_to_side_assemble(
		// 	const Elem &master,
		// 	FEType master_type,
		// 	const Elem &slave,
		// 	FEType slave_type,
		// 	const int slave_side_num,
		// 	std::vector<Matrix> &mat
		// ) override;

		inline const QMortarBuilder &get_q_builder() const
		{
			assert(q_builder);
			return *q_builder;
		}

		inline int n_forms() const override
		{
			return (assemble_mass_mat_)? 2 : 1;
		}

		inline Type type(const int index) const override
		{
			assert(index < n_forms());
			assert(index >= 0);
			
			return index == 0 ? MASTER_X_SLAVE : SLAVE_X_SLAVE;
		}

		void print_stats(std::ostream &os = std::cout) const override;

	private:
		int dim;
		bool use_biorth;
		bool must_compute_biorth;
		// QMortar composite_ir;
		QMortar q_trial;
		QMortar q_test;

		Matrix biorth_weights;

		std::shared_ptr<QMortarBuilder> q_builder;
		std::unique_ptr<libMesh::FEBase> trial_fe, test_fe;

		bool assemble_mass_mat_;
		int max_n_quad_points_;

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
