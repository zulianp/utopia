// #ifndef UTOPIA_PETSC_PROJECTED_GAUSS_SEIDEL_HPP
// #define UTOPIA_PETSC_PROJECTED_GAUSS_SEIDEL_HPP 

// #include "utopia_Base.hpp"

// #ifdef WITH_PASSO_EXTENSIONS

// #include "utopia_ForwardDeclarations.hpp"
// #include "utopia_BoxConstraints.hpp"
// #include "utopia_IterativeSolver.hpp"
// #include "utopia_Smoother.hpp"
// #include "utopia_ProjectedGaussSeidel.hpp"
// #include "utopia_petsc_Types.hpp"
// #include <memory>

// namespace utopia {

// 	template<>
// 	class ProjectedGaussSeidel<DSMatrixd, DVectord, PETSC> : public IterativeSolver<DSMatrixd, DVectord>, public Smoother<DSMatrixd, DVectord> {
// 	public:
// 		typedef utopia::DSMatrixd Matrix;
// 		typedef utopia::DVectord Vector;
// 		DEF_UTOPIA_SCALAR(Matrix)
// 		typedef Traits<DSMatrixd>::SizeType SizeType;

// 		typedef utopia::BoxConstraints<Vector>  BoxConstraints;

// 		virtual bool set_box_constraints(const BoxConstraints & box)
// 		{
// 			constraints_ = box;
// 			return true;
// 		}

// 		virtual void set_parameters(const Parameters params) override
// 		{
// 			IterativeSolver<Matrix, Vector>::set_parameters(params);
// 		}

// 		bool smooth(const Matrix &A, const Vector &b, Vector &x) override;
// 		bool apply(const Vector &b, Vector &x) override;

// 		void set_block_size(const PetscInt block_size)
// 		{
// 			block_size_ = block_size;
// 		}

// 		ProjectedGaussSeidel()
// 		: block_size_(1), use_line_search_(true)
// 		{}

// 		~ProjectedGaussSeidel();

// 		void init(const Matrix &A);

// 		void set_use_line_search(const bool val) 
// 		{
// 			use_line_search_ = val;
// 		}

// 		virtual void update(const std::shared_ptr<const Matrix> &op) override
// 		{
// 		    IterativeSolver<Matrix, Vector>::update(op);
// 		    init(*op);
// 		}
		
// 	private:
// 		PetscInt block_size_;
// 		bool use_line_search_;

// 		BoxConstraints constraints_;

// 		Vector active_set_;
// 		Vector prev_active_set_;
// 		Vector changed_;
// 		Vector old_x_;
		
// 		std::shared_ptr<void> ctx_;
// 	};
// }

// #endif //WITH_PASSO_EXTENSIONS
// #endif //UTOPIA_PETSC_PROJECTED_GAUSS_SEIDEL_HPP
