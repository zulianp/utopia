#ifndef UTOPIA_LINEAR_MULTI_LEVELHPP
#define UTOPIA_LINEAR_MULTI_LEVELHPP

#include "utopia_MultiLevelBase.hpp"
#include "utopia_Recorder.hpp"
#include "utopia_MatrixTransfer.hpp"
#include "utopia_MultiLevelMask.hpp"

#include <iostream>



namespace utopia
{

	/**
	 * @brief      Base class for all linear multilevel solvers. \n
	 *             Takes care of inializing multilevel hierarchy. \n
	 *             Different levels are created by interpolation and restriction operators.\n
	 *             Additionally, it provides stifness matrices on each level, created by using Galerkin assembly. \n
	 *
	 * @tparam     Matrix
	 * @tparam     Vector
	 */
	template<class Matrix, class Vector>
	class LinearMultiLevel : public MultiLevelBase<Matrix, Vector> {
		typedef UTOPIA_SCALAR(Vector)    Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
		typedef utopia::Level<Matrix, Vector> Level;
		typedef utopia::Transfer<Matrix, Vector> Transfer;
		typedef utopia::MatrixTransfer<Matrix, Vector> MatrixTransfer;
	public:

		using MultiLevelBase<Matrix, Vector>::set_transfer_operators;

		LinearMultiLevel(const Parameters params = Parameters())
		: MultiLevelBase<Matrix, Vector>(params),
		fix_semidefinite_operators_(false)
		{ }

		virtual ~LinearMultiLevel(){}

		/**
		 * @brief
		 * Function initializes restriction transfer operators.
		 * Operators need to be ordered FROM COARSE TO FINE.
		 *
		 * @param[in]  operators                The inteprolation operators.
		 *
		 */
		virtual bool set_transfer_operators(const std::vector<std::shared_ptr<Matrix>> &interpolation_operators)
		{
			this->transfers_.clear();
			for(auto I = interpolation_operators.begin(); I != interpolation_operators.end() ; ++I )
				this->transfers_.push_back(std::make_shared<MatrixTransfer>(*I));

			return true;
		}


		bool must_generate_masks()
		{
			return mask_.must_generate_masks();
		}

		void must_generate_masks(const bool must_generate_masks)
		{
			mask_.must_generate_masks(must_generate_masks);
		}

		void add_level(Level &&level)
		{
			levels_.push_back(std::move(level));
		}

		void add_level(const Level &level)
		{
			levels_.push_back(level);
		}

		static void fix_semidefinite_operator(Matrix &A)
		{

			Vector d;

			Size s = local_size(A);
			d = local_values(s.get(0), 1.);

			{
				Write<Vector> w_d(d);

				each_read(A,[&d](const SizeType i, const SizeType, const double) {
					d.set(i, 0.);
				});
			}

			A += Matrix(diag(d));
		}

		/**
		 * @brief
		 *        The function creates corser level operators provided by assembling on differnet levels of MG hierarchy
		 *		The first element of the vector is the coarset matrix and the last is the fienst
		 * @param[in]  stifness matrix for finest level
		 *
		 */
		inline bool set_linear_operators(const std::vector<std::shared_ptr<const Matrix>> &A)
		{
			levels_.clear();
			levels_.insert(levels_.begin(), A.begin(), A.end());
			return true;
		}

		inline void set_fix_semidefinite_operators(const bool val)
		{
			fix_semidefinite_operators_ = val;
		}


		inline const Level &level(const SizeType l) const
		{
			return levels_[l];
		}

		virtual void describe(std::ostream &os = std::cout) const override
		{
			SizeType i = 0;
			for(const auto &l : levels_) {
				const auto &A = l.A();
				auto s = size(A);

				os << "level: " << ++i << ", n_dofs: " << s.get(0) << std::endl;
			}
		}

		virtual void update_transfer(const SizeType level, const std::shared_ptr<Transfer> &t)
		{
			this->transfers_[level] = t;
		}


	protected:
		std::vector<Level>                  levels_;      /*!< vector of level operators     */
		MultiLevelMask<Matrix, Vector> mask_;
		bool fix_semidefinite_operators_;
		bool must_generate_masks_;

		/**
		 * @brief
		 *        The function creates corser level operators by using Galerkin assembly.
		 *
		 *        $\f J_{i-1} = R * J_{i} * I  $\f
		 *
		 *        Resulting operators are assumed to go from fines = 0 to coarse = numlevels_
		 *
		 * @param[in]  stifness matrix for finest level
		 *
		 */
		virtual bool galerkin_assembly(const std::shared_ptr<const Matrix> &A)
		{
			levels_.clear();
			SizeType t_s = this->transfers_.size();
			if(t_s <= 0)
				std::cerr<<"Provide interpolation operators first!  \n";

			levels_.push_back(Level(A));

			auto L = this->n_levels();

			mask_.generate_masks(*A, this->transfers_);


			for(SizeType i = 1; i < L; i++)
			{
				// J_{i-1} = R * J_{i} * I
				std::shared_ptr<Matrix> J_h = std::make_shared<Matrix>();
				this->transfer(t_s - i).restrict(levels_[i - 1].A(), *J_h);

				if(fix_semidefinite_operators_) {
					fix_semidefinite_operator(*J_h);
				}

				levels_.push_back(Level(J_h));
			}

			std::reverse(std::begin(levels_), std::end(levels_));
			return true;
		}

		virtual void apply_mask(const SizeType l, Vector &v) const
		{
			mask_.apply(l, v);
		}

	};

}



#endif //UTOPIA_LINEAR_MULTI_LEVELHPP
