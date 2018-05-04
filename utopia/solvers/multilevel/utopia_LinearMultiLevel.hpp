#ifndef UTOPIA_LINEAR_MULTI_LEVELHPP
#define UTOPIA_LINEAR_MULTI_LEVELHPP

#include "utopia_MultiLevelBase.hpp"
#include "utopia_Recorder.hpp"




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
	public:
	
		LinearMultiLevel(const Parameters params = Parameters())
		: MultiLevelBase<Matrix, Vector>(params),
		  fix_semidefinite_operators_(false),
		  must_generate_masks_(false)
		{ }
		
		virtual ~LinearMultiLevel(){}

		/**
		 * @brief
		 * Function initializes restriction transfer operators.
		 * Operators need to be ordered FROM COARSE TO FINE.
		 *
		 * @param[in]  operators                The restriction operators.
		 *
		 */
		virtual bool set_transfer_operators(const std::vector<std::shared_ptr<Matrix>> &interpolation_operators)
		{
			this->transfers_.clear();	
			for(auto I = interpolation_operators.begin(); I != interpolation_operators.end() ; ++I )
				this->transfers_.push_back(Transfer(*I));
			
			return true;
		}
		

		bool must_generate_masks()
		{
			return must_generate_masks_; 
		}

		void must_generate_masks(const bool must_generate_masks)
		{
			must_generate_masks_ = must_generate_masks; 
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

		inline Transfer &transfer(const SizeType level)
		{
			return this->transfers_[level];
		}

		inline const Transfer &transfer(const SizeType level) const
		{
			return this->transfers_[level];
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

		virtual void update_transfer(const SizeType level, Transfer &&t)
		{
			this->transfers_[level] = std::move(t);
		}

		virtual void update_transfer(const SizeType level, const Transfer &t)
		{
			this->transfers_[level] = t;
		}
		
	protected:
		std::vector<Level>                  levels_;      /*!< vector of level operators     */		
		std::vector<Vector>					masks_;
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

			if(must_generate_masks_) {
				this->generate_masks(*A);
			}
			
			for(SizeType i = 1; i < L; i++)
			{
				// J_{i-1} = R * J_{i} * I
				std::shared_ptr<Matrix> J_h = std::make_shared<Matrix>();
				this->transfers_[t_s - i].restrict(levels_[i - 1].A(), *J_h);
				
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
			if(masks_.empty()) return;

			v = e_mul(masks_[l], v);
		}

	private:
		static void generate_mask_from_matrix(const Matrix &A, Vector &mask, const Scalar on_value, const Scalar off_value)
		{
			const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon() * 1e6;

			auto ls = local_size(A);
			mask = local_values(ls.get(0), off_value);

			{
				Write<Vector> w_(mask);

				each_read(A, [&](const SizeType i, const SizeType j, const Scalar value) {
					if(i == j) return;

					if(std::abs(value) > off_diag_tol) {
						mask.set(i, on_value);
					}
				});
			}
		}

		virtual void generate_masks(const Matrix &A)
		{
			const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon() * 1e6;

			const auto L = this->n_levels();
			masks_.resize(L);
			auto &mask = masks_[L - 1];

			generate_mask_from_matrix(A, mask, 0., 1.);
			// UTOPIA_RECORD_VALUE("r_mask", mask);

			for(SizeType l = L - 1; l > 0; --l) {
				auto &mask_l = masks_[l - 1];
				this->transfer(l-1).boolean_restrict_or(masks_[l], mask_l);

				// UTOPIA_RECORD_VALUE("r_mask", mask_l);
			}	

			for(auto &m : masks_) {
				m = local_values(local_size(m).get(0), 1.) - m;
				// UTOPIA_RECORD_VALUE("mask", m);
			}
		}
	};
	
}



#endif //UTOPIA_LINEAR_MULTI_LEVELHPP
