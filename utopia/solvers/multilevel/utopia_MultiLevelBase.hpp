#ifndef UTOPIA_ML_BASE_HPP
#define UTOPIA_ML_BASE_HPP
#include "utopia_Transfer.hpp"
#include "utopia_Core.hpp"

#include <algorithm>
#include <vector>

namespace utopia {
	/**
	 * @brief      Base class for all multilevel solvers. \n
	 *             Takes care of inializing multilevel hierarchy. \n
	 *             Different levels are created by interpolation and restriction operators.\n
	 *             Additionally, it provides stifness matrices on each level, created by using Galerkin assembly. \n
	 *
	 * @tparam     Matrix
	 * @tparam     Vector
	 */
	template<class Matrix, class Vector>
	class MultiLevelBase {
		typedef UTOPIA_SCALAR(Vector)    Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
		typedef utopia::Level<Matrix, Vector> Level;
		typedef utopia::Transfer<Matrix, Vector> Transfer;
		typedef std::shared_ptr<Transfer> TransferPtr;

	public:

		MultiLevelBase(const Parameters params = Parameters())
		{
			set_parameters(params);
		}

		virtual ~MultiLevelBase(){}

		virtual void set_parameters(const Parameters params)
		{
			parameters_ = params;
			pre_smoothing_steps_ = params.pre_smoothing_steps();
			post_smoothing_steps_ = params.post_smoothing_steps();
			mg_type_ = params.mg_type();
			cycle_type_ = params.cycle_type();
			v_cycle_repetition_ = 1;  // TODO:: create option in params for this
		}

		/**
		 * @brief Returns number of levels in hierarchy.
		 */
		inline SizeType n_levels()
		{
			return transfers_.size() + 1;
		}

		/**
		 * @brief      Function sets type of cycle
		 */
		inline bool cycle_type(const int &type_in)
		{
			cycle_type_ = type_in;
			return true;
		}

		/**
		 * @brief    Sets amount of V-cycles inside of F-cycle
		 */
		inline bool v_cycle_repetition(const SizeType & v_cycle_repetition_in)
		{
			v_cycle_repetition_ = v_cycle_repetition_in;
			return true;
		}

		/**
		 * @brief      Setting number pre-smoothing steps.
		 *
		 * @param[in]  pre_smoothing_steps_in  Number of pre-smoothing steps.
		 */
		inline void pre_smoothing_steps(const SizeType & pre_smoothing_steps_in )
		{
			pre_smoothing_steps_ = pre_smoothing_steps_in;
		}

		/**
		 * @brief      Setting number of post-smoothing steps.
		 *
		 * @param[in]  post_smoothing_steps_in  Number of post-smoothing steps.
		 */
		inline void post_smoothing_steps(const SizeType & post_smoothing_steps_in )
		{
			post_smoothing_steps_ = post_smoothing_steps_in;
		}

		/**
		 * @brief      Setting type of MG:  1 goes for V_CYCLE, 2 for W-cycle.
		 *
		 * @param[in]  mg_type_in  Choice of MG cycle.
		 */
		inline void mg_type(const bool mg_type_in)
		{
			mg_type_ = mg_type_in;
		}

		/**
		 * @return     Number of pre-smoothing steps.
		 */
		inline SizeType pre_smoothing_steps() const
		{
			return pre_smoothing_steps_;
		}

		/**
		 * @return     Number of post-smoothing steps.
		 */
		inline SizeType post_smoothing_steps() const
		{
			return post_smoothing_steps_;
		}

		/**
		 * @return     Type of MG cycle.
		 */
		inline bool mg_type() const { return mg_type_; }

		/**
		 * @return     Type of MG cycle.
		 */
		inline int cycle_type() const { return cycle_type_; }


		/**
		 * @brief      Amount of V-cycles on each level during full-cycle
		 */
		inline SizeType v_cycle_repetition() const { return v_cycle_repetition_; }

		inline Transfer &transfer(const SizeType level)
		{
			assert(level < transfers_.size());
			assert(transfers_[level]);

			return *transfers_[level];
		}

		inline const Transfer &transfer(const SizeType level) const
		{
			assert(level < transfers_.size());
			assert(transfers_[level]);
			
			return *transfers_[level];
		}

		virtual void describe(std::ostream &os = std::cout) const
		{
			(void) os;
		}

		bool set_transfer_operators(const std::vector<std::shared_ptr<utopia::Transfer<Matrix, Vector>>> &transfers)
		{
			transfers_ = transfers;
			return true;
		}


	protected:
		std::vector<TransferPtr>               transfers_;   /*!< vector of transfer operators  */

		Parameters                          parameters_;

		SizeType                            pre_smoothing_steps_;
		SizeType                            post_smoothing_steps_;
		SizeType                            mg_type_;

		int                                 cycle_type_;
		SizeType                            v_cycle_repetition_;

		bool fix_semidefinite_operators_;
	};

}

#endif //UTOPIA_ML_BASE_HPP

