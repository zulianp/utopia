#ifndef UTOPIA_ML_BASE_HPP
#define UTOPIA_ML_BASE_HPP
#include "utopia_Transfer.hpp"
#include "utopia_Core.hpp"
#include "utopia_Input.hpp"

#include <algorithm>
#include <vector>

namespace utopia 
{
	// type of cycles used in Multilevel stuff 
	static const int MULTIPLICATIVE_CYCLE = 1;
	static const int ADDITIVE_CYCLE       = 2;
	static const int FULL_CYCLE           = 3;
	static const int NESTED_ITERATION     = 4;

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
	class MultiLevelBase : virtual public Configurable
	{
		typedef UTOPIA_SCALAR(Vector)    			Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) 			SizeType;
		typedef utopia::Level<Matrix, Vector> 		Level;
		typedef utopia::Transfer<Matrix, Vector> 	Transfer;
		typedef std::shared_ptr<Transfer> 			TransferPtr;

	public:

		MultiLevelBase(): 	pre_smoothing_steps_(3),
							post_smoothing_steps_(3),
							mg_type_(1), 
							cycle_type_(MULTIPLICATIVE_CYCLE), 
							v_cycle_repetition_(10),
							fix_semidefinite_operators_(false),
							num_levels_(-99)
		{

		}

		virtual ~MultiLevelBase(){}


		virtual void read(Input &in) override
        {
            in.get("pre_smoothing_steps", pre_smoothing_steps_);
            in.get("post_smoothing_steps", post_smoothing_steps_);
            in.get("mg_type", mg_type_);
            in.get("cycle_type", cycle_type_);
            in.get("v_cycle_repetition", v_cycle_repetition_);
            in.get("fix_semidefinite_operators", fix_semidefinite_operators_);
            in.get("num_levels", num_levels_);
        }


        virtual void print_usage(std::ostream &os) const override
        {
            this->print_param_usage(os, "pre_smoothing_steps", "int", "Number of pre-smoothing steps.", "3"); 
            this->print_param_usage(os, "post_smoothing_steps", "int", "Number of post-smoothing steps.", "3"); 
            this->print_param_usage(os, "mg_type", "int", "Multigrid type.", "1"); 
            this->print_param_usage(os, "cycle_type", "int", "Type of cycle used inside of multigrid method.", "MULTIPLICATIVE_CYCLE"); 
            this->print_param_usage(os, "v_cycle_repetition", "int", "Number of v-cycles used after one F-cycle.", "10"); 
            this->print_param_usage(os, "fix_semidefinite_operators", "bool", "Flag to fix semidefinite operators", "false"); 
            this->print_param_usage(os, "num_levels", "int", "Number of levels in ML hierarchy.", "-"); 
        }


		/**
		 * @brief Returns number of levels in hierarchy.
		 */
		inline SizeType n_levels() const 
		{
			return num_levels_; 
		}

		inline void n_levels(const SizeType & n_lev)
		{
			num_levels_ = n_lev; 
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

		inline void fix_semidefinite_operators(const bool val)
		{
			fix_semidefinite_operators_ = val;
		}

		inline bool fix_semidefinite_operators() const
		{
			return fix_semidefinite_operators_; 
		}

		bool set_transfer_operators(const std::vector<std::shared_ptr<utopia::Transfer<Matrix, Vector>>> &transfers)
		{
			if(num_levels_ <= 0) {
				num_levels_ = transfers.size() + 1; 
			}

			if(num_levels_ != transfers.size() + 1) {
				utopia_error("utopia::MultilevelBase:: number of levels and transfer operators do not match ... \n"); 
				std::cout << num_levels_ << " != " << (transfers.size() + 1);
			}

			transfers_ = transfers;
			return true;
		}


	protected:
		std::vector<TransferPtr> 	transfers_;   /*!< vector of transfer operators  */

		SizeType        			pre_smoothing_steps_;
		SizeType        			post_smoothing_steps_;
		SizeType        			mg_type_;

		SizeType         			cycle_type_;
		SizeType    				v_cycle_repetition_;

		bool 						fix_semidefinite_operators_;
		SizeType 					num_levels_; 
	};

}

#endif //UTOPIA_ML_BASE_HPP

