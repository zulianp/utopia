#ifndef UTOPIA_BOX_CONSTRAINTS_HPP
#define UTOPIA_BOX_CONSTRAINTS_HPP 

#include <memory>

namespace utopia {
	template<class Vector>
	class BoxConstraints {
	public:
		BoxConstraints(const std::shared_ptr<Vector> &lower_bound,
					   const std::shared_ptr<Vector> &upper_bound)
		: lower_bound_(lower_bound), upper_bound_(upper_bound), has_lower_(true), has_upper_(true)
		{}


		BoxConstraints(const std::shared_ptr<Vector> &bound, const std::string type = "upper")
		{
			if(type == "upper")
			{
				upper_bound_ = bound; 
				has_lower_   = false; 
				has_upper_	 = true; 
			}
			else if(type == "lower")
			{
				lower_bound_ = bound; 
				has_lower_   = true; 
				has_upper_	 = false; 
			}
			else
			{
				has_lower_   = false; 
				has_upper_	 = false; 
				std::cout<<"err: utopia_BoxConstraints:: check which kind of constraints are you setting .... \n"; 
			}
		}

		inline std::shared_ptr<Vector> &upper_bound()
		{
			return upper_bound_;
		}

		inline std::shared_ptr<const Vector> upper_bound() const
		{
			return upper_bound_;
		}

		inline std::shared_ptr<Vector> &lower_bound()
		{
			return lower_bound_;
		}

		inline std::shared_ptr<const Vector> lower_bound() const
		{
			return lower_bound_;
		}


		inline bool has_lower_bound() const 
		{
			return has_lower_; 
		}

		inline bool has_upper_bound() const 
		{
			return has_upper_; 
		}


	private:
		std::shared_ptr<Vector> lower_bound_;
		std::shared_ptr<Vector> upper_bound_;

		bool has_lower_; 
		bool has_upper_; 
		
	};

	template<class Vector>
	inline BoxConstraints<Vector> make_box_constaints(const std::shared_ptr<Vector> &lower_bound, 
											   const std::shared_ptr<Vector> &upper_bound)
	{
		return BoxConstraints<Vector>(lower_bound, upper_bound);
	}
	
}

#endif //UTOPIA_BOX_CONSTRAINTS_HPP
