#ifndef UTOPIA_BOX_CONSTRAINTS_HPP
#define UTOPIA_BOX_CONSTRAINTS_HPP 

#include <memory>

namespace utopia {

	template<class Vector>
	class BoxConstraints {
	public:
		BoxConstraints(const std::shared_ptr<Vector> &lower_bound,
					   const std::shared_ptr<Vector> &upper_bound)
		: lower_bound_(lower_bound), upper_bound_(upper_bound)
		{}

		BoxConstraints() {}

		inline std::shared_ptr<Vector> upper_bound()
		{
			return upper_bound_;
		}

		inline std::shared_ptr<const Vector> upper_bound() const
		{
			return upper_bound_;
		}

		inline std::shared_ptr<Vector> lower_bound()
		{
			return lower_bound_;
		}

		inline std::shared_ptr<const Vector> lower_bound() const
		{
			return lower_bound_;
		}

		inline bool has_lower_bound() const 
		{
			return static_cast<bool>(lower_bound_); 
		}

		inline bool has_upper_bound() const 
		{
			return static_cast<bool>(upper_bound_); 
		}

	private:
		std::shared_ptr<Vector> lower_bound_;
		std::shared_ptr<Vector> upper_bound_;		
	};

	template<class Vector>
	inline BoxConstraints<Vector> make_box_constaints(const std::shared_ptr<Vector> &lower_bound, 
											          const std::shared_ptr<Vector> &upper_bound)
	{
		return BoxConstraints<Vector>(lower_bound, upper_bound);
	}

	template<class Vector>
	inline BoxConstraints<Vector> make_lower_bound_constraints(const std::shared_ptr<Vector> &lower_bound)
	{
		return BoxConstraints<Vector>(lower_bound, nullptr);
	}
	
	template<class Vector>
	inline BoxConstraints<Vector> make_upper_bound_constraints(const std::shared_ptr<Vector> &upper_bound)
	{
		return BoxConstraints<Vector>(nullptr, upper_bound);
	}
}

#endif //UTOPIA_BOX_CONSTRAINTS_HPP
