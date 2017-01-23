#ifndef UTOPIA_NUMERICAL_TOLLERANGE_HPP
#define UTOPIA_NUMERICAL_TOLLERANGE_HPP   

namespace utopia {
	
	template<typename Scalar>
	class NumericalTollerance {
	public:
		NumericalTollerance(const Scalar absolute_tollerance,
							const Scalar relative_tollerance,
							const Scalar step_tollerance) 
		:
							absolute_tollerance_(absolute_tollerance),
							relative_tollerance_(relative_tollerance),
							step_tollerance_(step_tollerance)
		{}
		
		inline Scalar absolute_tollerance() const
		{
			return absolute_tollerance_;
		}

		inline Scalar relative_tollerance() const
		{
			return 	relative_tollerance_;
		}

		inline Scalar step_tollerance() const
		{
			return 	step_tollerance_;
		}

		inline void set_absolute_tollerance(const Scalar value) const
		{
			absolute_tollerance_ = value;
		}
		
		inline void set_relative_tollerance(const Scalar value) const
		{
			relative_tollerance_ = value;
		}

		inline void set_step_tollerance(const Scalar value) const
		{
			step_tollerance_ = value;
		}

		virtual void set_parameters(const Parameters params)
		{
			absolute_tollerance_ = params.atol();            
			relative_tollerance_ = params.rtol(); 
			step_tollerance_     = params.stol(); 
		}

	private:
		Scalar absolute_tollerance_;
		Scalar relative_tollerance_;
		Scalar step_tollerance_;
	};

}


#endif //UTOPIA_NUMERICAL_TOLLERANGE_HPP
