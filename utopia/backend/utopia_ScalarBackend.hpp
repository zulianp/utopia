#ifndef UTOPIA_SCALAR_BACKEND_HPP
#define UTOPIA_SCALAR_BACKEND_HPP 

namespace utopia {

	template<typename Scalar>
	class ScalarBackend {
	public:
		static inline bool apply(const Scalar left, const Scalar right, const Plus &op, Scalar &result)
		{
			aux_apply(left, right, op, result);
			return true;
		}

		static inline bool apply(const Scalar left, const Scalar right, const Minus &op, Scalar &result)
		{
			aux_apply(left, right, op, result);
			return true;
		}


		static inline bool apply(const Scalar left, const Scalar right, const Multiplies &op, Scalar &result)
		{
			aux_apply(left, right, op, result);
			return true;
		}

		static inline bool apply(const Scalar left, const Scalar right, const EMultiplies &op, Scalar &result)
		{
			aux_apply(left, right, op, result);
			return true;
		}

		static inline bool apply(const Scalar left, const Scalar right, const Divides &op, Scalar &result)
		{
			aux_apply(left, right, op, result);
			return true;
		}

		static inline bool apply(const Scalar left, const Scalar right, const AbsPlus &op, Scalar &result)
		{
			aux_apply(left, right, op, result);
			return true;
		}

		static inline bool apply(const Scalar left, const Scalar right, const ApproxEqual &op, bool &result)
		{
			aux_apply(left, right, op, result);
			return true;
		}

		static inline bool zaxpy(const Scalar alpha, const Scalar x, const Scalar y, Scalar &result)
		{
			result = alpha * x + y;
			return true;
		}

	protected:
		ScalarBackend() {}


		template<typename Left, typename Right, typename Operation, typename Result>
		static inline void aux_apply(const Left &left, const Right &right, const Operation &op, Result &result)
		{
			result = op.template apply<Right>(left, right);
		}
	};

}

#endif //UTOPIA_SCALAR_BACKEND_HPP
