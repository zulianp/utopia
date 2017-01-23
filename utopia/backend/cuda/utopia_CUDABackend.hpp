#ifndef UTOPIA_CUDA_BACKEND_HPP
#define UTOPIA_CUDA_BACKEND_HPP 

#include "utopia_Backend.hpp"
#include "utopia_CUDABackendTPL.hpp"

namespace utopia {

	template<>
	class Backend<double, utopia::CUDA> {
	public:
		typedef double Scalar;
		typedef utopia::CUDAVector<double> Vector;
		typedef utopia::CUDAMatrix<double> Matrix;

		static Backend &Instance() 
		{
			static Backend instance;
			return instance;
		}
		template<class Tensor>
		void describe(Tensor &&t)
		{
			cuda_double::describe(t);
		}

		template<class Left, class Right>
		void assign(Left &left, Right &&right)
		{
			left = std::forward<Right>(right);
		}

		void build(Matrix &m, const Size &size, const Identity &)
		{
			cuda_double::build_identity(size.get(0), size.get(1), m);
		}

		void build(Vector &v, const Size &size, const Values<double> &value)
		{
			cuda_double::build_values(size.get(0), value.value(), v);
		}

		template<typename VectorT>
		bool apply(const Matrix &m, VectorT &&v, const Multiplies &, Vector &result)
		{
			cuda_double::mat_vec_mul(m, v, result);
			return true;
		}

		template<class Tensor>
		double dot(const Tensor &left, const Tensor &right)
		{
			return cuda_double::dot(left, right);
		}

	private:
		Backend() {}
	};
}


#endif //UTOPIA_CUDA_BACKEND_HPP
