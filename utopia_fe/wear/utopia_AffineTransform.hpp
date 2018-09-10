#ifndef UTOPIA_AFFINE_TRANSFORM_HPP
#define UTOPIA_AFFINE_TRANSFORM_HPP

#include <array>
#include <algorithm>

namespace utopia {
	class AffineTransform {
	public:
		AffineTransform() {}

		template<std::size_t N>
		std::array<double, N> apply(const std::array<double, N> &p) const
		{
			std::array<double, N> ret;

			for(std::size_t i = 0; i < N; ++i) {
				ret[i] = translation[i];

				for(std::size_t k = 0; k < N; ++k) {
					ret[i] += get(i, k, linear_transformation) * p[k];
				}
			}

			return ret;
		}

		template<class Point>
		Point apply(const Point &p) const
		{
			Point ret;
			for(std::size_t i = 0; i < 3; ++i) {
				ret(i) = translation[i];

				for(std::size_t k = 0; k < 3; ++k) {
					ret(i) += get(i, k, linear_transformation) * p(k);
				}
			}

			return ret;
		}	

		void make_rotation(const int dim, const double angle, const char axis)
		{
			std::fill(begin(translation), end(translation), 0.0);
			std::fill(begin(linear_transformation), end(linear_transformation), 0.);

			switch(dim) {
				case 2:
				{
					assert(axis == 'y' || axis == 'Y');
					make_rotation_2(angle, linear_transformation);
					return;
				}
				case 3:
				{
					make_rotation_3(angle, axis, linear_transformation);
					return;
				}
				default:
				{
					assert(false);
					return;
				}
			}
		}

		void make_translation(const int dim, const double offset, const char axis)
		{
			std::fill(begin(translation), end(translation), 0.0);
			std::fill(begin(linear_transformation), end(linear_transformation), 0.);

			set(0, 0, 1., linear_transformation);
			set(1, 1, 1., linear_transformation);
			set(2, 2, 1., linear_transformation);

			switch(axis) {
				case 'X':
				case 'x':
				{
					translation[0] = offset;
					return;
				}

				case 'Y':
				case 'y':
				{
					translation[1] = offset;
					return;
				}

				case 'Z':
				case 'z':
				{
					translation[2] = offset;
					return;
				}

				default: {
					return;
				}
			}
		}

		AffineTransform operator*(const AffineTransform &other) const
		{
			AffineTransform ret;

			mat_mult(linear_transformation, other.linear_transformation, ret.linear_transformation);
			mat_mult(linear_transformation, other.translation, ret.translation);
			axpy(1., translation, ret.translation);
			return ret;
		}

		std::array<double, 3*3> linear_transformation;
		std::array<double, 3> translation;

	private:

		static void mat_mult(
			const std::array<double, 3*3> &A,
			const std::array<double, 3*3> &B,
			std::array<double, 3*3> &C)
		{

			std::fill(std::begin(C), std::end(C), 0.);

			for(int i = 0; i < 3; ++i) {
				auto i_offset = i * 3;

				for(int j = 0; j < 3; ++j) {
					for(int k = 0; k < 3; ++k) {
						auto k_offset = k * 3;
						C[i_offset + j] += A[i_offset + k] * B[k_offset + j];
					}
				}
			}
		}

		static void mat_mult(
			const std::array<double, 3*3> &A,
			const std::array<double, 3> &x,
			std::array<double, 3> &y)
		{

			std::fill(std::begin(y), std::end(y), 0.);

			for(int i = 0; i < 3; ++i) {
				auto i_offset = i * 3;

				for(int k = 0; k < 3; ++k) {
					y[i] += A[i_offset + k] * x[k];
				}
			}
		}

		static void axpy(
			const double alpha,
			const std::array<double, 3> &x,
			std::array<double, 3> &y)
		{
			std::fill(std::begin(y), std::end(y), 0.);

			for(int i = 0; i < 3; ++i) {
				y[i] += alpha * x[i];
			}
		}

		static void make_rotation_3(const double angle, const char axis, std::array<double, 9> &result)
		{
			set(0, 0, 1., result);
			set(1, 1, 1., result);
			set(2, 2, 1., result);

			if ((axis == 'x') || (axis == 'X')) {
				set(1, 1,  cos(angle), result);
				set(1, 2, -sin(angle), result);
				set(2, 1,  sin(angle), result);
				set(2, 2,  cos(angle), result);
			} else
				// set rotation around y axis:
			if ((axis == 'y') || (axis == 'Y')) {
				set(0, 0,  cos(angle), result);
				set(0, 2,  sin(angle), result);
				set(2, 0, -sin(angle), result);
				set(2, 2,  cos(angle), result);
			} else
					// set rotation around z axis:
			if ((axis == 'z') || (axis == 'Z')) {
				set(0, 0,  cos(angle), result);
				set(0, 1, -sin(angle), result);
				set(1, 0,  sin(angle), result);
				set(1, 1,  cos(angle), result);
			}
		}

		static void make_rotation_2(const double angle, std::array<double, 9> &result)
		{
			set(0, 0, cos(angle),  result);
			set(0, 1, -sin(angle), result);
			set(1, 0, sin(angle),  result);
			set(1, 1, cos(angle),  result);
		}

		inline static void set(
			const std::size_t i,
			const std::size_t j,
			const double value,
			std::array<double, 9> &mat)
		{
			mat[i * 3 + j] = value;
		}

		inline static double get(
			const std::size_t i,
			const std::size_t j,
			const std::array<double, 9> &mat)
		{
			return mat[i * 3 + j];
		}

	};
}


#endif //UTOPIA_AFFINE_TRANSFORM_HPP
