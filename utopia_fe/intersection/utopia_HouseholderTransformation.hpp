#ifndef UTOPIA_HOUSEHOLDER_TRANSFORMATION_HPP
#define UTOPIA_HOUSEHOLDER_TRANSFORMATION_HPP

#include "utopia_Intersect.hpp"

namespace utopia {

	class HouseholderTransformation {
	public:
		using Scalar = Intersector::Scalar;
		using Vector3 = Intersector::Vector3;
		using Vector2 = Intersector::Vector2;

		template<int Dim>
		using Vector = moonolith::Vector<Scalar, Dim>;

		HouseholderTransformation(const Vector3 &n)
		: trafo(3*3)
		{
			Scalar temp[3] = { n.x, n.y, n.z };
			Intersector::householder_reflection_3(&temp[0], &trafo[0]);
		}

		HouseholderTransformation(const Vector2 &n)
		: trafo(2*2)
		{
			Scalar temp[2] = { n.x, n.y };
			Intersector::householder_reflection_2(&temp[0], &trafo[0]);
		}

		template<int Dim>
		inline void apply(const Vector<Dim> &in, Vector<Dim> &out) const
		{
			for(int i = 0; i < Dim; ++i) {
				out[i] = 0.;

				for(int j = 0; j < Dim; ++j) {
					out[i] += trafo[i * Dim + j] * in[j];
				}
			}
		}

	private:
		std::vector<Scalar> trafo;
	};
}


#endif //UTOPIA_HOUSEHOLDER_TRANSFORMATION_HPP
