#ifndef UTOPIA_CONSTRUCTIBLE_HPP
#define UTOPIA_CONSTRUCTIBLE_HPP

#include "utopia_Size.hpp"

namespace utopia {
	template<typename Scalar_, typename SizeType_, int Order_>
	class Constructible {};

	template<typename Scalar_, typename SizeType_>
	class Constructible<Scalar_, SizeType_, 2> {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		virtual ~Constructible() {}
		virtual void identity(const Size &s, const Scalar diag = 1.0) = 0;
		virtual void zeros(const Size &s) { values(s, 0.0); }
		virtual void values(const Size &s, const Scalar val) = 0;
			

		///Specialize for sparse matrices
		virtual void sparse(const Size &s, const SizeType /*nnz*/) 
		{
			zeros(s);
		}

		virtual void local_zeros(const Size &s) { values(s, 0.0); }
		virtual void local_values(const Size &s, const Scalar val) { values(s, val); }

		///Specialize for sparse matrices
		virtual void local_sparse(const Size &s, const SizeType /*nnz*/) 
		{
			local_zeros(s);
		}

	};

	template<typename Scalar_, typename SizeType_>
	class Constructible<Scalar_, SizeType_, 1> {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		virtual ~Constructible() {}
		virtual void zeros(const Size &s) { values(s, 0.0); }
		virtual void values(const Size &s, const Scalar val) { values(s.get(0), val); }


		virtual void zeros(const SizeType &s) { values(s, 0.0); }
		virtual void values(const SizeType &s, const Scalar val) = 0;

		virtual void local_zeros(const Size &s) { values(s, 0.0); }
		virtual void local_values(const Size &s, const Scalar val) { values(s, val); }
	};
}

#endif //UTOPIA_CONSTRUCTIBLE_HPP
