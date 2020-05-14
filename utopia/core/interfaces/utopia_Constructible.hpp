#ifndef UTOPIA_CONSTRUCTIBLE_HPP
#define UTOPIA_CONSTRUCTIBLE_HPP

#include "utopia_Size.hpp"

namespace utopia {
	template<typename Scalar_, typename SizeType_, int Order_>
	class Constructible {};

	template<typename Scalar_, typename SizeType_>
	class SparseConstructible {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		virtual ~SparseConstructible() {}

		virtual void sparse(const Size &s, const SizeType &/*nnz*/) = 0;
		virtual void local_sparse(const Size &s, const SizeType &/*nnz*/) = 0;

		virtual void identity(const Size &s, const Scalar &diag = 1.0) = 0;
		virtual void local_identity(const Size &s, const Scalar &diag = 1.0) { identity(s, diag); }
	};


	template<typename Scalar_, typename SizeType_>
	class DenseConstructible {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;
		
		virtual ~DenseConstructible() {}

		virtual void zeros(const Size &s) { values(s, 0.0); }
		virtual void values(const Size &s, const Scalar &val) = 0;
		virtual void dense_identity(const Size &s, const Scalar &diag = 1.0)  = 0;

		virtual void local_zeros(const Size &s) { local_values(s, 0.0); }
		virtual void local_values(const Size &s, const Scalar &val) { values(s, val); }
		virtual void local_dense_identity(const Size &s, const Scalar &diag = 1.0) { dense_identity(s, diag); }
	};


	template<typename Scalar_, typename SizeType_>
	class Constructible<Scalar_, SizeType_, 2> : 
		public SparseConstructible<Scalar_, SizeType_>,
		public DenseConstructible<Scalar_, SizeType_>  {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

                ~Constructible() override {}
                // virtual void identity(const Size &s, const Scalar &diag = 1.0) = 0;

		///Specialize for sparse matrices
                void sparse(const Size &s, const SizeType & /*nnz*/) override { zeros(s); }

                ///Specialize for sparse matrices
                void local_sparse(const Size &s, const SizeType & /*nnz*/) override { local_zeros(s); }

                void local_identity(const Size &s, const Scalar &diag = 1.0) override { this->identity(s, diag); }

                void zeros(const Size &s) override { this->values(s, 0.0); }
                // virtual void values(const Size &s, const Scalar &val) = 0;
                void dense_identity(const Size &s, const Scalar &diag = 1.0) override { this->identity(s, diag); }

                void local_zeros(const Size &s) override { local_values(s, 0.0); }
                void local_values(const Size &s, const Scalar &val) override { this->values(s, val); }
                void local_dense_identity(const Size &s, const Scalar &diag = 1.0) override { dense_identity(s, diag); }
        };

	template<typename Scalar_, typename SizeType_>
	class Constructible<Scalar_, SizeType_, 1> {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		virtual ~Constructible() {}
		
		virtual void zeros(const SizeType &s) { values(s, 0.0); }
		virtual void values(const SizeType &s, const Scalar &val) = 0;

		virtual void local_zeros(const SizeType &s) { local_values(s, 0.0); }
		virtual void local_values(const SizeType &s, const Scalar &val) { values(s, val); }

		//comodity
		virtual void zeros(const Size &s) { values(s, 0.0); }
		virtual void values(const Size &s, const Scalar &val) { values(s.get(0), val); }

		virtual void local_zeros(const Size &s) { local_values(s, 0.0); }
		virtual void local_values(const Size &s, const Scalar &val) { local_values(s.get(0), val); }
	};
}

#endif //UTOPIA_CONSTRUCTIBLE_HPP
