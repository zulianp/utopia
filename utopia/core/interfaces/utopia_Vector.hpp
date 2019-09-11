#ifndef UTOPIA_VECTOR_HPP
#define UTOPIA_VECTOR_HPP

#include "utopia_Enums.hpp"

namespace utopia {

	template<typename Scalar_, typename SizeType_>
	class VectorBase {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		virtual ~VectorBase() {}

		//locks
		virtual void read_lock() 			 	= 0;
		virtual void write_lock(WriteMode mode) = 0;

		virtual void read_unlock()  			  = 0;
		virtual void write_unlock(WriteMode mode) = 0;

		//basic mutators
		virtual void set(const SizeType &i, const Scalar &value) = 0;
		virtual void add(const SizeType &i, const Scalar &value) = 0;
		virtual Scalar get(const SizeType &i) = 0;

		//print function
		virtual void describe() const = 0;

		//utility functions
		virtual bool empty() const = 0;
		virtual bool clear() const = 0;
	};

	template<typename Scalar_, typename SizeType_>
	class Vector : public VectorBase<Scalar_, SizeType_> {
	public:
		virtual ~Vector() {}
	};

	//parallel types, collective operations
	template<typename Scalar_, typename SizeType_>
	class DistributedVector : public VectorBase<Scalar_, SizeType_> {
	public:
		using Scalar = Scalar_;
		using SizeType = SizeType_;

		//basic collective mutators allowing to write on other processes (e.g. for FE assembly)
		virtual void c_set(const SizeType &i, const Scalar &value) = 0;
		virtual void c_add(const SizeType &i, const Scalar &value) = 0;

		virtual ~DistributedVector() {}
	};
}

#endif //UTOPIA_VECTOR_HPP
