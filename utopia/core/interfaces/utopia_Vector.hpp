#ifndef UTOPIA_VECTOR_HPP
#define UTOPIA_VECTOR_HPP

#include "utopia_Enums.hpp"
#include "utopia_Range.hpp"
#include "utopia_Layout.hpp"
#include "utopia_DistributedObject.hpp"
#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

	template<typename Scalar_, typename SizeType_>
	class Traits< DistributedVector<Scalar_, SizeType_> > {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;
		static const int Order = 1;
	};

	template<typename Scalar_, typename SizeType_>
	class VectorBase {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

                virtual ~VectorBase() = default;

                //locks
		virtual void read_lock() 			 	= 0;
		virtual void write_lock(WriteMode mode) = 0;

		virtual void read_unlock()  			  = 0;
		virtual void write_unlock(WriteMode mode) = 0;

		virtual void read_and_write_lock(WriteMode mode) = 0;
		virtual void read_and_write_unlock(WriteMode mode) = 0;

		//basic mutators
		virtual void set(const SizeType &i, const Scalar &value) = 0;
		virtual void add(const SizeType &i, const Scalar &value) = 0;
		virtual Scalar get(const SizeType &i) const = 0;

		//print function
		virtual void describe() const = 0;

		//utility functions
		virtual bool empty() const = 0;
		virtual void clear() = 0;
		virtual void set(const Scalar &val) = 0;

		virtual SizeType size() const = 0;
	};

	template<typename Scalar_, typename SizeType_>
	class Vector : public VectorBase<Scalar_, SizeType_> {
	public:
		using SizeType = SizeType_;

                ~Vector() override = default;

                // facade functions for treating a node vector and distributed vector in the same way
		inline Range range() const { return Range(0, this->size()); }

		// inline void layout(Layout<SizeType, 1> &l) const
		// {
		// 	l.local_size()  = local_size();
		// 	l.global_size() = local_size();
		// }

		inline SizeType local_size() const { return this->size(); }

	};

	//parallel types, collective operations
	template<typename Scalar_, typename SizeType_>
	class DistributedVector : public VectorBase<Scalar_, SizeType_>, public virtual DistributedObject {
	public:
		using Scalar   = Scalar_;
		using SizeType = SizeType_;

		//basic collective mutators allowing to write on other processes (e.g. for FE assembly)
		virtual void c_set(const SizeType &i, const Scalar &value) = 0;
		virtual void c_add(const SizeType &i, const Scalar &value) = 0;
		virtual Range range() const = 0;

		// virtual void layout(Layout<SizeType, 1> &l) const
		// {
		// 	l.local_size()  = local_size();
		// 	l.global_size() = this->size();
		// }

		virtual SizeType local_size() const = 0;

                ~DistributedVector() override = default;
        };

}

#endif //UTOPIA_VECTOR_HPP
