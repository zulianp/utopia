#ifndef UTOPIA_BLAS_VECTOR_HPP
#define UTOPIA_BLAS_VECTOR_HPP

#include "utopia_Interfaces.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_blas_Traits.hpp"
#include "utopia_BLAS_Operands.hpp"
#include "utopia_blas_Algorithms.hpp"

#include <vector>
#include <memory>
#include <iostream>

namespace utopia {
    template<typename T>
    class BlasVector : 
        public Vector<T, std::size_t>,
        public Tensor<BlasVector<T>, 2>,
        public BLAS1Tensor<BlasVector<T>>
    {
        typedef std::vector<T> Entries;
        typedef size_t SizeType;

	public:
        BlasVector(const BlasVector &other)
        : entries_(other.entries_)
        {}

        BlasVector(BlasVector &&other)
        : entries_(std::move(other.entries_))
        {}

        inline BlasVector &operator=(const BlasVector &other)
        {
        	if(this == &other) return *this;
        	copy(other);
        	return *this;
        }

        inline BlasVector &operator=(BlasVector &&other)
        {
        	if(this == &other) return *this;
        	entries_ = std::move(other.entries_);
        	return *this;
        }

        inline SizeType size()
        {
        	return entries_.size();
        }

        inline void resize(const SizeType n)
        {
        	entries_.resize(n);
        }

        inline T * ptr()
        {
        	return &entries_[0];
        }

        inline const T * ptr() const
        {
        	return &entries_[0];
        }


        auto begin() { return entries_.begin(); }
        auto end()   { return entries_.end(); }

        auto begin() const { return entries_.begin(); }
        auto end() const { return entries_.end(); }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR MatrixBase ///////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        //locks (No Op)
        inline void read_lock() override {}
        inline void write_lock(WriteMode) override {}
        inline void read_unlock() override {}
        inline void write_unlock(WriteMode) override {}

        //basic mutators
        inline void set(const SizeType &i, const T &value) override
        {
            assert(i < size());
            
            entries_[i] = value;
        }

        inline void add(const SizeType &i, const T &value) override
        {
            assert(i < size());
            
            entries_[i] += value;
        }

        //print function
        inline void describe() const override
        {
            auto &os = std::cout;
            describe(os);        }

        inline void describe(std::ostream &os) const
        {
        	for(auto e : entries_) {
        		os << e << " ";
        	}

        	os << std::endl;
        }

        //utility functions
        inline bool empty() const override
        {
            return entries_.empty();
        }

        inline void clear() override
        {
            entries_.clear();
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR DenseMatrix //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        inline T get(const SizeType &i) const override
        {
            assert(i < size());
            
            return entries_[i];
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR BLAS1Tensor //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////


        ///<T>SWAP - swap x and y
        inline void swap(BlasVector &x) override
        {
            std::swap(entries_, x.entries_);
        }

        ///<T>SCAL - x = a*x
        inline void scale(T &a) override
        {
            BLASAlgorithms<T>::scal(size(), a, ptr(), 1);
        }

        ///<T>COPY - copy x into y (this)
        inline void copy(const BlasVector &x) override
        {
            entries_.resize(x.size());
            BLASAlgorithms<T>::copy(x.size(), x.ptr(), 1, ptr(), 1);
        }

        ///<T>AXPY - y = a*x + y
        inline void axpy(const T &a, const BlasVector &x) override
        {
            assert(size() == x.size());
            BLASAlgorithms<T>::axpy(size(), a, x.ptr(), 1, ptr(), 1);
        }

        ///<T>DOT - dot product
        inline T dot(const BlasVector &other) const override
        {
            assert(size() == other.size());
            return BLASAlgorithms<T>::ddot(size(), ptr(), 1, other.ptr(), 1);
        }

        ///<T>NRM2 - Euclidean norm
        inline T norm2() const override
        {
            return BLASAlgorithms<T>::nrm2(size(), ptr(), 1);
        }

        ///<T>ASUM - sum of absolute values
        inline T asum() const override
        {
            return BLASAlgorithms<T>::asum(size(), ptr(), 1);
        }

        ///<T>AMAX - index of max abs value
        inline SizeType amax() const override
        {
            return BLASAlgorithms<T>::amax(size(), ptr(), 1);
        }

    private:
    	Entries entries_;

    };

}

#endif //UTOPIA_BLAS_VECTOR_HPP
