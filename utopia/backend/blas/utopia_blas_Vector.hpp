#ifndef UTOPIA_BLAS_VECTOR_HPP
#define UTOPIA_BLAS_VECTOR_HPP

#include "utopia_Interfaces.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_blas_Traits.hpp"
#include "utopia_BLAS_Operands.hpp"
#include "utopia_blas_Algorithms.hpp"
#include "utopia_Normed.hpp"
#include "utopia_Comparable.hpp"
#include "utopia_ElementWiseOperand.hpp"
#include "utopia_Transformable.hpp"
#include "utopia_Constructible.hpp"
#include "utopia_Reducible.hpp"
#include "utopia_blas_IndexSet.hpp"
#include "utopia_Allocations.hpp"
#include "utopia_Select.hpp"

#include <vector>
#include <memory>
#include <iostream>
#include <fstream>

namespace utopia {
    template<typename T>
    class BlasVector final :
        // Dynamic polymorphic types
        public Vector<T, std::size_t>,
        public Normed<T>,
        public Transformable<T>,
        public Reducible<T>,
        public Constructible<T, std::size_t, 1>,
        public ElementWiseOperand<T>,
        // Static polymorphic types
        public Tensor<BlasVector<T>, 1>,
        public Selectable<BlasVector<T>, 1>,
        public BLAS1Tensor<BlasVector<T>>,
        public Comparable<BlasVector<T>>,
        public ElementWiseOperand<BlasVector<T>>
    {
    public:
        typedef std::vector<T> Entries;

        using Scalar = T;
        using SizeType = std::size_t;

        using Constructible<T, std::size_t, 1>::values;
        using iterator = typename Entries::iterator;
        using const_iterator = typename Entries::iterator;

       ////////////////////////////////////////////////////////////////////
       ///////////////////////// BOILERPLATE CODE FOR EDSL ////////////////
       ////////////////////////////////////////////////////////////////////
        using Super = utopia::Tensor<BlasVector<T>, 1>;
        using Super::Super;

        template<class Expr>
        BlasVector(const Expression<Expr> &expr)
        {
            //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::construct_eval(expr.derived());
        }

        template<class Expr>
        inline BlasVector &operator=(const Expression<Expr> &expr)
        {
            Super::assign_eval(expr.derived());
            return *this;
        }

        void assign(const BlasVector &other) override
        {
            copy(other);
        }

        void assign(BlasVector &&other) override
        {
            entries_ = std::move(other.entries_);
        }

        ////////////////////////////////////////////////////////////////////



        BlasVector()
        {}

        BlasVector(std::initializer_list<T> args)
        : entries_(args)
        {
            UTOPIA_REPORT_ALLOC("BlasVector::BlasVector(std::initializer_list<T> args)");
        }

        BlasVector(const BlasVector &other)
        : entries_(other.entries_)
        {
            UTOPIA_REPORT_ALLOC("BlasVector::BlasVector(const BlasVector &other)");
        }

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

        inline SizeType size() const override
        {
        	return entries_.size();
        }

        inline void resize(const SizeType n)
        {
            UTOPIA_REPORT_ALLOC("BlasVector::resize(const SizeType n)");
        	entries_.resize(n);
        }

        inline void resize(const Size &s)
        {
            UTOPIA_REPORT_ALLOC("BlasVector::resize(const Size &s)");
            entries_.resize(s.get(0));
        }

        inline T * ptr()
        {
        	return &entries_[0];
        }

        inline const T * ptr() const
        {
        	return &entries_[0];
        }

        inline iterator begin() { return entries_.begin(); }
        inline iterator end()   { return entries_.end(); }

        inline const_iterator begin() const { return entries_.begin(); }
        inline const_iterator end()   const { return entries_.end(); }

        Entries &entries()
        {
        	return entries_;
        }

        const Entries &entries() const
        {
        	return entries_;
        }


        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR MatrixBase ///////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        //locks (No Op)
        inline void read_lock() override {}
        inline void write_lock(WriteMode) override {}
        inline void read_unlock() override {}
        inline void write_unlock(WriteMode) override {}

        inline void read_and_write_lock(WriteMode) override {}
        inline void read_and_write_unlock(WriteMode) override {}

        //basic mutators
        inline void set(const SizeType &i, const T &value) override
        {
            assert(i < size());

            entries_[i] = value;
        }

        virtual void set(const T &val) override
        {
            std::fill(std::begin(entries_), std::end(entries_), val);
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
        inline void scale(const T &a) override
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



        ///<T>AMAX - index of max abs value
        inline SizeType amax() const //override
        {
            return BLASAlgorithms<T>::amax(size(), ptr(), 1);
        }


        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Normed //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        inline T norm_infty() const override
        {
            using std::max_element;
            if(entries_.empty()) return 0.;
            return *max_element(entries_.begin(), entries_.end());
        }

        ///<T>NRM2 - Euclidean norm
        inline T norm2() const override
        {
            return BLASAlgorithms<T>::nrm2(size(), ptr(), 1);
        }

        ///<T>ASUM - sum of absolute values
        inline T norm1() const override
        {
            return BLASAlgorithms<T>::asum(size(), ptr(), 1);
        }
        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Comparable //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        bool equals(const BlasVector &other, const T &tol = 0.0) const override
        {
            const SizeType n = entries_.size();
            if(n != other.size()) return false;

            for(SizeType i = 0; i < n; ++i) {
                if(std::abs(entries_[i] - other.entries_[i]) > tol) {
                    return false;
                }
            }

            return true;
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR ElementWiseOperand //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        inline void e_mul(const BlasVector &other) override
        {
            const SizeType n = entries_.size();
            assert(n == other.size());

            for(SizeType i = 0; i < n; ++i) {
                entries_[i] *= other.entries_[i];
            }
        }

        inline void e_div(const BlasVector &other) override
        {
            const SizeType n = entries_.size();
            assert(n == other.size());

            for(SizeType i = 0; i < n; ++i) {
                entries_[i] /= other.entries_[i];
            }
        }

        inline void e_min(const BlasVector &other) override
        {
            const SizeType n = entries_.size();
            assert(n == other.size());

            for(SizeType i = 0; i < n; ++i) {
                entries_[i] = std::min(other.entries_[i], entries_[i]);
            }
        }

        inline void e_max(const BlasVector &other) override
        {
            const SizeType n = entries_.size();
            assert(n == other.size());

            for(SizeType i = 0; i < n; ++i) {
                entries_[i] = std::max(other.entries_[i], entries_[i]);
            }
        }

        /////////////////////////////////////////////

        inline void e_mul(const T &other) override
        {
            const SizeType n = entries_.size();

            for(SizeType i = 0; i < n; ++i) {
                entries_[i] *= other;
            }
        }

        inline void e_div(const T &other) override
        {
            const SizeType n = entries_.size();

            for(SizeType i = 0; i < n; ++i) {
                entries_[i] /= other;
            }
        }


        inline void e_min(const T &other) override
        {
            const SizeType n = entries_.size();

            for(SizeType i = 0; i < n; ++i) {
                entries_[i] = std::min(other, entries_[i]);
            }
        }

        inline void e_max(const T &other) override
        {
            const SizeType n = entries_.size();

            for(SizeType i = 0; i < n; ++i) {
                entries_[i] = std::max(other, entries_[i]);
            }
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Transformable //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////


        void transform(const Sqrt &op) override
        {
            aux_transform(op);
        }

        void transform(const Pow2 &op) override
        {
            aux_transform(op);
        }

        void transform(const Log &op) override
        {
            aux_transform(op);
        }

        void transform(const Exp &op) override
        {
            aux_transform(op);
        }

        void transform(const Cos &op) override
        {
            aux_transform(op);
        }

        void transform(const Sin &op) override
        {
            aux_transform(op);
        }

        void transform(const Abs &op) override
        {
            aux_transform(op);
        }


        void transform(const Pow &op) override
        {
            aux_transform(op);
        }

        void transform(const Reciprocal<T> &op) override
        {
            aux_transform(op);
        }

        void transform(const Minus &) override
        {
            for(auto &e : entries_) {
                e = -e;
            }
        }


        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Reducible //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        inline T reduce(const Min &op) const override
        {
            return aux_reduce(op);
        }

        inline T reduce(const Max &op) const override
        {
            return aux_reduce(op);
        }

        inline T reduce(const Plus &op) const override
        {
            return aux_reduce(op);
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Constructible //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        inline void values(const SizeType &s, const Scalar &val) override
        {
            resize(s);
            set(val);
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Selectable //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////


        inline void select(const BlasIndexSet &index, BlasVector &result) const override
        {
            const SizeType n = index.size();
            result.resize(n);

            for(SizeType i = 0; i < n; ++i) {
                result.set(i, get(index.get(i)));
            }
        }

        inline std::string get_class() const override {
            return "BlasVector";
        }


        inline bool has_nan_or_inf() const
        {
            for(const auto &e : entries_) {
                if(is_nan_or_inf(e)) return true;
            }

            return false;
        }

        inline bool is_alias(const BlasVector &other) const
        {
            return this == &other;
        }

        inline static SelfCommunicator &comm()
        {
            static SelfCommunicator instance_;
            return instance_;
        }

        inline bool write(const std::string &path) const
        {
            std::ofstream os(path.c_str());
            bool ok = os.good();

            if(ok) {
                const SizeType n = size();

                os << "x = [\n";

                for(SizeType i = 0; i < n; ++i) {
                    os << "\t" << get(i) << "\n";
                }

                os << "];\n";
            }

            os.close();
            return ok;
        }

    private:
    	Entries entries_;

        template<class Op>
        inline void aux_transform(const Op &op)
        {
            for(auto &e : entries_) {
                e = op.apply(e);
            }
        }

        template<class Op>
        inline T aux_reduce(const Op &op) const
        {
            if(entries_.empty()) return 0.0;

            const SizeType n = size();

            T ret = entries_[0];

            for(SizeType i = 1; i < n; ++i) {
                ret = op.apply(ret, entries_[i]);
            }

            return ret;
        }

    };

}

#endif //UTOPIA_BLAS_VECTOR_HPP
