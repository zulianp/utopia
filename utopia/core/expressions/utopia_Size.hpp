//
// Created by Patrick Zulian on 26/05/15.
//

#ifndef UTOPIA_UTOPIA_SIZE_HPP
#define UTOPIA_UTOPIA_SIZE_HPP


#include "utopia_Base.hpp"
#include <vector>
#include <ostream>
#include <iostream>

namespace utopia {
    class Size {
    public:
        // typedef std::vector<SizeType>::size_type SizeType;
       typedef utopia::SizeType SizeType;


        /**    @defgroup   size Size
        *      @brief      Manipulation with sizes of objects.
        *      @ingroup    base_functions
        */



        /**
         * @ingroup    size
         * @brief       Returns size of the object. \n
         *              To obtain size of the object, use following syntax: M.size().get(d), where M represents tensor and d dimension of interest. 
         *
         */
        inline SizeType get(SizeType index) const 
        {
            assert(index < _data.size());
            return _data[index];
        }

        Size()
        {}

        explicit Size(const int n) : _data(n, 0)
        {
            assert(n <= UTOPIA_MAX_TENSOR_ORDER);
        }

        Size(std::initializer_list<SizeType> args)
        {
            using std::copy;
            _data.resize(args.size());
            copy(args.begin(), args.end(), _data.begin());
        }

        template<typename IntegerT>
        Size(std::initializer_list<IntegerT> args)
        {
            using std::copy;
            _data.resize(args.size());
            copy(args.begin(), args.end(), _data.begin());
        }

        void setDims(const SizeType n)
        {
            assert(n <= UTOPIA_MAX_TENSOR_ORDER);
            _data.resize(n);
        }

        void set_dims(const SizeType n)
        {
            assert(n <= UTOPIA_MAX_TENSOR_ORDER);
            _data.resize(n);
        }

        
        SizeType nDims() const
        {
            return n_dims();
        }


        /**
         * @ingroup     size
         * @brief       Returns dimension of the object. \n
         *              Use following syntax: long d = M.size().n_dims(), where M represents tensor and d is its dimension. 
         *
         */
        inline SizeType n_dims() const
        {
            return _data.size();
        }

        bool empty() const
        {
            return _data.empty();
        }

        void set(const SizeType index, const SizeType value)
        {
            assert(index < _data.size());
            _data[index] = value;
        }

        const std::vector<SizeType> &data() const {
            return _data;
        }

        friend bool operator==(const Size &l, const Size &r)
        {
            if(l.nDims() != r.nDims()) return false;

            for(SizeType i = 0; i < l.nDims(); ++i) {
                if(l.get(i) != r.get(i)) return false;
            }

            return true;
        }

        friend bool operator!=(const Size &l, const Size &r)
        {
            return !(l == r);
        }

    private:

        std::vector<SizeType> _data;
    };

    inline void disp(const Size &size, std::ostream &os)
    {
        disp(size.data().begin(), size.data().end(), os);
    }

    inline void disp(const Size &size)
    {
        return disp(size, std::cout);
    }

    inline Size kron_prod(const Size &left, const Size &right)
    {
        Size ret(left.nDims() + right.nDims());

        for(Size::SizeType i = 0; i < left.nDims(); ++i) {
            ret.set(i, left.get(i));
        }

        for(Size::SizeType i = 0; i < right.nDims(); ++i) {
            ret.set(left.nDims() + i, right.get(i));
        }

        return ret;
    }
}

#endif //UTOPIA_UTOPIA_SIZE_HPP
