#ifndef UTOPIA_BACKEND_CONVERT_HPP
#define UTOPIA_BACKEND_CONVERT_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Range.hpp"
#include "utopia_Each.hpp"
#include "utopia_Factory.hpp"

#include <algorithm>
#include <numeric>

// namespace utopia {
// 	template<class TensorFrom, class TensorTo>
// 	void backend_convert_sparse(const Wrapper<TensorFrom, 2> &from, Wrapper<TensorTo, 2> &to)
// 	{
// 	    auto ls = local_size(from);
// 	    auto n_row_local = ls.get(0);
// 	    std::vector<int> nnzxrow(n_row_local, 0);
// 	    auto r = row_range(from);

// 	    each_read(from, [&nnzxrow,&r](const SizeType i, const SizeType j, const double val) {
// 	        ++nnzxrow[i - r.begin()];
// 	    });


// 	    auto nnz = *std::max_element(nnzxrow.begin(), nnzxrow.end());

// 	    //FIXME use nnzxrow instead
// 	    to = local_sparse(ls.get(0), ls.get(1), nnz);


// 	    Write<Wrapper<TensorTo, 2>> w_t(to);
// 	    each_read(from, [&to](const SizeType i, const SizeType j, const double val) {
// 	        to.set(i, j, val);
// 	    });
// 	}

// 	template<class TensorFrom, class TensorTo>
// 	void backend_convert(const Wrapper<TensorFrom, 1> &from, Wrapper<TensorTo, 1> &to)
// 	{
// 	    auto ls = local_size(from).get(0);
// 	    to = local_zeros(ls);

// 	    Write< Wrapper<TensorTo, 1> > w_t(to);
// 	    each_read(from, [&to](const SizeType i, const double val) {
// 	        to.set(i, val);
// 	    });
// 	}
// }


#endif //UTOPIA_BACKEND_CONVERT_HPP
