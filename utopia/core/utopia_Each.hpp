#ifndef UTOPIA_EACH_HPP
#define UTOPIA_EACH_HPP 

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Base.hpp"
#include "utopia_RowView.hpp"

namespace utopia {

	template<class Tensor, int Order = Tensor::Order, int FILL_TYPE = Tensor::FILL_TYPE>
	class Each {};

	template<class Tensor, int FILL_TYPE>
	class Each<Tensor, 1, FILL_TYPE> {
	public:
		template<class Fun>
		inline static void apply_read(const Tensor &v, Fun fun)
		{
			Range r = range(v);
			Read<Tensor> read_lock(v);

			for(auto i = r.begin(); i != r.end(); ++i) {
				fun(i, v.get(i));
			}
		}

		template<class Fun>
		inline static void apply_write(Tensor &v, Fun fun)
		{
			Range r = range(v);
			Write<Tensor> write_lock(v);

			for(auto i = r.begin(); i != r.end(); ++i) {
				v.set(i, fun(i));
			}
		}

		template<class Fun>
		inline static void apply_transform(const Tensor &in, Tensor &out, Fun fun)
		{
			assert(&in != &out && "in and out cannot be the same object");
			
			Range r = range(in);
			out = zeros(size(in));

			Read<Tensor>  read_lock(in);
			Write<Tensor> write_lock(out);

			for(auto i = r.begin(); i != r.end(); ++i) {
				out.set(i, fun(i, in.get(i)));
			}
		}
	};	

	template<class Tensor>
	class Each<Tensor, 2, FillType::DENSE> {
	public:
		template<class Fun>
		inline static void apply_read(const Tensor &m, Fun fun)
		{
			Range r = row_range(m);
			Range c = col_range(m);

			Read<Tensor> read_lock(m);

			for(auto i = r.begin(); i != r.end(); ++i) {
				for(auto j = c.begin(); j != c.end(); ++j) {
					fun(i, j, m.get(i, j));
				}
			}
		}

		template<class Fun>
		inline static void apply_write(Tensor &m, Fun fun)
		{
			Range r = row_range(m);
			Range c = col_range(m);
			Write<Tensor> write_lock(m);

			for(auto i = r.begin(); i != r.end(); ++i) {
				for(auto j = c.begin(); j != c.end(); ++j) {
					m.set(i, j, fun(i, j));
				}
			}
		}
	};

	template<class Tensor>
	class Each<Tensor, 2, FillType::SPARSE> {
	public:
		template<class Fun>
		inline static void apply_read(const Tensor &m, Fun fun)
		{
			Range r = row_range(m);
			for(auto i = r.begin(); i != r.end(); ++i) {
				RowView<const Tensor> row_view(m, i);
				for(auto index = 0; index < row_view.n_values(); ++index) {
					fun(i, row_view.col(index), row_view.get(index));
				}
			}
		}

		// template<class Fun>
		// inline static void apply_write(Tensor &m, Fun fun)
		// {
		// 	Range r = row_range(m);

		// 	//It will not be rewritten anyhow so this is safe
		// 	for(auto i = r.begin(); i != r.end(); ++i) {
		// 		RowView<Tensor> row_view(m, i);

		// 		for(auto index = 0; index < row_view.n_values(); ++index) {
		// 			row_view.set_value_at(index, fun(i, row_view.get_col_at(index)));
		// 		}
		// 	}
		// }

		// template<class Fun>
		// inline static void apply_transform(Tensor &m, Fun fun)
		// {
		// 	Range r = row_range(m);

		// 	//It will not be rewritten anyhow so this is safe
		// 	for(auto i = r.begin(); i != r.end(); ++i) {
		// 		RowView<Tensor> row_view(m, i);

		// 		for(auto index = 0; index < row_view.n_values(); ++index) {
		// 			auto val = fun(i, row_view.get_col_at(index), row_view.get_value_at(index)); 
		// 			row_view.set_value_at(index, val);
		// 		}
		// 	}
		// }
	};


     /** 	@defgroup element_acess Element Acess
     * 		@ingroup read_write
     *  	@brief  Actual acess to the elements of tensor
     */




	/**
	 * @ingroup element_acess
	 * @brief      Creates read lock on the tensor, iterates over all elements and applies provided function on them. \n
	 * 			   Example usage: Printing vector v. 
	 *
	 	\code{.cpp} 
	 	each_read(v, [](const SizeType i, const double entry) { std::cout << "v(" << i << ") = " << entry << std::endl;  });
	 	\endcode
	 * 
	 *
	 * @param[in]  v       The tensor. 
	 * @param[in]  fun     The  function with desirable action. 
	 */
	template<class Tensor, class Fun>
	inline void each_read(const Tensor &v, Fun fun) 
	{
		Each<Tensor>::apply_read(v, fun);
	}

	
	/**
	 * @ingroup element_acess
	 * @brief      Creates write lock on the tensor, iterates over all elements and applies provided function on them. \n
	 * 			   Example usage: Writing prescribed value to all elements of vector, but the first and the last. 
	 	\code{.cpp} 
	 	Vector v = zeros(10);
    		const double value = 6.0;

    		each_write(rhs, [value](const SizeType i) -> double {
				//The returned value will be written in the vector
				if(i == 0 || i == 10 - 1) {
					return 0;}
				return value;
			});   
		\endcode
	 * @warning    If tensor is a sparse matrix, it will iterate only following the sparsity pattern. 
	 * @param[in]  v       The tensor. 
	 * @param[in]  fun     The  function with desirable action. 
	 */
	template<class Tensor, class Fun>
	inline void each_write(Tensor &v, Fun fun) 
	{
		Each<Tensor>::apply_write(v, fun);
	}


	/**
	 * @ingroup element_acess
	 * @brief      Creates read lock on the tensor a and write lock on the tensor b, then applies provided function. 
	 * 			  
	 * 			  
	 * 			   Example usage: Applying a filter to the content of a and writing it into b. 
	 	
	 	\code{.cpp} 
		 	{ 
		 	const double n = 10; 
		        VectorT a = zeros(n);
		        VectorT b = zeros(n);
		        VectorT c = values(n, 0.5);  

		        //writing i/n in the vector
		        each_write(a, [](const SizeType i) -> double  { return i/double(n); }   );

		        {
		            //if another vector is needed, just provide a lock and pass it to the lambda functor
		            Read<VectorT> r(c);

		            //applying a filter to the content of a and writing it into b. a cannot be equal to b (for the moment)
		            each_transform(a, b, [&c](const SizeType i, const double entry) -> double  { return exp(-entry) * c.get(i); }    );
		        }
		    }
		\endcode
	 * @warning    Tensor a cannot be equal to the tensor b (for the moment). 
	 * @param[in]  a       The tensor to be red from/ input. 
	 * @param[in]  b       The tensor to be write into/ output. 
	 * @param[in]  fun     The  function with desirable action. 
	 */
	template<class Tensor, class Fun>
	inline void each_transform(const Tensor &a, Tensor &b, Fun fun) 
	{
		Each<Tensor>::apply_transform(a, b, fun);
	}
}

#endif //UTOPIA_EACH_HPP
