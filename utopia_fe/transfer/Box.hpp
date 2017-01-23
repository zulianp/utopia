#ifndef MFEM_L2P_BOX_HPP
#define MFEM_L2P_BOX_HPP 

#include "libmesh/libmesh.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include <libmesh/point.h>

namespace utopia {

	
	void get_row(const int row, const libMesh::DenseMatrix<libMesh::Real> &mat, libMesh::DenseVector<libMesh::Real> &result);
	void max_col(const libMesh::DenseMatrix<libMesh::Real> &mat, libMesh::DenseVector<libMesh::Real> &vec, bool include_vec_elements);
	void min_col(const libMesh::DenseMatrix<libMesh::Real> &mat, libMesh::DenseVector<libMesh::Real> &vec, bool include_vec_elements);

	class Box {
	public:

		Box(const int n);
		Box();
		
		virtual ~Box();

		void reset(const int n);
		void reset();

		Box & operator += (const libMesh::Point &point);
		Box & operator += (const Box &box);
		
		bool intersects(const Box &other) const;
		bool intersects(const Box &other, const double tol) const;
		
		void enlarge(const double value);

		void print(std::ostream &os = std::cout) const;

		inline double get_min(const int index) const
		{
			return min_(index);
		}

		inline double get_max(const int index) const
		{
			return max_(index);
		}

		inline const libMesh::DenseVector<libMesh::Real> &get_min() const
		{
			return min_;
		}

		inline const libMesh::DenseVector<libMesh::Real> &get_max() const
		{
			return max_;
		}
		
		inline libMesh::DenseVector<libMesh::Real> &get_min() 
		{
			return min_;
		}

		inline libMesh::DenseVector<libMesh::Real> &get_max() 
		{
			return max_;
		}
		
		inline int get_dims() const
		{
			return min_.size();
		}

		inline bool empty() const
		{
			if(min_.size() == 0) return true;
			return get_min(0) > get_max(0);
		}

	private:
		libMesh::DenseVector<libMesh::Real> min_;
		libMesh::DenseVector<libMesh::Real> max_;
	};

}

#endif //MFEM_L2P_BOX_HPP
