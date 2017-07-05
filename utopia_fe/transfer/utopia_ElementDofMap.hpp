#ifndef UTOPIA_ELEMENT_DOF_MAP_HPP
#define UTOPIA_ELEMENT_DOF_MAP_HPP 

#include "cutk_Serializable.hpp"
#include "cutk_InputStream.hpp"
#include "cutk_OutputStream.hpp"

#include <vector>

namespace utopia {

	template<typename T>
	inline void write_vector(
		const std::vector<T> &v, 
		cutk::OutputStream &os)
	{
		int n = v.size();
		os << n;
		os.write(&v[0], n);
	}

	template<typename T>
	inline void read_vector(
		std::vector<T> &v, 
		cutk::InputStream &is)
	{
		int n;
		is >> n;
		v.resize(n);
		is.read(&v[0], n);
	}

	class ElementDofMap : public cutk::Serializable {
	public:
		inline void read(cutk::InputStream &is) override
		{
			read_vector(global, is);
		}

		inline void write(cutk::OutputStream &os) const override
		{
			write_vector(global, os);
		}

		inline bool empty() const
		{
			return global.empty();
		}

		std::vector<long> global;
	    //for(int side_number = 0; side_number < e.n_sides(); ++side_number)
	        //boundary_face_dof_map[side_number] gives the dofs associated to 
	        //if(boundary_face_dof_map[side_number].empty()) { means is not a boundary face }
	    //the boundary face
	    //std::vector< std::vector<long> > boundary_face_dof_map;

	    //get_side_dof_query(side, std::vector<bool> &is_boundary)
	};
}

#endif //UTOPIA_ELEMENT_DOF_MAP_HPP
