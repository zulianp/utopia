#ifndef MFEML2P_MORTAR_ASSEMBLER_HPP
#define MFEML2P_MORTAR_ASSEMBLER_HPP 


#include "utopia_fe_core.hpp"
#include "utopia_LibMeshBackend.hpp"
#include <libmesh/sparse_matrix.h>
#include "moonolith_predicate.hpp"

namespace utopia {

	class MortarAssembler {
	public:

		MortarAssembler(
			const std::shared_ptr<LibMeshFESpaceBase> &master, 
			const std::shared_ptr<LibMeshFESpaceBase> &slave);

		 bool assemble(DSMatrixd &B);

		 	DSMatrixd D;

		 	void set_use_biorthogonal_multipliers(const bool use_biorth)
		 	{
		 		use_biorth_ = use_biorth;
		 	}

		private:
			std::shared_ptr<LibMeshFESpaceBase> master_;
			std::shared_ptr<LibMeshFESpaceBase> slave_;	

			bool use_biorth_;
	};

	static const ushort MASTER     = 1;
	static const ushort SLAVE      = 2;
	static const ushort UNASSIGNED = 0;
	static const ushort REMOVED    = 3;

	class Contact {
	public:
		typedef libMesh::Real RealT;

		long parent_element_master;
		int side_number_master;
		long id_master;

		long parent_element_slave;
		int side_number_slave;
		long id_slave;

		RealT relative_area;
		RealT isect_area;
		
		RealT avg_gap;

		libMesh::DenseMatrix<RealT> coupling;
		libMesh::DenseVector<RealT> gap;
		libMesh::DenseMatrix<RealT> normals;

		bool is_valid; 

		friend bool operator<(const Contact &left, const Contact &right) {
			if(!left.is_valid) {
				return false;
			}

			if(!right.is_valid) {
				return true;
			}

			if(std::abs(left.avg_gap) < std::abs(right.avg_gap)) {
				return true;
			}

			if(std::abs(left.avg_gap) > std::abs(right.avg_gap)) {
				return false;
			}

			if(left.id_slave < right.id_slave) {
				return true;
			}

			if(left.id_slave > right.id_slave) {
				return false;
			}

			return left.id_master < right.id_master;
		}

		void describe(std::ostream &os = std::cout) const
		{
			os << "-------------------------------------\n";
			os << id_master << ", " << id_slave << "\n";
			os << "gap:\n";
			gap.print(os);
			os << "\n";
			os << "relative_area: " << relative_area << "\n";
			os << "normals:\n";
			normals.print(os);
			os << "\n";

			os << "coupling:\n";
			coupling.print(os);
			os << "\n";
			os << "-------------------------------------\n";
		}

		void finalize();
	};

	static void build_dag(std::vector< std::shared_ptr<Contact> > &contacts, std::vector< std::vector<long> > &dag, std::vector<long> &ordering);
	static void assign_master_and_slave_roles(const std::vector< std::vector<long> > &dag, const std::vector<long> &ordering, const std::vector< std::vector<long> > &adj_list, std::vector<ushort> &role);

	class MortarContactAssembler {
	public:

		//pass to true to accept pairs whose average distance is smaller than the gap
		void set_strict_gap_policy(const bool strict_gap_policy)
		{
			this->strict_gap_policy = strict_gap_policy;
		}
		

		MortarContactAssembler(const std::shared_ptr<LibMeshFESpaceBase> &space);
		bool assemble(DSMatrixd &coupling, DVectord &gap, DVectord &normals, DSMatrixd &orthogonal_trafos, std::vector<bool> &is_contact_node, const libMesh::Real search_radius, const std::shared_ptr<moonolith::Predicate> &predicate = std::shared_ptr<moonolith::Predicate>());

		private:
			std::shared_ptr<LibMeshFESpaceBase> space_;
			bool strict_gap_policy;

			
	};

}


#endif //MFEML2P_MORTAR_ASSEMBLER_HPP
