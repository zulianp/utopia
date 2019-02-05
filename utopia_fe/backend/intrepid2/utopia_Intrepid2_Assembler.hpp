#ifndef UTOPIA_INTREPID_2_ASSEMBLER_HPP
#define UTOPIA_INTREPID_2_ASSEMBLER_HPP

#include "utopia_libmesh_AssemblyContext.hpp"
#include "utopia_make_unique.hpp"

// Intrepid2 includes
#include <Intrepid2_FunctionSpaceTools.hpp>
#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_ArrayTools.hpp>
#include <Intrepid2_HGRAD_HEX_C1_FEM.hpp>
#include <Intrepid2_HGRAD_TRI_C1_FEM.hpp>
#include <Intrepid2_RealSpaceTools.hpp>
#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_Utils.hpp>
#include <KokkosSparse_CrsMatrix.hpp>

// Teuchos includes
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_RCP.hpp>
// #include <Teuchos_BLAS.hpp>
// #include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>
// #include <Teuchos_XMLParameterListCoreHelpers.hpp>

// // Shards includes
// #include "Shards_CellTopology.hpp"

// Tpetra includes
#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>

namespace utopia {
	class Intrepid2Assembler {
	public:
		typedef utopia::USparseMatrix GlobalMatrix;
		typedef utopia::UVector GlobalVector;
		typedef UTOPIA_SCALAR(GlobalVector) Scalar;

		typedef Kokkos::DefaultExecutionSpace Host;
		typedef Kokkos::DefaultExecutionSpace Device;
		typedef Intrepid2::FunctionSpaceTools<Device> FST;
		typedef shards::CellTopology CellTopology;

		typedef Kokkos::DynRankView<Scalar, Device>    DeviceRankView;
		typedef Kokkos::DynRankView<Scalar, Host>    HostRankView;
		
		template<typename T>
		using TensorView2 = Kokkos::DualView<T**, Kokkos::LayoutRight>;

		template<typename T>
		using TensorView3 = Kokkos::DualView<T***, Kokkos::LayoutRight>;

		template<typename T, typename Dev>
		static void disp(const Kokkos::DynRankView<T, Dev> &v) {
			auto size = v.size();

			for(auto i = 0; i < size; ++i) {
				std::cout << i << " " << v(i) << std::endl;
			}
		}

		class IntrepidCubature {
		public:
			Intrepid2::DefaultCubatureFactory  cub_factory_;

			inline auto get_cubature(CellTopology &cell, const int order) -> decltype( this->cub_factory_.create<Device, Scalar, Scalar>(cell, order) )
			{
				return cub_factory_.create<Device, Scalar, Scalar>(cell, order);
			}

			inline void init(CellTopology &cell, const int order)
			{
				q = get_cubature(cell, order);
				q_points_  = utopia::make_unique<DeviceRankView>("q_points",  q->getNumPoints(), q->getDimension());
				q_weights_ = utopia::make_unique<DeviceRankView>("q_weights", q->getNumPoints());
				q->getCubature(*q_points_, *q_weights_);
			}

			Teuchos::RCP<Intrepid2::Cubature<Device, Scalar, Scalar> > q;
			std::unique_ptr<DeviceRankView> q_points_;
			std::unique_ptr<DeviceRankView> q_weights_;
		};

	
		

		class IntrepidFunctionSpace {
		public:

			void init(const LibMeshFunctionSpace &space) {
				auto &m           = space.mesh();
				auto dim 	      = m.mesh_dimension();
				auto num_nodes    = m.n_local_nodes(); //FIXME include ghost nodes
				num_elements      = m.n_active_local_elem();

				cell_ptr_     = std::make_shared<CellTopology>(shards::getCellTopologyData<shards::Triangle<3>>());
				points_ 	  = TensorView3<Scalar>("points", num_elements, cell_ptr_->getNodeCount(), dim);
				elem_to_node_ = TensorView2<int>("elem2node", num_elements, cell_ptr_->getNodeCount());

				std::vector<libMesh::dof_id_type> indices;

				auto e2n = elem_to_node_.template view<Kokkos::HostSpace>();
				auto pts = points_.template view<Kokkos::HostSpace>();

				int index = 0;
				for(auto e_it = m.active_local_elements_begin(); e_it != m.active_local_elements_end(); ++e_it, ++index) {
					const auto &e = **e_it;

					space.dof_map().dof_indices(&e, indices);
					assert(cell_ptr_->getNodeCount() == indices.size());
					
					for(auto i = 0; i < e.n_nodes(); ++i) {
						const auto &node = e.node_ref(i);
						e2n(index, i) = indices[i];

						for(auto d = 0; d < dim; ++d) {
							pts(index, i, d) = node(d);
						}
					}
				}

				points_.sync<Device>();
				elem_to_node_.sync<Device>();
			}

			IntrepidCubature &get_cubature(const int order)
			{
				cub_.init(*cell_ptr_, order);
				return cub_;
			}

			std::shared_ptr<CellTopology> cell_ptr_;
			
			TensorView3<Scalar> points_; 
			TensorView2<int> elem_to_node_; 

			IntrepidCubature cub_;
			int num_elements;
		};


		class IntrepidFE {
		public:
			Intrepid2::Basis_HGRAD_TRI_C1_FEM<Device, Scalar, Scalar> fe;

			void init(int order, IntrepidFunctionSpace &space)
			{
				try {
					auto &q = space.get_cubature(order);

					int num_elements 	= space.num_elements;
					int num_fields 		= fe.getCardinality();
					int num_quad_points = q.q->getNumPoints();
					int dim 			= q.q->getDimension();


					std::cout << "num_elements: "    << num_elements    << std::endl;
					std::cout << "num_fields: "      << num_fields      << std::endl;
					std::cout << "num_quad_points: " << num_quad_points << std::endl;
					std::cout << "dim: "             << dim             << std::endl;
					std::cout << std::flush;

					fun           = DeviceRankView("fun", num_fields, num_quad_points);
					ref_grad      = DeviceRankView("ref_grad", num_fields, num_quad_points, dim);
					
					grad     	  = DeviceRankView("grad", num_elements, num_fields, num_quad_points, dim);
					test_weighted = DeviceRankView("test_weighted", num_elements, num_fields, num_quad_points, dim);


					jacobian     = DeviceRankView("jacobian", num_elements, num_quad_points, dim, dim);
					jacobian_inv = DeviceRankView("jacobian_inv", num_elements, num_quad_points, dim, dim);
					jacobian_det = DeviceRankView("jacobian_det", num_elements, num_quad_points);
					dx 			 = DeviceRankView("dx", num_elements, num_quad_points);

					auto &q_points  = *q.q_points_;
					auto &q_weights = *q.q_weights_;

					fe.getValues(fun, 	   q_points, Intrepid2::OPERATOR_VALUE);
					fe.getValues(ref_grad, q_points, Intrepid2::OPERATOR_GRAD);	

					DeviceRankView nodes = space.points_.template view<Device>();

					Intrepid2::CellTools<Device>::setJacobian(jacobian, q_points, nodes, *space.cell_ptr_);
					Intrepid2::CellTools<Device>::setJacobianInv(jacobian_inv, jacobian);
					Intrepid2::CellTools<Device>::setJacobianDet(jacobian_det, jacobian);
					FST::HGRADtransformGRAD<double>(grad, jacobian_inv, ref_grad);
					FST::computeCellMeasure<double>(dx, jacobian_det, q_weights);

					disp(jacobian_det);
					// disp(dx);
				} catch(const std::exception &ex) {
					std::cerr << ex.what() << std::endl;
					abort();
				}
			}

			DeviceRankView fun;
			DeviceRankView ref_grad;
			DeviceRankView grad;

			DeviceRankView jacobian;
			DeviceRankView jacobian_inv;
			DeviceRankView jacobian_det;

			DeviceRankView dx;
			DeviceRankView test_weighted;
		};

		Intrepid2Assembler()
		: verbose_(Utopia::instance().verbose())
		{}

		class ElementMatrix {
		public:
			void init(IntrepidFunctionSpace &space, IntrepidFE &fe)
			{
				int num_elements 	= space.num_elements;
				int num_fields 		= fe.fe.getCardinality();
				mat = DeviceRankView("element_matrix", num_elements, num_fields, num_fields);
			}

			DeviceRankView mat;
		};

		class ElementVector {
		public:
			void init(IntrepidFunctionSpace &space, IntrepidFE &fe)
			{
				int num_elements 	= space.num_elements;
				int num_fields 		= fe.fe.getCardinality();
				vec = DeviceRankView("element_matrix", num_elements, num_fields);
			}

			DeviceRankView vec;
		};

		class QuadratureData {
		public:
			void init(IntrepidFunctionSpace &space)
			{
				int num_elements 	= space.num_elements;
				int num_quad_points = space.cub_.q->getNumPoints();
				data = DeviceRankView("element_matrix", num_elements, num_quad_points);
			}

			void set(IntrepidFunctionSpace &space, const Scalar value) 
			{
				int num_quad_points = space.cub_.q->getNumPoints();
				Kokkos::parallel_for(space.num_elements, KOKKOS_LAMBDA(int k) {
					for(int i = 0; i < num_quad_points; ++i) {
						data(k, i) = value;
					}
				});
			}

			DeviceRankView data;
		};

		//FIXME put in utopia
		template<class T>
		static bool is_ghosted(const Wrapper<T, 1> &vec)
		{
			return vec.implementation().has_ghosts();
		}

		template<class Expr>
		bool assemble(const Expr &expr, Scalar &val)
		{
			//perf
			Chrono c;
			c.start();

			typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
			
			typedef typename TraitsT::Vector ElementVector;

			static const int Backend = TraitsT::Backend;

			const auto &space = find_space<LibMeshFunctionSpace>(expr);
			const auto &dof_map = space.dof_map();
			auto &m = space.mesh();

			i_space.init(space);

			IntrepidFE fe;
			fe.init(0, i_space);


			val = 0.;

			//assemble energy

			//global redu e	
			m.comm().sum(val);

			//perf
			c.stop();

			if(verbose_) {
				std::cout << "assemble: value" << std::endl;
				std::cout << c << std::endl;
			}

			return false;
		}


		template<class Expr>
		bool assemble(const Expr &expr, GlobalMatrix &mat)
		{
			//perf
			Chrono c;
			c.start();

			typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
			
			static const int Backend = TraitsT::Backend;

			const auto &space = find_space<LibMeshFunctionSpace>(expr);
			const auto &dof_map = space.dof_map();
			auto &m = space.mesh();

			i_space.init(space);

			IntrepidFE fe;
			fe.init(0, i_space);


			ElementMatrix elem_mat;
			elem_mat.init(i_space, fe);

			auto s_m = size(mat);
			if(empty(mat) || s_m.get(0) != dof_map.n_dofs() || s_m.get(1) != dof_map.n_dofs()) {
				SizeType nnz_x_row = 0;
				if(!dof_map.get_n_nz().empty()) {
					nnz_x_row = 
					*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()) + 
					*std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end());
				}

				mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

				std::vector<libMesh::dof_id_type> indices;
				typename TraitsT::Matrix zero;

				Write<GlobalMatrix> w_(mat);
				for(auto e_it = m.active_local_elements_begin(); e_it != m.active_local_elements_end(); ++e_it) {
					dof_map.dof_indices(*e_it, indices);
					zero = zeros(indices.size(), indices.size());
					add_matrix(zero.implementation(), indices, indices, mat);
				}
				//BUILD the crs graph
			} else {
				mat *= 0.;
			}

			//assemble Hessian
			FST::multiplyMeasure<Scalar>(fe.test_weighted, fe.dx, fe.grad);
			FST::integrate(elem_mat.mat, fe.grad, fe.test_weighted);

			auto local_mat = raw_type(mat)->getLocalMatrix();

			int num_fields = fe.fe.getCardinality();

			auto elem_to_node = i_space.elem_to_node_.template view<Device>();

			// i_space.elem_to_node_.print();
			disp(elem_mat.mat);

			Kokkos::parallel_for(i_space.num_elements, KOKKOS_LAMBDA(const int k) {
			  for (int row = 0; row < num_fields; row++){
			      for (int col = 0; col < num_fields; col++){
			         
			          int row_index = elem_to_node(k, row);
			          int col_index = elem_to_node(k, col);

			          Scalar val = elem_mat.mat(k, row, col);

			          local_mat.sumIntoValues(
			            row_index,
			            &col_index,
			            1,
			            &val,
			            true,
			            true
			          );
			      }
			  	}
			});

			//perf
			c.stop();

			if(verbose_) {
				std::cout << "assemble: lhs" << std::endl;
				std::cout << c << std::endl;
			}

			return true;
		}


		template<class Expr>
		bool assemble(const Expr &expr, GlobalVector &vec, const bool apply_constraints = false)
		{
			//perf
			Chrono c;
			c.start();

			typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;


			static const int Backend = TraitsT::Backend;

			const auto &space = find_space<LibMeshFunctionSpace>(expr);
			const auto &dof_map = space.dof_map();
			auto &m = space.mesh();

			i_space.init(space);

			IntrepidFE fe;
			fe.init(1, i_space);

			ElementVector elem_vec;
			elem_vec.init(i_space, fe);

			QuadratureData quad_data;
			quad_data.init(i_space);
			quad_data.set(i_space, 1.);

			if(empty(vec) || size(vec).get(0) != dof_map.n_dofs() || !is_ghosted(vec)) {
				// vec = local_zeros(dof_map.n_local_dofs());
				vec = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list()); 
			} else {
				vec *= 0.;
			}

			FST::multiplyMeasure<Scalar>(fe.test_weighted, fe.dx, fe.fun);
			FST::integrate(elem_vec.vec, quad_data.data, fe.test_weighted);

			//assemble gradient
			auto local_vec = raw_type(vec)->getLocalView<Device>();
			int num_fields = fe.fe.getCardinality();

			auto elem_to_node = i_space.elem_to_node_.template view<Device>();
			// i_space.elem_to_node_.print();

			Kokkos::parallel_for(i_space.num_elements, KOKKOS_LAMBDA(const int k) {
			  for(int row = 0; row < num_fields; row++){
			      int row_index = elem_to_node(k, row);
			      Kokkos::atomic_fetch_add (&local_vec(row_index, 0),  elem_vec.vec(k, row) );
			  }
			});

			//perf
			c.stop();

			if(verbose_) {
				std::cout << "assemble: rhs" << std::endl;
				std::cout << c << std::endl;
			}

			return true;
		}

	private:
		IntrepidFunctionSpace i_space;
		IntrepidFE i_fe;
		bool verbose_;
	};
}

#endif //UTOPIA_INTREPID_2_ASSEMBLER_HPP
