#include "utopia_ISO14243_3.hpp"

#include "utopia_AffineTransform.hpp"

#include "utopia_ui.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

namespace utopia {

	void ISO14243_3::read(InputStream &is)
	{
		flexion_extension_angle_axis_ = 'z';
		char axial_force_axis_char    = 'z';
		ap_motion_axis_ 			  = 'z';
		tibial_rotation_axis_         = 'z';

		std::string temp;
		is.read("flexion-extension-angle-axis", temp); flexion_extension_angle_axis_ = temp[0];
		is.read("axial-force-axis", temp); 		       axial_force_axis_char = temp[0];
		is.read("ap-motion-axis", temp); 			   ap_motion_axis_ = temp[0];
		is.read("tibial-rotation-axis", temp); 	       tibial_rotation_axis_ = temp[0];


		switch(axial_force_axis_char) {
			case 'x':
			case 'X':
			{
				axial_force_axis_ = 0; break;
			}

			case 'y':
			case 'Y':
			{
				axial_force_axis_ = 1; break;
			}

			case 'z':
			case 'Z':
			{
				axial_force_axis_ = 2; break;
			}

			default:
			{
				assert(false);
				break;
			}
		}


		femural_block_ = -1;
		tibial_block_ = -1;
		axial_force_side_ = -1;

		is.read("femural-block", femural_block_); assert(femural_block_ != -1);
		is.read("tibial-block", tibial_block_); assert(tibial_block_ != -1);
		is.read("axial-force-side", axial_force_side_); assert(axial_force_side_ != -1);

		dt_ = 0.1;
		is.read("dt", dt_);

		is.read("csv", temp);

		bool ok = read(temp); assert(ok);
	}
	
	void ISO14243_3::update(const int time_step) {
		assert(time_step >= 0);
		file_.get(time_step, 0, percentage_of_time_cycle_);
		file_.get(time_step, 1, flexion_extension_angle_);
		file_.get(time_step, 2, axial_force_);
		file_.get(time_step, 3, ap_motion_);
		file_.get(time_step, 4, tibial_int_ext_rotation_);


		std::cout << "time_step " << time_step << std::endl;
		describe_params(std::cout);
		//dt_?
	}

	void ISO14243_3::displacement_and_forces(
		ProductFunctionSpace<LibMeshFunctionSpace> &space,
		DVectord &displacement,
		DVectord &forces) const
	{
		const libMesh::MeshBase &mesh  = space[0].mesh();
		const libMesh::DofMap &dof_map = space[0].dof_map();
		auto dim = mesh.mesh_dimension();

		AffineTransform femoral_component;
		femoral_component.make_rotation(dim, flexion_extension_angle_ * (M_PI/180), flexion_extension_angle_axis_);

		AffineTransform tibial_component;
		{
			AffineTransform rot, tr;
			rot.make_rotation(dim, tibial_int_ext_rotation_ * (M_PI/180), tibial_rotation_axis_);
			tr.make_translation(dim, ap_motion_, ap_motion_axis_);
			tibial_component = tr * rot;
		}

		const int main_system_number = 0;

		if(empty(displacement)) {
			displacement = local_zeros(dof_map.n_local_dofs());
		}

		auto r = range(displacement);
		Write<DVectord> w_d(displacement);

		for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
			const auto &e = **e_it;
			// if(e.subdomain_id() != block_id_rot && e.subdomain_id() != block_id_trasl) continue;

			for(std::size_t i = 0; i < e.n_nodes(); ++i) {
				const auto &node = e.node_ref(i);

				Point3d p3{ node(0), node(1), node(2) };

				if(e.subdomain_id() == femural_block_) {
					p3 = femoral_component.apply(p3);
				} else if(e.subdomain_id() == tibial_block_) {
					p3 = tibial_component.apply(p3);
				}

				for(unsigned int d = 0; d < dim; ++d) {
					unsigned int dof = node.dof_number(main_system_number, d, 0);

					if(r.inside(dof)) {
						displacement.set(dof, p3[d] - node(d));
					}
				}
			}
		}


		forces = local_zeros(dof_map.n_local_dofs());
		auto v = test(space);


		auto l_form = surface_integral(
			inner(coeff(axial_force_), v[axial_force_axis_]), 
			axial_force_side_
		);

		utopia::assemble(l_form, forces);
	}

	void ISO14243_3::describe_params(std::ostream &os) const
	{
		os << "percentage_of_time_cycle: " << percentage_of_time_cycle_ << "\n";
		//degrees
		os << "flexion_extension_angle: " << flexion_extension_angle_ << "\n";

		//Newton
		os << "axial_force: " << axial_force_ << "\n";


		//mm
		os << "ap_motion: " << ap_motion_ << "\n";

		//degrees
		os << "tibial_int_ext_rotation: " << tibial_int_ext_rotation_ << "\n";

	}

	void ISO14243_3::describe(std::ostream &os) const
	{
		file_.write(os);

		os << "flexion_extension_angle_axis: " << flexion_extension_angle_axis_ << "\n";			
		os << "ap_motion_axis: " << ap_motion_axis_ << "\n";
		os << "tibial_rotation_axis: " << tibial_rotation_axis_ << "\n";
		os << "axial_force_axis: " << axial_force_axis_ << "\n";
		
		//blocks and side-sets
		os << "femural_block: " << femural_block_ << "\n";
		os << "tibial_block: " << tibial_block_ << "\n";
		os << "axial_force_side: " << axial_force_side_ << "\n";

		os << "dt: " << dt_ << "\n";

		os << "n_steps: " << n_steps()  << "\n";


		describe_params(os);
		
	}
}
