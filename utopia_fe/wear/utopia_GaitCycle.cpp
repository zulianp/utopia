#include "utopia_GaitCycle.hpp"
#include "utopia_AffineTransform.hpp"

#include "utopia_ui.hpp"
#include "utopia_libmesh.hpp"

namespace utopia {
	void GaitCycle::read(InputStream &is)
	{
		is.read("rotation", [this](InputStream &is) {
			//TODO
			// is.read("block", rot.block);
			// is.read("axis",  rot.axis);
			// is.read("begin", rot.begin_angle_degree);
			// is.read("end",   rot.end_angle_degree);

			// rot.finalize();
		});
	}

	GaitCycle::GaitCycle()
	{
		n_time_steps = 50;
		t_end = 10.;
		angle_degree = 4.;
		start_angle_degree = -5.;
		init();
	}

	void GaitCycle::set_time_step(const std::size_t time_step)
	{
		t = time_step * dt;
		if(time_step > n_time_steps/2) {
			negative_dir = true;
		}
	}

	void GaitCycle::toggle_dir()
	{
		negative_dir = !negative_dir;
	}

	void GaitCycle::init()
	{
		t = 0.;
		dt = (t_end - t)/(n_time_steps - 1);
		angle_radian = (angle_degree/180 * M_PI);
		start_angle_radian = (start_angle_degree/180 * M_PI);
		d_angle = angle_radian/(t_end - t);
		negative_dir = false;

		rotate2 = [this](const Point2d &p) -> Point2d {
			AffineTransform trafo;
			trafo.make_rotation(2, this->t * this->d_angle, 'y');
			trafo.translation[0] = -p[0];
			trafo.translation[1] = -p[1];

			return trafo.apply(p);
		};

		rotate3 = [this](const Point3d &p) -> Point3d {
			AffineTransform trafo;
			trafo.make_rotation(3, this->start_angle_radian + this->t * this->d_angle, 'x');
			trafo.translation[0] = -p[0];
			trafo.translation[1] = -p[1] + 2.;
			trafo.translation[2] = -p[2] + this->t * 0.1;

			return trafo.apply(p);
		};

		zero2 = [](const Point2d &p) -> Point2d {
			return {0., 0.};
		};

		zero3  = [](const Point3d &p) -> Point3d {
			return {0., 0., 0.};
		};

		translate2_y = [this](const Point2d &p) -> Point2d {
			return {0., std::min(0.5 + 2*this->dt*0.1, 0.52) };
		};

		translate3_z = [this](const Point3d &p) -> Point3d {
			return {0., 0., 0.};
		};

		bc34 = [this](const Point3d &p) -> Point3d {
			return { 0., 0., 2.2 };
		};
	}

	void GaitCycle::override_displacement(
		const libMesh::MeshBase &mesh,
		const libMesh::DofMap &dof_map,
		const int block_id_rot,
		const int block_id_trasl,
		DVectord &displacement) const
	{
		//FIXME
		const int main_system_number = 0;
		// displacement = local_zeros(dof_map.n_local_dofs());

		const auto dim = mesh.mesh_dimension();
		const bool is_3d = dim == 3;

		auto r = range(displacement);
		Write<DVectord> w_d(displacement);

		for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
			const auto &e = **e_it;
			if(e.subdomain_id() != block_id_rot && e.subdomain_id() != block_id_trasl) continue;

			for(std::size_t i = 0; i < e.n_nodes(); ++i) {
				const auto &node = e.node_ref(i);

				Point2d p2{ node(0), node(1) };
				Point3d p3{ node(0), node(1), node(2) };

				if(e.subdomain_id() == block_id_rot){
					if(is_3d) {
						p3 = rotate3(p3);
					} else {
						p2 = rotate2(p2);
					}
				} else if(e.subdomain_id() == block_id_trasl) {
					if(is_3d) {
						p3 = translate3_z(p3);
					} else {
						p2 = translate2_y(p2);
					}
				}

				for(unsigned int d = 0; d < dim; ++d) {
					unsigned int dof = node.dof_number(main_system_number, d, 0);

					if(r.inside(dof)) {
						if(is_3d) {
							displacement.set(dof, p3[d]);
						} else {
							displacement.set(dof, p2[d]);
						}
					}
				}
			}
		}
	}
}