#include "utopia_GaitCycle.hpp"
#include "utopia_AffineTransform.hpp"

#include "utopia_ui.hpp"
#include "utopia_libmesh.hpp"

namespace utopia {
	void GaitCycle::read(InputStream &is)
	{
		t = 0.;
		t_end = 1.;
		n_time_steps = 1;

		is.read("t0", t);
		is.read("t_end", t_end);
		is.read("steps", n_time_steps);

		dt = (t_end - t)/n_time_steps;

		is.read("rotations", [this](InputStream &is) {
			is.read_all([this](InputStream &is) {
				Rotation rot;

				std::string axis;

				is.read("block", rot.block);
				is.read("axis",  axis);
				is.read("begin", rot.begin_angle_degree);
				is.read("end",   rot.end_angle_degree);

				rot.axis = axis[0];

				rotations.push_back(rot);
			});
		});

		is.read("translations", [this](InputStream &is) {
			is.read_all([this](InputStream &is) {
				Translation tr;

				std::string axis;

				is.read("block", tr.block);
				is.read("axis",  axis);
				is.read("begin", tr.begin_offset);
				is.read("end",   tr.end_offset);

				tr.axis = axis[0];

				translations.push_back(tr);
			});
		});
	}

	GaitCycle::GaitCycle()
	{
		n_time_steps = 50;
		t_end = 10.;
	}

	void GaitCycle::set_time_step(const std::size_t time_step)
	{
		t = time_step * dt;
		// if(time_step > n_time_steps/2) {
		// 	negative_dir = true;
		// }

		update_trafos(time_step, t);
	}

	void GaitCycle::update_trafos(const int step, const double t)
	{
		for(auto &r : rotations) {
			r.update(step, t);
		}

		for(auto &tr : translations) {
			tr.update(step, t);
		}
	}

	void GaitCycle::init_trafos(const int n_dims, const int n_steps)
	{
		for(auto &r : rotations) {
			r.init(n_dims, n_steps);
		}

		for(auto &tr : translations) {
			tr.init(n_dims, n_steps);
		}
	}

	void GaitCycle::toggle_dir()
	{
		// negative_dir = !negative_dir;
	}

	void GaitCycle::init(int n_dims)
	{
		init_trafos(n_dims, n_time_steps);
	}

	void GaitCycle::override_displacement(
		const libMesh::MeshBase &mesh,
		const libMesh::DofMap &dof_map,
		DVectord &displacement) const
	{

		const int main_system_number = 0;

		if(empty(displacement)) {
			displacement = local_zeros(dof_map.n_local_dofs());
		}

		const auto dim = mesh.mesh_dimension();
		const bool is_3d = dim == 3;



		auto r = range(displacement);
		Write<DVectord> w_d(displacement);

		for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
			const auto &e = **e_it;
			// if(e.subdomain_id() != block_id_rot && e.subdomain_id() != block_id_trasl) continue;

			for(std::size_t i = 0; i < e.n_nodes(); ++i) {
				const auto &node = e.node_ref(i);

				Point3d p3{ node(0), node(1), node(2) };

				for(auto &r : rotations) {
					if(e.subdomain_id() == r.block) {
						p3 = r.trafo.apply(p3);
					}
				}

				for(auto &t : translations) {
					if(e.subdomain_id() == t.block) {
						p3 = t.trafo.apply(p3);
					}
				}

				for(unsigned int d = 0; d < dim; ++d) {
					unsigned int dof = node.dof_number(main_system_number, d, 0);

					if(r.inside(dof)) {
						displacement.set(dof, p3[d] - node(d));
					}
				}
			}
		}
	}

	GaitCycle::Rotation::Rotation()
	{
		block = 1;
		n_dims = 3;
		axis = 'z';
		begin_angle_degree = 0;
		end_angle_degree = 90;
		d_angle = end_angle_degree - begin_angle_degree;
	}

	void GaitCycle::Rotation::init(
		const int n_dims,
		const int n_steps)
	{
		double range = (end_angle_degree - begin_angle_degree) * (M_PI/180.);
		d_angle = range/n_steps;
	}

	void GaitCycle::Rotation::update(
		const int step,
		const double t)
	{
		auto begin_angle_radian = begin_angle_degree * (M_PI/180.);
		trafo.make_rotation(n_dims, begin_angle_radian + step * d_angle, axis);
	}


	GaitCycle::Translation::Translation()
	{
		block = 1;
		n_dims = 3;
		axis = 'z';
		begin_offset = 0;
		end_offset   = 1;
		d_offset = end_offset - begin_offset;
	}

	void GaitCycle::Translation::init(
		const int n_dims,
		const int n_steps)
	{
		double range = (end_offset - begin_offset);
		d_offset = range/n_steps;
	}

	void GaitCycle::Translation::update(
		const int step,
		const double t
		)
	{
		trafo.make_translation(n_dims, begin_offset + step * d_offset, axis);
	}
}
