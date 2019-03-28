#include "utopia_libmesh.hpp"
#include "utopia_OldConfiguration.hpp"

#include "utopia_libmesh_FunctionSpace.hpp"

namespace utopia {

    OldConfiguration::OldConfiguration()
    {
        n_time_steps = 50;
        t_end = 10.;
    }

    void OldConfiguration::read(Input &is)
    {
        t = 0.;
        t_end = 1.;
        n_time_steps = 1;

        is.get("t0", t);
        is.get("t_end", t_end);
        is.get("steps", n_time_steps);

        dt_ = (t_end - t)/n_time_steps;

        is.get("rotations", [this](Input &is) {
            is.get_all([this](Input &is) {
                Rotation rot;

                std::string axis;

                is.get("block", rot.block);
                is.get("axis",  axis);
                is.get("begin", rot.begin_angle_degree);
                is.get("end",   rot.end_angle_degree);
                is.get("from",  rot.from_step);
                is.get("to",    rot.to_step);

                rot.axis = axis[0];

                rotations.push_back(rot);
            });
        });

        is.get("translations", [this](Input &is) {
            is.get_all([this](Input &is) {
                Translation tr;

                std::string axis;

                is.get("block", tr.block);
                is.get("axis",  axis);
                is.get("begin", tr.begin_offset);
                is.get("end",   tr.end_offset);

                is.get("from",  tr.from_step);
                is.get("to",    tr.to_step);

                tr.axis = axis[0];

                translations.push_back(tr);
            });
        });
    }

    double OldConfiguration::dt() const
    {
        return dt_;
    }

    void OldConfiguration::update(const int time_step)
    {
        t = time_step * dt_;
        for(auto &r : rotations) {
            r.update(time_step, t);
        }

        for(auto &tr : translations) {
            tr.update(time_step, t);
        }
    }

    void OldConfiguration::displacement_and_forces(
        ProductFunctionSpace<LibMeshFunctionSpace> &space,
        UVector &displacement,
        UVector &forces) const
    {
        const libMesh::MeshBase &mesh  = space[0].mesh();
        const libMesh::DofMap &dof_map = space[0].dof_map();
        auto dim = mesh.mesh_dimension();

        const int main_system_number = 0;

        if(empty(displacement)) {
            displacement = local_zeros(dof_map.n_local_dofs());
        }

        if(empty(forces)) {
            forces = local_zeros(dof_map.n_local_dofs());
        }

        auto r = range(displacement);
        Write<UVector> w_d(displacement);

        for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
            const auto &e = **e_it;
            // if(e.subdomain_id() != block_id_rot && e.subdomain_id() != block_id_trasl) continue;

            for(std::size_t i = 0; i < e.n_nodes(); ++i) {
                const auto &node = e.node_ref(i);

                Point3d p3{ node(0), node(1), node(2) };

                for(auto &r : rotations) {
                    if(!r.active) continue;

                    if(e.subdomain_id() == r.block) {
                        p3 = r.trafo.apply(p3);
                    }
                }

                for(auto &t : translations) {
                    if(!t.active) continue;

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

    void OldConfiguration::init(ProductFunctionSpace<LibMeshFunctionSpace> &space)
    {
        auto dim = space[0].mesh().mesh_dimension();

        for(auto &r : rotations) {
            r.init(dim, n_time_steps);
        }

        for(auto &tr : translations) {
            tr.init(dim, n_time_steps);
        }
    }


    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////


    OldConfiguration::Rotation::Rotation()
    {
        block = 1;
        n_dims = 3;
        axis = 'z';
        begin_angle_degree = 0;
        end_angle_degree = 90;
        d_angle = end_angle_degree - begin_angle_degree;

        from_step = -1;
        to_step   = -1;

        active = true;
    }

    void OldConfiguration::Rotation::init(
        const int n_dims,
        const int n_steps)
    {
        double range = (end_angle_degree - begin_angle_degree) * (M_PI/180.);

        if(from_step != -1 && to_step != -1) {
            d_angle = range/(to_step - from_step);
        } else {
            d_angle = range/n_steps;
        }
    }

    void OldConfiguration::Rotation::update(
        const int step,
        const double t)
    {
        auto begin_angle_radian = begin_angle_degree * (M_PI/180.);

        if(from_step == -1 || to_step == -1) {
            trafo.make_rotation(n_dims, begin_angle_radian + step * d_angle, axis);
            return;
        }

        active = from_step <= step && step < to_step;
        trafo.make_rotation(n_dims, begin_angle_radian + (step - from_step) * d_angle, axis);
    }


    ///////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////

    OldConfiguration::Translation::Translation()
    {
        block = 1;
        n_dims = 3;
        axis = 'z';
        begin_offset = 0;
        end_offset   = 1;
        d_offset = end_offset - begin_offset;

        from_step = -1;
        to_step   = -1;

        active = true;
    }

    void OldConfiguration::Translation::init(
        const int n_dims,
        const int n_steps)
    {
        double range = (end_offset - begin_offset);

        if(from_step != -1 && to_step != -1) {
            d_offset = range/(to_step - from_step);
        } else {
            d_offset = range/n_steps;
        }
    }

    void OldConfiguration::Translation::update(
        const int step,
        const double t
        )
    {
        if(from_step == -1 || to_step == -1) {
            trafo.make_translation(n_dims, begin_offset + step * d_offset, axis);
            return;
        }

        active = from_step <= step && step < to_step;
        trafo.make_translation(n_dims, begin_offset + (step - from_step) * d_offset, axis);
    }


}
