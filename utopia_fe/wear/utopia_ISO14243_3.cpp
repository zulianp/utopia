#include "utopia_libmesh.hpp"
#include "utopia_ISO14243_3.hpp"

#include "utopia_AffineTransform.hpp"

#include "utopia_ui.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"

namespace utopia {

    static int char_to_dim(const char dim)
    {
        switch(dim) {
            case '1':
            case 'x':
            case 'X':
            {
                return 0;
            }
            case '2':
            case 'y':
            case 'Y':
            {
                return 1;
            }
            case '3':
            case 'z':
            case 'Z':
            {
                return 2;
            }
            default:
            {
                break;
            }
        }

        assert(false);
        return -1;
    }

    void ISO14243_3::read(Input &is)
    {
        flexion_extension_angle_axis_ = 'z';
        char axial_force_axis_char    = 'z';
        ap_motion_axis_ 			  = 'z';
        tibial_rotation_axis_         = 'z';
        normalize_axial_force_by_area_ = 0;

        std::string temp;
        is.get("flexion-extension-angle-axis", temp); flexion_extension_angle_axis_ = temp[0];
        is.get("axial-force-axis", temp); 		       axial_force_axis_char = temp[0];
        is.get("ap-motion-axis", temp); 			   ap_motion_axis_ = temp[0];
        is.get("tibial-rotation-axis", temp); 	       tibial_rotation_axis_ = temp[0];
        is.get("normalize-axial-force-by-area", normalize_axial_force_by_area_);

        axial_force_axis_ = char_to_dim(axial_force_axis_char);


        femural_block_ = -1;
        tibial_block_ = -1;
        axial_force_side_ = -1;

        is.get("femural-block", femural_block_); assert(femural_block_ != -1);
        is.get("tibial-block", tibial_block_); assert(tibial_block_ != -1);
        is.get("axial-force-side", axial_force_side_); assert(axial_force_side_ != -1);
        is.get("tibial-rotation-offset", tibial_rotation_offset_);

        is.get("tibial-rotation-offset-axis", temp);
        tibial_rotation_offset_axis_ = char_to_dim(temp[0]);


        dt_ = 0.1;
        is.get("dt", dt_);

        is.get("csv", temp);

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
        UVector &displacement,
        UVector &forces) const
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


        CompositeAffineTransform composite_tibial_component;
        composite_tibial_component.compose_left(tibial_component);

        if(tibial_rotation_offset_ != 0.)
        {
            AffineTransform tibial_offset;
            tibial_offset.make_translation(dim, tibial_rotation_offset_, tibial_rotation_offset_axis_);
            composite_tibial_component.compose_right(tibial_offset);

            tibial_offset.make_translation(dim, -tibial_rotation_offset_, tibial_rotation_offset_axis_);
            composite_tibial_component.compose_left(tibial_offset);
        }




        const int main_system_number = 0;

        if(empty(displacement)) {
            displacement = local_zeros(dof_map.n_local_dofs());
        }

        auto r = range(displacement);
        Write<UVector> w_d(displacement);

        for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
            const auto &e = **e_it;
            // if(e.subdomain_id() != block_id_rot && e.subdomain_id() != block_id_trasl) continue;

            for(std::size_t i = 0; i < e.n_nodes(); ++i) {
                const auto &node = e.node_ref(i);

                Point3d p3{ node(0), node(1), node(2) };

                if(e.subdomain_id() == femural_block_) {
                    p3 = femoral_component.apply(p3);
                } else if(e.subdomain_id() == tibial_block_) {
                    p3 = composite_tibial_component.apply(p3);
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

        double area = 1.;
        if(normalize_axial_force_by_area_) {
            area = surface_area(space[0], axial_force_side_);
            std::cout << "normalizing axial force by area: " << area << std::endl;
        }

        auto l_form = surface_integral(
            (1./area) * inner(coeff(axial_force_), v[axial_force_axis_]),
            axial_force_side_
        );

        utopia::assemble(l_form, forces);

        double sum_axial_force = sum(forces);
        std::cout << "sum(axial-force): " << sum_axial_force << std::endl;

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

        os << "normalize_axial_force_by_area: " << normalize_axial_force_by_area_ << "\n";


        os << "tibial-rotation-offset: " << tibial_rotation_offset_ << "\n";
        os << "tibial-rotation-offset-axis: " << tibial_rotation_offset_axis_ << "\n";

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
