#ifndef UTOPIA_VOXEL_TO_ELEMENT_HPP
#define UTOPIA_VOXEL_TO_ELEMENT_HPP

namespace utopia {

    class Voxel2Element {
    public:
        inline static std::unique_ptr<libMesh::SerialMesh> build(
            const libMesh::Parallel::Communicator &comm,
            const Grid<2> &grid,
            Grid<2>::Integer cell_index)
        {
            auto index_min = grid.element_index(cell_index);
            auto index_max = index_min;

            for(auto &i : index_max) {
                ++i;
            }

            auto p_min = grid.point(index_min);
            auto p_max = grid.point(index_max);

            auto mesh = make_unique<libMesh::SerialMesh>(comm, 2);
            mesh->reserve_nodes(4);

            libMesh::Point p;
            //point 0
            p(0) = p_min.x;
            p(1) = p_min.y;
            mesh->add_point(p);

            //point 1
            p(0) = p_max.x;
            p(1) = p_min.y;
            mesh->add_point(p);

            //point 2
            p(0) = p_max.x;
            p(1) = p_max.y;
            mesh->add_point(p);

            //point 3
            p(0) = p_min.x;
            p(1) = p_max.y;
            mesh->add_point(p);

            auto elem = libMesh::Elem::build(libMesh::QUAD4);

            for (int i = 0; i < 4; ++i) {
                elem->set_node(i) = & mesh->node(i);
            }


            mesh->add_elem(elem.release());
            return std::move(mesh);
        }

        inline static std::unique_ptr<libMesh::SerialMesh> build(
            const libMesh::Parallel::Communicator &comm,
            const Grid<3> &grid,
            Grid<3>::Integer cell_index)
        {

            auto index_min = grid.element_index(cell_index);
            auto index_max = index_min;

            for(auto &i : index_max) {
                ++i;
            }

            auto p_min = grid.point(index_min);
            auto p_max = grid.point(index_max);

                /*
                * HEX8:   7        6
                *         o--------o
                *        /:       /|
                *       / :      / |
                *    4 /  :   5 /  |
                *     o--------o   |
                *     |   o....|...o 2
                *     |  .3    |  /
                *     | .      | /
                *     |.       |/
                *     o--------o
                *     0        1
                */

            auto mesh = make_unique<libMesh::SerialMesh>(comm, 3);
            mesh->reserve_nodes(8);

            libMesh::Point p;
            //point 0
            p(0) = p_min.x;
            p(1) = p_min.y;
            p(2) = p_min.z;
            mesh->add_point(p);

            //point 1
            p(0) = p_max.x;
            p(1) = p_min.y;
            p(2) = p_min.z;
            mesh->add_point(p);

            //point 2
            p(0) = p_max.x;
            p(1) = p_max.y;
            p(2) = p_min.z;
            mesh->add_point(p);

            //point 3
            p(0) = p_min.x;
            p(1) = p_max.y;
            p(2) = p_min.z;
            mesh->add_point(p);

            //point 4
            p(0) = p_min.x;
            p(1) = p_min.y;
            p(2) = p_max.z;
            mesh->add_point(p);

            //point 5
            p(0) = p_max.x;
            p(1) = p_min.y;
            p(2) = p_max.z;
            mesh->add_point(p);

            //point 6
            p(0) = p_max.x;
            p(1) = p_max.y;
            p(2) = p_max.z;
            mesh->add_point(p);

            //point 7
            p(0) = p_min.x;
            p(1) = p_max.y;
            p(2) = p_max.z;
            mesh->add_point(p);

            auto elem = libMesh::Elem::build(libMesh::HEX8);

            for (int i = 0; i < 8; ++i) {
                elem->set_node(i) = & mesh->node(i);
            }

            mesh->add_elem(elem.release());
            return std::move(mesh);
        }
    };

}

#endif //UTOPIA_VOXEL_TO_ELEMENT_HPP
