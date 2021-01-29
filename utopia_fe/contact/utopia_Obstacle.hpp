#ifndef UTOPIA_OBSTACLE_HPP
#define UTOPIA_OBSTACLE_HPP

#include "utopia_ui.hpp"

#include "utopia_IContact.hpp"
#include "utopia_LibMeshFunctionSpaceAdapter.hpp"
#include "utopia_NodeBlackLister.hpp"

namespace utopia {
    class ObstacleParams : public Configurable {
    public:
        void read(Input &in) override;

        int variable_number{0};
        double max_gap{0.1};
    };

    class Obstacle {
    public:
        bool assemble(libMesh::MeshBase &mesh,
                      libMesh::DofMap &dof_map,
                      const ObstacleParams &params = ObstacleParams());

        bool init(libMesh::MeshBase &obstacle_mesh);
        void transform(const USparseMatrix &in, USparseMatrix &out);
        void transform(const UVector &in, UVector &out);
        void inverse_transform(const UVector &in, UVector &out);

        class Impl;
        Obstacle();
        ~Obstacle();

        class Output {
        public:
            UVector gap;
            UVector normals;
            USparseMatrix orthogonal_trafo;

            USparseMatrix basis_trafo;
            UVector mass_vector;
            UVector inverse_mass_vector;
            UVector is_contact;
        };

        inline Output &output() { return output_; }

    private:
        std::unique_ptr<Impl> impl_;
        Output output_;
    };
}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_HPP
