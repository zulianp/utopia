#ifndef UTOPIA_BOUNDING_BOX_COARSENER_HPP
#define UTOPIA_BOUNDING_BOX_COARSENER_HPP

#include <memory>
#include <iostream>
#include <vector>

namespace libMesh {
    class MeshBase;
    class UnstructuredMesh;
}

namespace utopia {
    class BoundingBoxCoarsener {
    public:
        class Impl;

        BoundingBoxCoarsener();
        ~BoundingBoxCoarsener();
        bool init(const int n_coarsening_levels, const libMesh::MeshBase &mesh);
        void describe(std::ostream &os = std::cout) const;

        const std::shared_ptr<libMesh::UnstructuredMesh> &get_mesh() const;
        const std::vector< std::pair<int, int> > &get_tags() const;


        std::shared_ptr<Impl> impl_;
    };
}

#endif //UTOPIA_BOUNDING_BOX_COARSENER_HPP
