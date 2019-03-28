#include "utopia_SurfUtils.hpp"
#include "libmesh/quadrature_gauss.h"

namespace utopia {

    void SurfUtils::avg_normal(
        const libMesh::Elem &trial,
        const libMesh::FEType &type,
        libMesh::Point &normal)
    {
        auto fe = libMesh::FEBase::build(trial.dim(), type);
        auto &n = fe->get_normals();

        libMesh::QGauss q(trial.dim(), type.order);
        q.init(trial.type());

        fe->attach_quadrature_rule(&q);
        fe->reinit(&trial);

        normal.zero();

        for(const auto &n_i : n) {
            normal += n_i;
        }

        normal /= n.size();
    }

}
