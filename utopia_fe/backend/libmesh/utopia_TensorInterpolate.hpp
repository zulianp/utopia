#ifndef UTOPIA_TENSOR_INTERPOLATE_HPP
#define UTOPIA_TENSOR_INTERPOLATE_HPP

#include "utopia_Traits.hpp"
#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_Path.hpp"
#include "libmesh/mesh.h"

#include <fstream>
#include <string>

namespace utopia {

    template<class FunctionSpace, class Vector>
    class TensorInterpolate : public Expression< TensorInterpolate<FunctionSpace, Vector> > {
    public:
        typedef typename Traits<Vector>::Scalar Scalar;
        typedef typename Traits<Vector>::SizeType SizeType;
        static const int Order = 2;

        TensorInterpolate(
            const FunctionSpace &space,
            const Vector &vec,
            const SizeType rows,
            const SizeType cols)
        : space_(space), vec_(vec), rows_(rows), cols_(cols)
        {}

        const FunctionSpace &space() const
        {
            return space_;
        }

        const Vector &vec() const
        {
            return vec_;
        }

        const SizeType &rows() const { return rows_; }
        const SizeType &cols() const { return cols_; }

    private:
        const FunctionSpace &space_;
        const Vector &vec_;
        const SizeType rows_;
        const SizeType cols_;
    };

    template<class FunctionSpace, class Vector>
    class Traits< TensorInterpolate<FunctionSpace, Vector> > : public Traits<FunctionSpace> {
    public:
        static const int FILL_TYPE = FillType::DENSE;
    };

    template<class FunctionSpace, class Vector>
    inline TensorInterpolate<FunctionSpace, Vector> tensor_interpolate(
        const FunctionSpace &space,
        const Vector &vec,
        const int rows,
        const int cols)
    {
        return TensorInterpolate<FunctionSpace, Vector>(space, vec, rows, cols);
    }

    template<class FunctionSpace, class Vector, class Visitor>
    inline static int traverse(const TensorInterpolate<FunctionSpace, Vector>  &expr, Visitor &visitor)
    {
        return visitor.visit(expr);
    }

    template<class FunctionSpace, class Vector, class Traits, int Backend, int IsQuadData>
    class FEEval<TensorInterpolate<FunctionSpace, Vector>, Traits, Backend, IsQuadData> {
    public:
        using Expr = utopia::TensorInterpolate<FunctionSpace, Vector>;

        inline static auto apply(
            const Expr &expr,
            AssemblyContext<Backend> &ctx) -> std::vector<USerialMatrix>
        {
            const auto &s        = expr.space();
            const auto &mesh     = s.mesh();
            const auto &elem_ptr = mesh.elem(ctx.current_element());
            const auto &dof_map  = s.dof_map();
            const auto n_entries = expr.rows() * expr.cols();

            std::vector<libMesh::dof_id_type> scalar_indices, indices;
            dof_map.dof_indices(elem_ptr, scalar_indices, s.subspace_id());

            const std::size_t n_dofs = scalar_indices.size();
            indices.resize(n_dofs * n_entries);

            std::size_t index = 0;
            for(std::size_t i = 0; i < n_dofs; ++i) {
                for(std::size_t j = 0; j < n_entries; ++j, ++index) {
                    indices[index] = scalar_indices[i] * n_entries + j;
                }
            }

            USerialVector element_values = zeros(indices.size());
            Read<UVector> r(expr.vec());

            typename utopia::Traits<UVector>::IndexSet u_index;
            u_index.insert(u_index.end(), indices.begin(), indices.end());

            assert( expr.vec().has_ghosts() || expr.vec().comm().size() == 1);
            expr.vec().get(u_index, element_values.entries());

            ////////////////////////////////////////
            const auto &f = ctx.fe()[0]->get_phi();

            const std::size_t n_fun = f.size();
            const std::size_t n_qp  = f[0].size();

            assert(n_fun == n_dofs);

            const std::size_t rows = expr.rows();
            const std::size_t cols = expr.cols();

            std::vector<USerialMatrix> tensors(n_fun, zeros(rows, cols));

            index = 0;
            for(std::size_t l = 0; l < n_fun; ++l) {
                for(std::size_t i = 0; i < rows; ++i) {
                    for(std::size_t j = 0; j < cols; ++j, ++index) {
                        tensors[l].set(i, j, element_values.get(index));
                    }
                }
            }

            std::vector<USerialMatrix> ret(n_qp, zeros(rows, cols));

            for(std::size_t k = 0; k < n_qp; ++k) {
                auto &mat = ret[k];
                for(std::size_t i = 0; i < n_fun; ++i) {
                    mat += f[i][k] * tensors[i];
                }
            }

            return ret;
        }
    };

    inline void build_tensor_ghosted(
        const Traits<UVector>::SizeType &n_local_dofs,
        const Traits<UVector>::SizeType &n_dofs,
        const Traits<UVector>::SizeType &dims,
        const UIndexSet &ghost_nodes,
        UVector &vec)
    {
        using SizeType = Traits<UVector>::SizeType;

        if(vec.comm().size() == 1) {
            if(empty(vec)) {
                vec = local_zeros(n_local_dofs * dims);
            }

            return;
        }

        UVector temp;
        if(!empty(vec)) {
            temp = vec;
        }

        UIndexSet tensor_ghost_nodes(ghost_nodes.size() * dims);

        const SizeType n_ghosts = ghost_nodes.size();

        for(SizeType i = 0, index = 0; i < n_ghosts; ++i) {
            for(SizeType d = 0; d < dims; ++d) {
                tensor_ghost_nodes[index++] = ghost_nodes[i] * dims + d;
            }
        }

        vec = ghosted(n_local_dofs * dims, n_dofs * dims, tensor_ghost_nodes);
        if(!empty(temp)) {
            //copy ghosted entries
            vec = temp;
            synchronize(vec);
        }
    }

    inline bool read_tensor_data(
        const utopia::Path &path,
        const USizeType &n_local_nodes,
        const USizeType &n,
        UVector &vec)
    {
        if(empty(vec) || !vec.comm().conjunction(n_local_nodes * n == vec.local_size())) {
            vec = local_zeros(n_local_nodes * n);
        }

        std::ifstream is(path.c_str());
        if(!is.good()) { return false; }

        auto r = range(vec);
        Write<UVector> w(vec);

        UScalar val = 0.0;
        SizeType i = 0;
        while(is.good()) {
            is >> val;

            if(i >= r.begin()) {
                vec.set(i, val);
            }

            i++;

            if(i >= r.end()) {
                break;
            }
        }

        is.close();
        return true;
    }
}


#endif //UTOPIA_TENSOR_INTERPOLATE_HPP
