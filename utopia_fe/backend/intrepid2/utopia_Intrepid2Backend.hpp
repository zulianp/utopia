#ifndef UTOPIA_INTREPID2_BACKEND_HPP
#define UTOPIA_INTREPID2_BACKEND_HPP

#include "utopia_fe_kokkos_fix.hpp"

//#include "utopia_LibMeshLambdaAssembly.hpp"
#include "utopia_fe_core.hpp"

// Define the Finite Element object.
/*#include "libmesh/fe.h"

// Define Gauss quadrature rules.
#include "libmesh/quadrature_gauss.h"*/

// Define the DofMap, which handles degree of freedom
// indexing.
#include "libmesh/dof_map.h"

// Define useful datatypes for finite element
// matrix and vector components.

// Define the PerfLog, a performance logging utility.
// It is useful for timing events in a code and giving
// you an idea where bottlenecks lie.
/*#include "libmesh/perf_log.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// To impose Dirichlet boundary conditions
#include "libmesh/analytic_function.h"
#include "libmesh/dirichlet_boundaries.h"

#include "libmesh/string_to_enum.h"
// #include "libmesh/getpot.h"

#include "libmesh/const_function.h"

#include "utopia.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_LambdaFunction.hpp"*/
#include "utopia_libmesh_Types.hpp"

namespace utopia {

    void apply_boundary_conditions(libMesh::DofMap &dof_map, TSUSerialMatrix &mat, TUSerialVector &vec);
    void apply_boundary_conditions(libMesh::DofMap &dof_map, TUSerialVector &vec);

    template <class Vector>
    void mark_constrained_dofs(libMesh::DofMap &dof_map, Wrapper<Vector, 1> &vec) {
        vec = local_zeros(dof_map.n_local_dofs());
        const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

        Write<Wrapper<Vector, 1>> w_v(vec);

        if (has_constaints) {
            Range r = range(vec);
            for (SizeType i = r.begin(); i < r.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    auto value = 1.;
                    vec.set(i, value);
                }
            }
        }
    }

    template <class FEFunction, class Matrix, class Vector>
    void apply_boundary_conditions(FEFunction &u, Matrix &mat, Vector &vec) {
        apply_boundary_conditions(u.dof_map(), mat, vec);
    }

    template <class DofMap, class Vector>
    void apply_zero_boundary_conditions(DofMap &dof_map, Vector &vec) {
        const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

        Write<Vector> w_v(vec);

        if (has_constaints) {
            Range r = range(vec);
            for (SizeType i = r.begin(); i < r.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    vec.set(i, 0.0);
                }
            }
        }
    }
    /*
        template<class DofMap, class Matrix>
        void set_identity_at_constraint_rows(DofMap &dof_map, Matrix &mat)
        {
            bool has_constaints = true;
            if(dof_map.constraint_rows_begin() == dof_map.constraint_rows_end()) {
                has_constaints = false;
            }

            using SizeT = UTOPIA_SIZE_TYPE(Matrix);

            std::vector<SizeT> rows;
            rows.reserve(local_size(mat).get(0));

            auto rr = row_range(mat);

            if(has_constaints) {
                for(auto i = rr.begin(); i < rr.end(); ++i) {
                    dof_map.is_constrained_dof(i);
                    rows.push_back(i);
                }
            }

            set_zero_rows(mat, rows, 1.);
        }

        template<class DofMap, class Matrix>
        void set_zero_at_constraint_rows(DofMap &dof_map, Matrix &mat)
        {
            bool has_constaints = true;
            if( dof_map.constraint_rows_begin() == dof_map.constraint_rows_end()) {
                has_constaints = false;
            }

            Size s = size(mat);
            Matrix temp = mat;

            {
                Write<Matrix> w_t(mat);

                each_read(temp, [&](const SizeType i, const SizeType j, const libMesh::Real value) {
                    if(has_constaints && dof_map.is_constrained_dof(i)) {
                        mat.set(i, j, 0.0);
                    }
                });
            }
        }

        inline void convert(libMesh::NumericVector<libMesh::Number> &lm_vec, PetscVector &utopia_vec)
        {
            using namespace libMesh;
            Vec p_vec = cast_ptr< libMesh::PetscVector<libMesh::Number> *>(&lm_vec)->vec();
            utopia::convert(p_vec, utopia_vec);
        }



        inline void convert(libMesh::SparseMatrix<libMesh::Number> &lm_mat, PetscMatrix &utopia_mat) {
            using namespace libMesh;

            Mat p_mat = cast_ptr< libMesh::PetscMatrix<libMesh::Number> *>(&lm_mat)->mat();
            utopia::convert(p_mat, utopia_mat);
        }

    #ifdef WITH_TRILINOS
        inline void convert(libMesh::NumericVector<libMesh::Number> &lm_vec, TUSerialVector &utopia_vec)
        {
            //FIXME inefficient
            PetscVector temp;
            utopia::convert(lm_vec, temp);
            utopia::backend_convert(temp, utopia_vec);
        }


        inline void convert(libMesh::SparseMatrix<libMesh::Number> &lm_mat, TSUSerialMatrix &utopia_mat) {
            using namespace libMesh;

            Mat p_mat = cast_ptr< libMesh::PetscMatrix<libMesh::Number> *>(&lm_mat)->mat();

            //FIXME inefficient
            PetscMatrix temp;
            utopia::convert(p_mat, temp);
            backend_convert(temp, utopia_mat);
        }

    #endif //WITH_TRILINOS


        inline void convert(USparseMatrix &utopia_mat, libMesh::SparseMatrix<libMesh::Number> &lm_mat) {
            using namespace libMesh;
            using namespace utopia;


            //Does not work
            // Mat p_mat = cast_ptr< libMesh::PetscMatrix<libMesh::Number> *>(&lm_mat)->mat();
            // utopia::convert(utopia_mat, p_mat);

            each_read(utopia_mat, [&lm_mat](const SizeType i, const SizeType j, const double value) -> void {
                lm_mat.set(i, j, value);
            });
        }

        inline void convert(UVector &utopia_vec, libMesh::NumericVector<libMesh::Number> &lm_vec)
        {
            {
                Read<UVector> w_s(utopia_vec);
                Range r = range(utopia_vec);
                for(long i = r.begin() ; i < r.end(); ++i) {
                    lm_vec.set(i, utopia_vec.get(i) );
                }
            }
        }

        template<class Space>
        inline static bool has_constrained_dofs(const Space &space, const libMesh::Elem &elem)
        {
            std::vector<libMesh::dof_id_type> indices;
            space.dof_map().dof_indices(&elem, indices, space.var_num());

            for(auto i : indices) {
                if(space.dof_map().is_constrained_dof(i)) {
                    return true;
                }
            }

            return false;
        }


        inline std::ostream &logger()
        {
            return std::cout;
        }



        inline bool is_symmetric(const libMesh::DenseMatrix<libMesh::Real> &mat)
        {
            if(mat.m() != mat.n()) return false;

            for(int i = 0; i < mat.m(); ++i) {
                for(int j = i+1; j < mat.n(); ++j) {
                    if(fabs(mat(i, j) - mat(j, i)) > 1e-14) {
                        return false;
                    }
                }
            }

            return true;
        }

        inline bool is_approx_equal(const libMesh::DenseMatrix<libMesh::Real> &left, const
    libMesh::DenseMatrix<libMesh::Real> &right)
        {
            if(left.m() != right.m()) return false;
            if(left.n() != right.n()) return false;

            for(int i = 0; i < left.m(); ++i) {
                for(int j = i+1; j < left.n(); ++j) {
                    if(fabs(left(i, j) - right(i, j)) > 1e-14) {
                        return false;
                    }
                }
            }

            return true;
        }

        inline bool is_approx_equal_tr(const libMesh::DenseMatrix<libMesh::Real> &left, const
    libMesh::DenseMatrix<libMesh::Real> &right)
        {
            if(left.m() != right.n()) return false;
            if(left.n() != right.m()) return false;

            for(int i = 0; i < left.m(); ++i) {
                for(int j = i+1; j < left.n(); ++j) {
                    if(fabs(left(i, j) - right(j, i)) > 1e-14) {
                        return false;
                    }
                }
            }

            return true;
        }

        inline libMesh::MeshTools::BoundingBox bounding_box(const libMesh::MeshBase &mesh) {
    #if LIBMESH_VERSION_LESS_THAN(1, 3, 0)
                    libMesh::MeshTools::BoundingBox bb = libMesh::MeshTools::bounding_box(mesh);
    #else
                    libMesh::MeshTools::BoundingBox bb = libMesh::MeshTools::create_bounding_box(mesh);
    #endif

            return bb;
        }*/
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_BACKEND_HPP
