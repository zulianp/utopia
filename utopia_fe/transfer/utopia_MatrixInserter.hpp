#ifndef UTOPIA_MATRIX_INSERTER_HPP
#define UTOPIA_MATRIX_INSERTER_HPP

#include "moonolith_communicator.hpp"
#include "moonolith_sparse_matrix.hpp"
#include "moonolith_redistribute.hpp"

#include "utopia_fe_base.hpp"

#include <vector>

namespace utopia {

    //matrix proxy for utopia
    class MatrixInserter {
    public:
        MatrixInserter(MPI_Comm mpi_comm, const bool use_add = true) :
          comm(mpi_comm),
          m_matrix(comm),
          redist(comm),
          use_add(use_add)
        {}

        void finalize(const int n_local_rows, const int n_local_cols)
        {
            ownership_ranges_rows.resize(comm.size() + 1);
            ownership_ranges_cols.resize(comm.size() + 1);

            std::fill(ownership_ranges_rows.begin(), ownership_ranges_rows.end(), 0);
            std::fill(ownership_ranges_cols.begin(), ownership_ranges_cols.end(), 0);

            ownership_ranges_rows[comm.rank() + 1] = n_local_rows;
            ownership_ranges_cols[comm.rank() + 1] = n_local_cols;

            comm.all_reduce(&ownership_ranges_rows[0], ownership_ranges_rows.size(), moonolith::MPISum());
            comm.all_reduce(&ownership_ranges_cols[0], ownership_ranges_cols.size(),  moonolith::MPISum());

            std::partial_sum(ownership_ranges_rows.begin(), ownership_ranges_rows.end(), ownership_ranges_rows.begin());
            std::partial_sum(ownership_ranges_cols.begin(), ownership_ranges_cols.end(), ownership_ranges_cols.begin());

            if(use_add) {
                redist.apply(ownership_ranges_rows, m_matrix, moonolith::AddAssign<double>());
            } else {
                redist.apply(ownership_ranges_rows, m_matrix, moonolith::Assign<double>());
            }
        }

        void finalize(const int n_local_rows)
        {
            ownership_ranges_rows.resize(comm.size() + 1);
            std::fill(ownership_ranges_rows.begin(), ownership_ranges_rows.end(), 0);
            ownership_ranges_rows[comm.rank() + 1] = n_local_rows;
            comm.all_reduce(&ownership_ranges_rows[0], ownership_ranges_rows.size(), moonolith::MPISum());
            std::partial_sum(ownership_ranges_rows.begin(), ownership_ranges_rows.end(), ownership_ranges_rows.begin());

            if(use_add) {
                redist.apply(ownership_ranges_rows, m_matrix, moonolith::AddAssign<double>());
            } else {
                redist.apply(ownership_ranges_rows, m_matrix, moonolith::Assign<double>());
            }
        }

        template<typename IDX, class ElementMatrix>
        void add(
            const std::vector<IDX> &rows,
            const std::vector<IDX> &cols,
            ElementMatrix &mat
            )
        {   
            std::size_t n_rows = rows.size();
            std::size_t n_cols = cols.size();

            for(std::size_t i = 0; i < n_rows; ++i) {
                auto dof_I = rows[i];
                for(std::size_t j = 0; j < n_cols; ++j) {
                    auto dof_J = cols[j];
                    m_matrix.add(dof_I, dof_J, mat(i, j));
                }
            }
        }


        template<typename IDX>
        void insert(const IDX &idx, const double val)
        {
            if(use_add) {
                m_matrix.add(idx, 0, val);
            } else {
                if(val != 0.) {
                    m_matrix.set(idx, 0, val);
                }
            }
        }

        template<typename IDX>
        void insert(const std::vector<IDX> &idx, const double val)
        {
            if(use_add) {

                for(auto i : idx) {
                    m_matrix.add(i, 0, val);
                }
            } else {
                if(val != 0.) {
                    for(auto i : idx) {
                        m_matrix.set(i, 0, val);
                    }
                }
            }
        }

        template<typename IDX, class ElementMatrix>
        void insert(
            const std::vector<IDX> &rows,
            const std::vector<IDX> &cols,
            ElementMatrix &mat
            )
        {
            if(use_add) {
                add(rows, cols, mat);
            } else {
                set_non_zero(rows, cols, mat);
            }
        }

        template<typename IDX, class ElementMatrix>
        void set(
            const std::vector<IDX> &rows,
            const std::vector<IDX> &cols,
            ElementMatrix &mat
            )
        {   
            std::size_t n_rows = rows.size();
            std::size_t n_cols = cols.size();

            for(std::size_t i = 0; i < n_rows; ++i) {
                auto dof_I = rows[i];
                for(std::size_t j = 0; j < n_cols; ++j) {
                    auto dof_J = cols[j];
                    m_matrix.set(dof_I, dof_J, mat(i, j));
                }
            }
        }

        template<typename IDX, class ElementMatrix>
        void set_non_zero(
            const std::vector<IDX> &rows,
            const std::vector<IDX> &cols,
            ElementMatrix &mat
            )
        {   
            std::size_t n_rows = rows.size();
            std::size_t n_cols = cols.size();

            for(std::size_t i = 0; i < n_rows; ++i) {
                auto dof_I = rows[i];
                for(std::size_t j = 0; j < n_cols; ++j) {
                    auto dof_J = cols[j];

                    if(std::abs(mat(i, j)) != 0.) {
                        m_matrix.set(dof_I, dof_J, mat(i, j));
                    }
                }
            }
        }

        template<typename IDX, class ElementVector>
        void set(
            const std::vector<IDX> &rows,
            ElementVector &vec
            )
        {   
            std::size_t n_rows = rows.size();

            for(std::size_t i = 0; i < n_rows; ++i) {
                auto dof_I = rows[i];
                m_matrix.set(dof_I, 0, vec(i));
            }
        }


        template<typename IDX, class ElementVector>
        void set_non_zero(
            const std::vector<IDX> &rows,
            ElementVector &vec
            )
        {   
            std::size_t n_rows = rows.size();

            for(std::size_t i = 0; i < n_rows; ++i) {
                auto dof_I = rows[i];

                if(std::abs(vec(i)) != 0.) {
                    m_matrix.set(dof_I, 0, vec(i));
                }
            }
        }

        template<typename IDX>
        void set_non_zero(
            const std::vector<IDX> &rows,
            std::vector<double> &vec
            )
        {   
            std::size_t n_rows = rows.size();

            for(std::size_t i = 0; i < n_rows; ++i) {
                auto dof_I = rows[i];

                if(std::abs(vec[i]) != 0.) {
                    m_matrix.set(dof_I, 0, vec[i]);
                }
            }
        }

        template<typename IDX, class ElementVector>
        void add(
            const std::vector<IDX> &rows,
            ElementVector &vec
            )
        {   
            std::size_t n_rows = rows.size();

            for(std::size_t i = 0; i < n_rows; ++i) {
                auto dof_I = rows[i];
                m_matrix.add(dof_I, 0, vec(i));
            }
        }

        template<typename IDX>
        void add(
            const std::vector<IDX> &rows,
            std::vector<double> &vec
            )
        {   
            std::size_t n_rows = rows.size();

            for(std::size_t i = 0; i < n_rows; ++i) {
                auto dof_I = rows[i];
                m_matrix.add(dof_I, 0, vec[i]);
            }
        }

        template<typename IDX, class ElementVector>
        void insert(
            const std::vector<IDX> &rows,
            ElementVector &vec
            )
        {   
            if(use_add) {
                add(rows, vec);
            } else {
                set_non_zero(rows, vec);
            }
        }

        template<typename IDX, class ElementVector>
        void add_tensor_product_idx(
            const std::vector<IDX> &rows,
            const int &tensor_dim,
            ElementVector &vec
            )
        {   
            std::size_t n_rows = rows.size();

            IDX idx = 0;
            for(std::size_t i = 0; i < n_rows; ++i) {
                auto dof_I = rows[i];
                auto dof_I_x_d = dof_I * tensor_dim;
               
                for(IDX k = 0; k < tensor_dim; ++k, ++idx) {
                    m_matrix.add(dof_I_x_d + k, 0, vec(idx));
                }
            }
        }

        template<typename IDX, class ElementVector>
        void set_non_zero_tensor_product_idx(
            const std::vector<IDX> &rows,
            const int &tensor_dim,
            ElementVector &vec
            )
        {   
            std::size_t n_rows = rows.size();

            IDX idx = 0;
            for(std::size_t i = 0; i < n_rows; ++i) {
                auto dof_I = rows[i];
                auto dof_I_x_d = dof_I * tensor_dim;
               
                for(IDX k = 0; k < tensor_dim; ++k, ++idx) {
                    if(std::abs(vec(idx)) != 0.0) {
                        m_matrix.set(dof_I_x_d + k, 0, vec(idx));
                    }
                }
            }
        }

        template<typename IDX, class ElementVector>
        void insert_tensor_product_idx(
            const std::vector<IDX> &rows,
            const int &tensor_dim,
            ElementVector &vec
        )
        {
            if(use_add) {
                add_tensor_product_idx(rows, tensor_dim, vec);
            } else {
                set_non_zero_tensor_product_idx(rows, tensor_dim, vec);
            }
        } 

        void fill(USparseMatrix &mat)
        {
            auto nnz = m_matrix.local_max_entries_x_col();
            auto n_local_rows = ownership_ranges_rows[comm.rank() + 1] - ownership_ranges_rows[comm.rank()];
            auto n_local_cols = ownership_ranges_cols[comm.rank() + 1] - ownership_ranges_cols[comm.rank()];
            mat = local_sparse(n_local_rows, n_local_cols, nnz);

            {
                utopia::Write<utopia::USparseMatrix> write(mat);
                for (auto it = m_matrix.iter(); it; ++it) {
                    mat.set(it.row(), it.col(), *it);
                }
            }
        }

        void fill(UVector &vec)
        {
            auto n_local_rows = ownership_ranges_rows[comm.rank() + 1] - ownership_ranges_rows[comm.rank()];
            vec = local_zeros(n_local_rows);
            {
                Write<UVector> w_g(vec);

                for(auto it = m_matrix.iter(); it; ++it) {
                    assert(it.col() == 0);
                    vec.set(it.row(), *it);
                }
            }
        }

        //remove row variants (incomplete intersections)
        void fill(const std::vector<bool> &remove_row, USparseMatrix &mat)
        {
            auto nnz = m_matrix.local_max_entries_x_col();
            auto n_local_rows = ownership_ranges_rows[comm.rank() + 1] - ownership_ranges_rows[comm.rank()];
            auto n_local_cols = ownership_ranges_cols[comm.rank() + 1] - ownership_ranges_cols[comm.rank()];
            mat = local_sparse(n_local_rows, n_local_cols, nnz);

            utopia::Write<utopia::USparseMatrix> write(mat);

            for (auto it = m_matrix.iter(); it; ++it) {

                const SizeType index = it.row() - ownership_ranges_rows[comm.rank()];
                assert(index < remove_row.size());

                if(!remove_row[index]) {
                    mat.set(it.row(), it.col(), *it);
                }
            }
        }

        void fill(const std::vector<bool> &remove_row, UVector &vec)
        {
            auto n_local_rows = ownership_ranges_rows[comm.rank() + 1] - ownership_ranges_rows[comm.rank()];

            vec = local_zeros(n_local_rows);
            
            {
                Write<UVector> w_g(vec);

                for (auto it = m_matrix.iter(); it; ++it) {
                    const SizeType index = it.row() - ownership_ranges_rows[comm.rank()];
                    assert(index < remove_row.size());
                    assert(it.col() == 0);

                    if(!remove_row[index]) {
                        vec.set(it.row(), *it);
                    }
                }
            }
        }

        moonolith::Communicator comm;
        moonolith::SparseMatrix<double> m_matrix;
        moonolith::Redistribute< moonolith::SparseMatrix<double> > redist;
        std::vector<moonolith::Integer> ownership_ranges_rows, ownership_ranges_cols;
        bool use_add;
    };

}

#endif //UTOPIA_MATRIX_INSERTER_HPP
