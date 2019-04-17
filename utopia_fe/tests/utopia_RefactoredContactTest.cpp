#include "utopia_RefactoredContactTest.hpp"


#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_Newmark.hpp"
#include "utopia_LibMeshBackend.hpp"
#include "utopia_ContactStabilizedNewmark.hpp"
#include "utopia_ui.hpp"
#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"
#include "utopia_UIContactParams.hpp"
#include "utopia_UIMaterial.hpp"
#include "utopia_UIScalarSampler.hpp"
#include "utopia_InputParameters.hpp"
#include "utopia_polymorphic_QPSolver.hpp"

#include "libmesh/mesh_refinement.h"

#include "utopia_ContactAssembler.hpp"
#include "utopia_LibMeshShape.hpp"
#include "moonolith_affine_transform.hpp"
#include "moonolith_contact.hpp"
#include "moonolith_sparse_matrix.hpp"
#include "moonolith_redistribute.hpp"
#include "moonolith_assign_functions.hpp"

#include <vector>
#include <memory>

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

    class ContactData {
    public:
        MatrixInserter B, D;
        MatrixInserter gap, normal;

        USparseMatrix B_x, D_x;
        UVector weighted_gap_x, gap_x, normal_xyz;
        UVector diag_inv_D_x;

        inline MPI_Comm comm() const
        {
            return B.comm.get();
        }

        void finalize(
            const SizeType n_local_dofs,
            const SizeType spatial_dim)
        {
            B.finalize(n_local_dofs, n_local_dofs);
            D.finalize(n_local_dofs, n_local_dofs);
            gap.finalize(n_local_dofs);
            normal.finalize(n_local_dofs * spatial_dim);

            B.fill(B_x);
            D.fill(D_x);
            gap.fill(weighted_gap_x);
            normal.fill(normal_xyz);

            double sum_B_x = sum(B_x);
            double sum_D_x = sum(D_x);

            UVector d = sum(D_x, 1);
            diag_inv_D_x = local_zeros(local_size(d));

            {
                Write<UVector> w(diag_inv_D_x);

                each_read(d, [this](const SizeType i, const double value) {
                    if(std::abs(value) > 1e-15) {
                        diag_inv_D_x.set(i, 1./value);
                    }

                });
            }

            gap_x = e_mul(diag_inv_D_x, weighted_gap_x);

            std::cout << sum_B_x << " == " << sum_D_x << std::endl;
        }

        ContactData(MPI_Comm comm) : B(comm), D(comm), gap(comm), normal(comm) {}
    private:
        ContactData(const ContactData &other) : B(other.comm()), D(other.comm()), gap(other.comm()), normal(other.comm()) {}
    };


    template<int Dim>
    class SurfaceQuadratureConverter {
    public:
        using SubVector = moonolith::Vector<double, Dim-1>;
        using Vector    = moonolith::Vector<double, Dim>;

        //from libmesh to moonolith
        SubVector point_shift;
        double    point_rescale;
        double    weight_rescale;

        //from moonolith to libmesh
        double    trial_weight_rescale;
        double    test_weight_rescale;

        int current_order;

        void convert_master(const moonolith::Quadrature<double, Dim-1> &in, QMortar &out)
        {
            convert(in, trial_weight_rescale, out);
        }

        void convert_slave(const moonolith::Quadrature<double, Dim-1> &in, QMortar &out)
        {
            convert(in, test_weight_rescale, out);
        }

        SurfaceQuadratureConverter()
        : current_order(-1)
        {}

        void init(
            const libMesh::Elem &trial,
            const int trial_order,
            const libMesh::Elem &test,
            const int test_order,
            moonolith::Quadrature<double, Dim-1> &q)
        {
            const int order = order_for_l2_integral(Dim-1, trial, trial_order, test, test_order);

            if(order != current_order) {
                moonolith::fill(point_shift, 0.);

                if(Dim == 2) {  
                    libMesh::QGauss ir(1, libMesh::Order(order));

                    if(order <= 2) {
                        ir.init(libMesh::EDGE2);
                    } else {
                        ir.init(libMesh::EDGE4);
                    }

                    point_shift.x = 1;
                    point_rescale = 0.5;
                    weight_rescale = 0.5;

                    convert(ir, point_shift, point_rescale, weight_rescale, q);
                } else if(Dim == 3) {

                    libMesh::QGauss ir(2, libMesh::Order(order));

                    if(order <= 2) {
                        ir.init(libMesh::TRI3);
                    } else {
                        ir.init(libMesh::TRI6);
                    }

                    weight_rescale = 2.0;

                    convert(ir, point_shift, point_rescale, weight_rescale, q);

                } else {
                    assert(false);
                }

                current_order = order;
            }

            trial_weight_rescale = ref_volume(trial.type());
            test_weight_rescale  = ref_volume(test.type());
        }

    };

    template<int Dim>
    class ProjectionAlgorithm {
    public:
        using Trafo     = moonolith::AffineTransform<double, Dim-1, Dim>;
        using Shape     = moonolith::Shape<double, Dim-1, Dim>;
        using SubVector = moonolith::Vector<double, Dim-1>;
        using Vector    = moonolith::Vector<double, Dim>;

        ContactData &data;

        //algorithms
        moonolith::AffineContact<double, Dim> affine_contact;
        moonolith::WarpedContact<double, Dim> warped_contact;

        //buffers
        std::shared_ptr<Trafo> trafo_m, trafo_s;
        std::shared_ptr<Shape> shape_m, shape_s;
        SurfaceQuadratureConverter<Dim> converter;
        std::vector<PetscInt> dofs_petsc_s, dofs_petsc_m;

        //libmesh-buffers
        QMortar lm_q_master, lm_q_slave;
        std::unique_ptr<libMesh::FEBase> master_fe, slave_fe;
        libMesh::DenseMatrix<libMesh::Real> b_elmat, d_elmat;
        libMesh::DenseVector<libMesh::Real> gap_vec, normal_vec;

        double area = 0.;

        ProjectionAlgorithm(ContactData &data)
        : data(data), lm_q_master(Dim-1), lm_q_slave(Dim-1)
        {
            trafo_m = std::make_shared<Trafo>();
            trafo_s = std::make_shared<Trafo>();
            affine_contact.trafo_master = trafo_m;
            affine_contact.trafo_slave  = trafo_s;

            if(Dim == 2) {
                //somehow the order of surface-line elements is clockwise
                warped_contact.invert_plane_dir = true;

                //TODO check for the affine contact too
            }


        }

        void init_fe(int master_order, int slave_order)
        {
            if(!master_fe) {
                master_fe = libMesh::FEBase::build(Dim-1, libMesh::Order(master_order));
                slave_fe  = libMesh::FEBase::build(Dim-1, libMesh::Order(slave_order));

                master_fe->get_phi();
                slave_fe->get_phi();
                slave_fe->get_JxW();

                //initialize biorth if needed
            }
        }

        template<class Adapter>
        void assemble(
            const Adapter &master,
            const Adapter &slave,
            const moonolith::Quadrature<double, Dim-1> &q_master,
            const moonolith::Quadrature<double, Dim-1> &q_slave,
            const moonolith::Vector<double, Dim> &normal,
            const moonolith::Storage<double> &gap
        )
        {

            ///////////////////////////////////////////////////////
            //set-up

            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            auto &e_m = master.elem();
            auto &e_s = slave.elem();

            auto &dofs_m = master.dofs().global;
            auto &dofs_s = slave.dofs().global;

            converter.convert_master(q_master, lm_q_master);
            converter.convert_slave(q_slave, lm_q_slave);

            init_fe(m_m.fe_type(0).order, m_s.fe_type(0).order);

            master_fe->attach_quadrature_rule(&lm_q_master);
            master_fe->reinit(&e_m);

            slave_fe->attach_quadrature_rule(&lm_q_slave);
            slave_fe->reinit(&e_s);

            ///////////////////////////////////////////////////////
            //assemble local element matrices

            b_elmat.zero();
            mortar_assemble(*master_fe, *slave_fe, b_elmat);

            d_elmat.zero();
            mortar_assemble(*slave_fe, *slave_fe, d_elmat);

            gap_vec.zero();
            integrate_scalar_function(
                *slave_fe,
                gap,
                gap_vec
            );

            normal_vec.zero();
            l2_project_normal(*slave_fe, normal, normal_vec);


            ///////////////////////////////////////////////////////
            //local to global

            data.B.insert(dofs_s, dofs_m, b_elmat);
            data.D.insert(dofs_s, dofs_s, d_elmat);
            data.gap.insert(dofs_s, gap_vec);
            data.normal.insert_tensor_product_idx(dofs_s, Dim, normal_vec);


            ///////////////////////////////////////////////////////
            //check on the area

            auto partial_sum = std::accumulate(
                b_elmat.get_values().begin(),
                b_elmat.get_values().end(),
                libMesh::Real(0.0)
            );

            assert(partial_sum > 0.);
            area += partial_sum;
        }

        template<class Adapter>
        bool apply(Adapter &master, Adapter &slave)
        {
            auto &m_m = master.collection();
            auto &m_s = slave.collection();

            auto &e_m = master.elem();
            auto &e_s = slave.elem();

            // const bool is_affine = e_m.has_affine_map() && e_s.has_affine_map();

            //force usage of non-affine code
            const bool is_affine = false;

            if(is_affine) {
                //AFFINE CONTACT
                make(e_m, affine_contact.master);
                make(e_s, affine_contact.slave);

                make_transform(e_m, *trafo_m);
                make_transform(e_s, *trafo_s);

                converter.init(
                   e_m,
                   m_m.fe_type(0).order,
                   e_s,
                   m_s.fe_type(0).order,
                   affine_contact.q_rule
                );

                if(affine_contact.compute()) {
                    auto slave_area = moonolith::measure(affine_contact.slave);
                    
                    assemble(
                        master,
                        slave,
                        affine_contact.q_master,
                        affine_contact.q_slave,
                        affine_contact.normal(),
                        affine_contact.gap
                    );

                    return true;
                } else {
                    return false;
                }

            } else {
                //WARPED CONTACT
                bool use_newton = false;
                // warped_contact.shape_master= std::make_shared<LibMeshShape<double, Dim>>(e_m, m_m.libmesh_fe_type(0), use_newton);               
                auto lm_shape_master= std::make_shared<LibMeshShape<double, Dim>>(e_m, m_m.libmesh_fe_type(0), use_newton);               
                // lm_shape_master->verbose(true);
                warped_contact.shape_master = lm_shape_master;
                warped_contact.shape_slave = std::make_shared<LibMeshShape<double, Dim>>(e_s, m_s.libmesh_fe_type(0), use_newton);

                make_non_affine(e_m, warped_contact.master);
                make_non_affine(e_s, warped_contact.slave);

                converter.init(
                   e_m,
                   m_m.fe_type(0).order,
                   e_s,
                   m_s.fe_type(0).order,
                   warped_contact.q_rule
                );

                if(warped_contact.compute()) {

                   assemble(
                       master,
                       slave,
                       warped_contact.q_master,
                       warped_contact.q_slave,
                       warped_contact.normal(),
                       warped_contact.gap
                   );

                    return true;
                } else {
                    return false;
                }
            }

            return false;
        }
    };

    template<int Dim>
    bool run_contact(
        const ContactParams &params,
        LibMeshFunctionSpaceAdapter &adapter,
        ContactData &contact_data)
    {
        using AlogrithmT = moonolith::SingleCollectionOneMasterOneSlaveAlgorithm<Dim, LibMeshFunctionSpaceAdapter>;
        using Adapter    = typename AlogrithmT::Adapter;

        auto cm = std::make_shared<LibMeshCollectionManagerT>(adapter.comm()); 

        moonolith::Communicator m_comm(cm->comm.get());

        moonolith::SearchSettings s;
        // s.verbosity_level = 3;
        // s.disable_redistribution = true;
        AlogrithmT algo(m_comm, cm, s);

        algo.init(
            adapter,
            params.contact_pair_tags,
            params.search_radius
        );

        ProjectionAlgorithm<Dim> contact_algo(contact_data);

        // Write<UVector> w(contact_data.is_contact);

        algo.compute([&](const Adapter &master, const Adapter &slave) -> bool {
            return contact_algo.apply(master, slave);
        });

        m_comm.all_reduce(&contact_algo.area, 1, moonolith::MPISum());
        std::cout << "area: " << contact_algo.area << std::endl;


        auto n_local_dofs = adapter.n_local_dofs();
        contact_data.finalize(n_local_dofs, Dim);

        return contact_algo.area > 0.;
    }

    void RefactoredContactTest::run(Input &in) {

        using ProductSpaceT    = utopia::ProductFunctionSpace<LibMeshFunctionSpace>;
        using MaterialT        = utopia::UIMaterial<ProductSpaceT, USparseMatrix, UVector>;
        using ForcingFunctionT = UIForcingFunction<ProductSpaceT, UVector>;

        in.get("contact-problem", [&](Input &in) { 
            UIMesh<libMesh::DistributedMesh> mesh(this->comm());
            UIFunctionSpace<LibMeshFunctionSpace> space(make_ref(mesh));
            UIContactParams params;
            std::shared_ptr< ElasticMaterial<USparseMatrix, UVector> > model;

            in.get("mesh", mesh);
            in.get("space", space);
            in.get("contact", params);

            model = make_unique<MaterialT>(space.space());

            auto &V = space.space().subspace(0);


            LibMeshFunctionSpaceAdapter adapter;

            bool is_volume = V.mesh().spatial_dimension() == V.mesh().mesh_dimension();

            if(is_volume) {
           
                adapter.extract_surface_init(
                    make_ref(V.mesh()),
                    V.dof_map(),
                    params.contact_params.variable_number
                );
            } else {
                //shell mesh
                adapter.init(
                    make_ref(V.mesh()),
                    V.dof_map(),
                    params.contact_params.variable_number
                );
            }

            adapter.print_tags();

            ContactData contact_data(V.mesh().comm().get());

            bool found_contact = false;
            if(V.mesh().spatial_dimension() == 2) {
                found_contact = run_contact<2>(params.contact_params, adapter, contact_data);
            } 
            else if(V.mesh().spatial_dimension() == 3) {
                found_contact = run_contact<3>(params.contact_params, adapter, contact_data);
            }

            assert(found_contact);



            if(is_volume) {
                UVector x = (*adapter.permutation()) * contact_data.gap_x;
                write("warped.e", V, x);
            }

           
        });
        
    }
}

