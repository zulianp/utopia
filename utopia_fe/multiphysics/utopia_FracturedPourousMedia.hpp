#ifndef UTOPIA_FRACTURED_POUROUS_MEDIA_HPP
#define UTOPIA_FRACTURED_POUROUS_MEDIA_HPP

#include "utopia_PourousMatrix.hpp"
#include "utopia_DiscreteFractureNetwork.hpp"
#include "utopia_NewTransferAssembler.hpp"
#include "utopia_TransferUtils.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"

#include "utopia_fe_base.hpp"

namespace utopia {

    template<class Matrix, class Vector>
    class LagrangeMultiplier : public Configurable {
    public:
        using FunctionSpace = LibMeshFunctionSpace;

        void read(Input &in) override
        {
            in.get("type", type_);
            in.get("use-interpolation", use_interpolation_);
            
            if(is_dual()) {
                space_ = nullptr;
                return;
            } else if(is_same()) {
                space_ = nullptr;
            } else if(is_different()) {
                in.get("space", ui_space_);
                space_ = ui_space_.space().subspace_ptr(0);
            }
        }

        inline bool is_dual() const
        {
            return type_ == "dual";
        }

        ///@brief same as fracture network
        inline bool is_same() const
        {
            return type_ == "same";
        }

        inline bool is_different() const
        {
            return type_ == "different";
        }

        LagrangeMultiplier(libMesh::Parallel::Communicator &comm)
        : type_("dual"), use_interpolation_(false), ui_mesh_(comm), ui_space_(make_ref(ui_mesh_))
        {}

        bool init(
            FunctionSpace &pourous_matrix,
            FunctionSpace &fracture_newtork)
        {
            coupling_matrix_ = std::make_shared<Matrix>();
            mass_matrix_     = std::make_shared<Matrix>();

            if(is_dual()) {
                NewTransferAssembler transfer_assembler;
                transfer_assembler.remove_incomplete_intersections(false);

                if(!transfer_assembler.assemble(
                    pourous_matrix.mesh(),
                    pourous_matrix.dof_map(), 
                    fracture_newtork.mesh(),
                    fracture_newtork.dof_map())) {
                    return false;
                }

                auto op = transfer_assembler.build_operator();
                coupling_matrix_ = op->matrix();

                Vector d = sum(*coupling_matrix_, 1);
                e_pseudo_inv(d, d, 1e-16);

                *mass_matrix_ = diag(d);

            } else if(is_same()) {
               
                if(use_interpolation_) {
                    
                    assemble_interpolation(
                        pourous_matrix,
                        fracture_newtork,
                        *coupling_matrix_,
                        *mass_matrix_
                    );

                } else {
                   
                    assemble_projection(
                        pourous_matrix,
                        fracture_newtork,
                        *coupling_matrix_,
                        *mass_matrix_
                    );
                }

            } else {

                assemble_projection(
                    pourous_matrix,
                    fracture_newtork,
                    *space_,
                    *coupling_matrix_,
                    *mass_matrix_
                );
            }

            (*mass_matrix_) *= -1.0;

            coupling_matrix_t_ = std::make_shared<Matrix>();
            mass_matrix_t_     = std::make_shared<Matrix>();


            (*coupling_matrix_t_) = transpose(*coupling_matrix_);
            (*mass_matrix_t_)     = transpose(*mass_matrix_);

            set_zero_at_constraint_rows(pourous_matrix.dof_map(),   *coupling_matrix_t_);
            set_zero_at_constraint_rows(fracture_newtork.dof_map(), *mass_matrix_t_);
            return true;
        }

        inline std::shared_ptr<Matrix> transfer_matrix() const
        {
            return coupling_matrix_;
        }

        inline std::shared_ptr<Matrix> coupling_matrix() const
        {
            return coupling_matrix_;
        }

        inline std::shared_ptr<Matrix> mass_matrix() const
        {
            return mass_matrix_;
        }

        inline std::shared_ptr<Matrix> coupling_matrix_t() const
        {
            return coupling_matrix_t_;
        }

        inline std::shared_ptr<Matrix> mass_matrix_t() const
        {
            return mass_matrix_t_;
        }

    private:
        std::string type_;
        bool use_interpolation_;
        std::shared_ptr<FunctionSpace> space_;

        UIMesh<libMesh::DistributedMesh> ui_mesh_;
        UIFunctionSpace<FunctionSpace> ui_space_;


        std::shared_ptr<Matrix> coupling_matrix_;
        std::shared_ptr<Matrix> mass_matrix_;

        std::shared_ptr<Matrix> coupling_matrix_t_;
        std::shared_ptr<Matrix> mass_matrix_t_;

    };

    template<class Matrix, class Vector>
    class FracturedPourousMedia : /*public Model<Matrix, Vector>*/ public Configurable {
    public:
        using PourousMatrix           = utopia::PourousMatrix<Matrix, Vector>;
        using DiscreteFractureNetwork = utopia::DiscreteFractureNetwork<Matrix, Vector>;
        using LagrangeMultiplier      = utopia::LagrangeMultiplier<Matrix, Vector>;

        void read(Input &in) override
        {
            //read in the matrix model
            in.get("pourous-matrix",   pourous_matrix_);
    
            //read in the fracture-network models (1 or more)
            in.get("fracture-networks", [this](Input &in) {
                in.get_all([this](Input &in) {
                    auto dfn = std::make_shared<DiscreteFractureNetwork>(this->comm_);
                    dfn->read(in);
                    fracture_network_.push_back(dfn);

                    auto lm = std::make_shared<LagrangeMultiplier>(this->comm_); 

                    in.get("multiplier", *lm);
                    lagrange_multiplier_.push_back(lm);

                    lm->init(pourous_matrix_.space(), dfn->space());
                });;
            });

            in.get("assembly-strategy", assembly_strategy_);

            std::cout << "n fracture networks " << fracture_network_.size() << std::endl;

        }

        inline bool compute_flow()
        {
            Matrix A;
            Vector rhs;

            if(!assemble_flow(A, rhs)) {
                std::cerr << "[Error] failed to assemble" << std::endl;
                return false;
            }

            // if(use_mg_) {
            //      //TODO
            // }

            Vector x = local_zeros(local_size(rhs));
            Factorization<Matrix, Vector> solver;
            
            if(!solver.solve(A, rhs, x)) {
                std::cerr << "[Error] failed to solve" << std::endl;
                return false;
            }

            disassemble_flow(x);


            export_flow();
            return true;
        }

        inline bool assemble_flow(Matrix &A, Vector &rhs)
        {
            A_m   = std::make_shared<Matrix>();
            rhs_m = std::make_shared<Vector>(); 
            x_m   = std::make_shared<Vector>(); 

            auto &dof_map_m = pourous_matrix_.space().dof_map();

            *x_m = local_zeros(dof_map_m.n_local_dofs());
            pourous_matrix_.assemble_flow(*x_m, *A_m, *rhs_m);

            apply_boundary_conditions(dof_map_m, *A_m, *rhs_m);

            if(assembly_strategy_ == "static-condensation") {
                A   = *A_m;
                rhs = *rhs_m;

                const std::size_t n_dfn = fracture_network_.size();

                assert(n_dfn == lagrange_multiplier_.size());

                A_f.resize(n_dfn);
                rhs_f.resize(n_dfn);
                x_f.resize(n_dfn);

                for(std::size_t i = 0; i < n_dfn; ++i) {
                    
                    auto dfn = fracture_network_[i];
                    auto &dof_map_f = dfn->space().dof_map();

                    A_f[i]   = std::make_shared<Matrix>();
                    x_f[i]   = std::make_shared<Vector>();
                    rhs_f[i] = std::make_shared<Vector>();

                    *x_f[i] = local_zeros(dof_map_f.n_local_dofs());
                    dfn->assemble_flow(*x_f[i], *A_f[i], *rhs_f[i]);

                    apply_boundary_conditions(dof_map_f, *A_f[i], *rhs_f[i]);

                    auto T = lagrange_multiplier_[i]->transfer_matrix();

                    Matrix A_temp = transpose(*T) * (*A_f[i]) * (*T);
                    A += A_temp;
                    rhs += transpose(*T) * (*rhs_f[i]);
                }

            } else //if(assembly_strategy_ == "monolithic") 
            {
                assert(false && "TODO");
                return false;
            }


            return true;
        }

        inline void disassemble_flow(const Vector &x)
        {
            if(assembly_strategy_ == "static-condensation") {
                *x_m = x;

                const std::size_t n_dfn = fracture_network_.size();

                for(std::size_t i = 0; i < n_dfn; ++i) {
                    auto T = lagrange_multiplier_[i]->transfer_matrix();
                    *x_f[i] = (*T) * x;
                }
            }  else //if(assembly_strategy_ == "monolithic") 
            {
                assert(false && "TODO");
            }
        }

        inline bool export_flow()
        {
            write("porous_matrix.e", pourous_matrix_.space(), *x_m);

            const std::size_t n_dfn = fracture_network_.size();

            for(std::size_t i = 0; i < n_dfn; ++i) {
                write("fracture_network_" + std::to_string(i) + ".e", fracture_network_[i]->space(), *x_f[i]);
            }

            return true;
        }

        FracturedPourousMedia(libMesh::Parallel::Communicator &comm)
        : comm_(comm), pourous_matrix_(comm), assembly_strategy_("monolithic"), use_mg_(false)
        {}

    private:
        libMesh::Parallel::Communicator &comm_;
        PourousMatrix pourous_matrix_;
        std::vector<std::shared_ptr<DiscreteFractureNetwork>> fracture_network_;
        std::vector<std::shared_ptr<LagrangeMultiplier>> lagrange_multiplier_;

        std::shared_ptr<Matrix> A_m;
        std::shared_ptr<Vector> rhs_m;
        std::shared_ptr<Vector> x_m;

        std::vector<std::shared_ptr<Matrix>> A_f;
        std::vector<std::shared_ptr<Vector>> rhs_f;
        std::vector<std::shared_ptr<Vector>> x_f;

        std::string assembly_strategy_;
        bool use_mg_;
    };
    
}

#endif //UTOPIA_FRACTURED_POUROUS_MEDIA_HPP
