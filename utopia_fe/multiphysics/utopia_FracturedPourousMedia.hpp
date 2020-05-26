#ifndef UTOPIA_FRACTURED_POUROUS_MEDIA_HPP
#define UTOPIA_FRACTURED_POUROUS_MEDIA_HPP

#include "utopia_DiscreteFractureNetwork.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_FractureFlowUtils.hpp"
#include "utopia_NewTransferAssembler.hpp"
#include "utopia_PourousMatrix.hpp"
#include "utopia_TransferUtils.hpp"

#include "utopia_fe_base.hpp"

namespace utopia {

    void remove_constrained_dofs(libMesh::DofMap &dof_map, USparseMatrix &mat, UVector &vec) {
        if (utopia::Utopia::instance().verbose()) {
            std::cout << "apply_boundary_conditions begin: " << std::endl;
        }

        Chrono c;
        c.start();

        assert(!empty(mat));
        assert(!empty(vec));

        using SizeType = Traits<UVector>::SizeType;

        const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

        Size ls = local_size(mat);
        Size s = size(mat);

        std::vector<SizeType> index;

        Range rr = range(vec);

        if (has_constaints) {
            for (SizeType i = rr.begin(); i < rr.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    index.push_back(i);
                }
            }
        }

        set_zero_rows(mat, index, 0.);

        Write<UVector> w_v(vec);

        if (has_constaints) {
            libMesh::DofConstraintValueMap &rhs_values = dof_map.get_primal_constraint_values();

            Range r = range(vec);
            for (SizeType i = r.begin(); i < r.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    // auto valpos = rhs_values.find(i);
                    vec.set(i, 0.0);
                }
            }
        }

        c.stop();

        if (utopia::Utopia::instance().verbose()) {
            std::cout << "apply_boundary_conditions end: " << c << std::endl;
        }
    }

    // void remove_constrained_dofs(libMesh::DofMap &dof_map, USparseMatrix &mat)
    // {
    //     if(utopia::Utopia::instance().verbose()) {
    //         std::cout << "apply_boundary_conditions begin: "  << std::endl;
    //     }

    //     Chrono c;
    //     c.start();

    //     assert(!empty(mat));

    //     using SizeType = Traits<UVector>::SizeType;

    //     const bool has_constaints = dof_map.constraint_rows_begin() !=
    //     dof_map.constraint_rows_end();

    //     Size ls = local_size(mat);
    //     Size s = size(mat);

    //     std::vector<SizeType> index;

    //     Range rr = row_range(mat);

    //     if(has_constaints) {
    //         for(SizeType i = rr.begin(); i < rr.end(); ++i) {
    //             if( dof_map.is_constrained_dof(i) ) {
    //                 index.push_back(i);
    //             }
    //         }
    //     }

    //     set_zero_rows(mat, index, 0.);

    //     c.stop();

    //     if(utopia::Utopia::instance().verbose()) {
    //         std::cout << "apply_boundary_conditions end: " << c << std::endl;
    //     }
    // }

    template <class Matrix, class Vector>
    class LagrangeMultiplier : public Configurable {
    public:
        using FunctionSpace = LibMeshFunctionSpace;
        using Scalar = UTOPIA_SCALAR(Vector);

        void read(Input &in) override {
            in.get("type", type_);
            in.get("use-interpolation", use_interpolation_);
            opts_.read(in);

            if (is_dual()) {
                space_ = nullptr;
                return;
            } else if (is_same()) {
                space_ = nullptr;
            } else if (is_different()) {
                in.get("space", ui_space_);
                space_ = ui_space_.space().subspace_ptr(0);
            }
        }

        inline bool is_dual() const { return type_ == "dual"; }

        ///@brief same as fracture network
        inline bool is_same() const { return type_ == "same"; }

        inline bool is_different() const { return type_ == "different"; }

        LagrangeMultiplier(libMesh::Parallel::Communicator &comm)
            : type_("dual"),
              use_interpolation_(false),
              handle_adaptivity_(false),
              ui_mesh_(comm),
              ui_space_(make_ref(ui_mesh_)) {}

        bool init(FunctionSpace &pourous_matrix,
                  const std::shared_ptr<USparseMatrix> &constraint_matrix_pm,
                  FunctionSpace &fracture_newtork) {
            coupling_matrix_ = std::make_shared<Matrix>();
            mass_matrix_ = std::make_shared<Matrix>();

            if (is_dual()) {
                NewTransferAssembler transfer_assembler;
                transfer_assembler.remove_incomplete_intersections(false);
                transfer_assembler.handle_adaptive_refinement(handle_adaptivity_);

                if (constraint_matrix_pm) {
                    transfer_assembler.constraint_matrix_from(constraint_matrix_pm);
                }

                if (constraint_matrix_pm && !empty(*constraint_matrix_pm)) {
                    m_utopia_status("using constraint matrix from mortar in porous-matrix");
                }

                if (!transfer_assembler.assemble(pourous_matrix.mesh(),
                                                 pourous_matrix.dof_map(),
                                                 fracture_newtork.mesh(),
                                                 fracture_newtork.dof_map(),
                                                 opts_)) {
                    return false;
                }

                auto op = transfer_assembler.build_operator();
                coupling_matrix_ = op->matrix();

                inverse_mass_vector_ = sum(*coupling_matrix_, 1);
                e_pseudo_inv(inverse_mass_vector_, inverse_mass_vector_, 1e-16);

                *mass_matrix_ = diag(inverse_mass_vector_);

            } else if (is_same()) {
                if (use_interpolation_) {
                    assemble_interpolation(pourous_matrix, fracture_newtork, *coupling_matrix_, *mass_matrix_);

                } else {
                    assemble_projection(pourous_matrix, fracture_newtork, *coupling_matrix_, *mass_matrix_);
                }

            } else {
                assemble_projection(pourous_matrix, fracture_newtork, *space_, *coupling_matrix_, *mass_matrix_);
            }

            (*mass_matrix_) *= -1.0;

            coupling_matrix_t_ = std::make_shared<Matrix>();
            mass_matrix_t_ = std::make_shared<Matrix>();

            (*coupling_matrix_t_) = transpose(*coupling_matrix_);
            (*mass_matrix_t_) = transpose(*mass_matrix_);

            set_zero_at_constraint_rows(pourous_matrix.dof_map(), *coupling_matrix_t_);
            // set_zero_at_constraint_rows(fracture_newtork.dof_map(), *mass_matrix_t_);
            m_utopia_warning(
                "set_zero_at_constraint_rows(fracture_newtork.dof_map(), "
                "*mass_matrix_t_)");

            is_constrained_ = sum(*transfer_matrix(), 1);
            return true;
        }

        inline std::shared_ptr<Matrix> transfer_matrix() const { return coupling_matrix_; }

        inline std::shared_ptr<Matrix> coupling_matrix() const { return coupling_matrix_; }

        inline std::shared_ptr<Matrix> mass_matrix() const { return mass_matrix_; }

        inline std::shared_ptr<Matrix> coupling_matrix_t() const { return coupling_matrix_t_; }

        inline std::shared_ptr<Matrix> mass_matrix_t() const { return mass_matrix_t_; }

        inline const UVector &is_constrained() const { return is_constrained_; }

        inline void remove_mass(UVector &in_out) { in_out = in_out / inverse_mass_vector_; }

        inline void handle_adaptivity(const bool val) { handle_adaptivity_ = val; }

    private:
        std::string type_;
        bool use_interpolation_;
        bool handle_adaptivity_;

        std::shared_ptr<FunctionSpace> space_;

        UIMesh<libMesh::DistributedMesh> ui_mesh_;
        UIFunctionSpace<FunctionSpace> ui_space_;

        std::shared_ptr<Matrix> coupling_matrix_;
        std::shared_ptr<Matrix> mass_matrix_;

        std::shared_ptr<Matrix> coupling_matrix_t_;
        std::shared_ptr<Matrix> mass_matrix_t_;

        TransferOptions opts_;

        UVector is_constrained_;
        UVector inverse_mass_vector_;
    };

    template <class Matrix>
    class MatrixPostProcessor : public Configurable {
    public:
        virtual ~MatrixPostProcessor() {}
        virtual void post_process(const Matrix &mat) = 0;
    };

    template <class Matrix>
    class MatrixExporter final : public MatrixPostProcessor<Matrix> {
    public:
        inline void post_process(const Matrix &mat) override { write(path_, mat); }

        inline void read(Input &in) override {
            // in.get("name", name_);
            in.get("path", path_);
        }

        MatrixExporter()
            : path_("./mat.m")  //, name_("A")
        {}

    private:
        std::string path_;
        // std::string name_;
    };

    template <class Matrix>
    class ConditionNumberPrinter final : public MatrixPostProcessor<Matrix> {
    public:
        inline void post_process(const Matrix &mat) override {
            std::cout << "condition-number: ";
#ifdef WITH_SLEPC
            auto c = cond(mat);
            std::cout << c << std::endl;
#else
            std::cout << "[not computed, missing library]" << std::endl;
#endif  // WITH_SLEPC
        }

        inline void read(Input &) override {}

        ConditionNumberPrinter() {}
    };

    template <class Matrix, class Vector>
    class DFMReport final : public MatrixPostProcessor<Matrix> {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);
        using SizeType = UTOPIA_SIZE_TYPE(Vector);

        inline void post_process(const Matrix &mat) override {
            if (skip_cond_) {
                std::cout << "condition-number: [skipped]" << std::endl;
            } else {
                std::cout << "condition-number: ";
#ifdef WITH_SLEPC
                condition_number_ = cond(mat);
                std::cout << condition_number_ << std::endl;
#else
                std::cout << "[not computed, missing library]" << std::endl;
#endif  // WITH_SLEPC
            }
            n_dofs_ = size(mat).get(0);
        }

        inline void read(Input &in) override {
            in.get("path", path_);
            in.get("print-header", print_header_);
            in.get("skip-cond", skip_cond_);
        }

        DFMReport() : path_("report.csv"), print_header_(true), condition_number_(-1), n_dofs_(-1), skip_cond_(false) {}

        bool save() const {
            std::ofstream os(path_.c_str());

            if (!os.good()) {
                os.close();
                return false;
            }

            if (print_header_) {
                os << "n_dofs,cond_num";

                for (const auto &s : stats_) {
                    os << ",";
                    s.print_header(os);
                }

                os << "\n";
            }

            os << std::setprecision(9);
            os << n_dofs_ << "," << condition_number_;

            for (const auto &s : stats_) {
                os << ",";
                s.print(os);
            }

            os << "\n";

            os.close();
            return true;
        }

        void add_stat(const LibMeshFunctionSpace &space, const Matrix &mat, const std::string &name = "") {
            Stats s;
            if (name.empty()) {
                s.name = space.equation_system().name();
            } else {
                s.name = name;
            }

            auto it = space.mesh().active_local_elements_begin();
            if (it != space.mesh().active_local_elements_end()) {
                s.element_type = libMesh::Utility::enum_to_string((*it)->type());
            }
            s.n_dofs = space.dof_map().n_dofs();
            s.n_elements = space.mesh().n_active_elem();
            s.nnz = nnz(mat);

            stats_.push_back(s);
        }

        ~DFMReport() {}

    private:
        std::string path_;
        bool print_header_;
        Scalar condition_number_;
        SizeType n_dofs_;
        SizeType nnz_;
        bool skip_cond_;

        class Stats {
        public:
            void print_header(std::ostream &os) const {
                os << "elem_type_" << name << ",";
                os << "n_dofs_" << name << ",";
                os << "n_elements_" << name << ",";
                os << "nnz_" << name;
            }

            void print(std::ostream &os) const {
                os << element_type << ",";
                os << n_dofs << ",";
                os << n_elements << ",";
                os << nnz;
            }

            std::string name;
            std::string element_type;

            SizeType n_dofs;
            SizeType n_elements;
            SizeType nnz;

            Stats() : name("undefined"), element_type("undefined"), n_dofs(-1), n_elements(-1), nnz(-1) {}
        };

        std::vector<Stats> stats_;
    };

    template <class Matrix, class Vector>
    class FracturedPourousMedia : /*public Model<Matrix, Vector>*/ public Configurable {
    public:
        using PourousMatrix = utopia::PourousMatrix<Matrix, Vector>;
        using DiscreteFractureNetwork = utopia::DiscreteFractureNetwork<Matrix, Vector>;
        using LagrangeMultiplier = utopia::LagrangeMultiplier<Matrix, Vector>;
        using Scalar = UTOPIA_SCALAR(Vector);

        void read(Input &in) override {
            in.get("rescale", rescale_);

            pourous_matrix_.rescale(rescale_);

            in.get("adaptivity", adaptivity_);

            // read in the matrix model
            in.get("pourous-matrix", pourous_matrix_);

            // read in the fracture-network models (1 or more)
            in.get("fracture-networks", [this](Input &in) {
                in.get_all([this](Input &in) {
                    auto dfn = std::make_shared<DiscreteFractureNetwork>(this->comm_);
                    dfn->rescale(rescale_);

                    dfn->read(in);
                    fracture_network_.push_back(dfn);

                    auto lm = std::make_shared<LagrangeMultiplier>(this->comm_);
                    lm->handle_adaptivity(adaptivity_);
                    in.get("multiplier", *lm);
                    lagrange_multiplier_.push_back(lm);

                    // lm->init(pourous_matrix_.space(), dfn->space());
                });
                ;
            });

            in.get("assembly-strategy", assembly_strategy_);
            in.get("remove-constrained-dofs", remove_constrained_dofs_);

            std::cout << "n fracture networks " << fracture_network_.size() << std::endl;

            in.get("matrix-processors", [this](Input &in) {
                in.get_all([this](Input &in) {
                    std::string type;
                    in.get("type", type);

                    if (type == "export") {
                        auto exporter = std::make_shared<MatrixExporter<Matrix>>();
                        exporter->read(in);
                        matrix_processors_.push_back(std::move(exporter));
                    } else if (type == "condition-number") {
                        matrix_processors_.push_back(std::make_shared<ConditionNumberPrinter<Matrix>>());
                    }

                    // other post-processors....
                });
            });

            in.get("report", [this](Input &in) {
                report_ = std::make_shared<DFMReport<Matrix, Vector>>();
                report_->read(in);

                // matrix_processors_.push_back(report_);
            });

            std::string solver_type;
            in.get("solver-type", solver_type);

            if (solver_type == "ksp") {
                auto ksp = std::make_shared<KSPSolver<Matrix, Vector>>();
                ksp->max_it(6000);
                ksp->ksp_type("bcgs");
                ksp->pc_type("bjacobi");
                // ksp->set_preconditioner(std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>());
                solver_ptr = ksp;
            } else if (solver_type == "mg") {
                int n_levels = 2;
                in.get("n_levels", n_levels);

                auto linear_solver = std::make_shared<Factorization<Matrix, Vector>>();
                auto smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();

                auto mg = std::make_shared<SemiGeometricMultigrid>(smoother, linear_solver);
                mg->algebraic().fix_semidefinite_operators(true);

                mg->init(pourous_matrix_.space(), n_levels);
                mg->max_it(1);

                auto cg = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();
                cg->verbose(true);
                cg->set_preconditioner(mg);
                solver_ptr = cg;
            } else if (solver_type == "cggs") {
                auto cg = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();
                cg->max_it(20000);
                auto pgs = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
                pgs->max_it(40);

                cg->set_preconditioner(pgs);
                solver_ptr = cg;
            }

            in.get("solver", *solver_ptr);
        }

        inline bool init() {
            const std::size_t n_dfn = fracture_network_.size();

            assert(n_dfn == lagrange_multiplier_.size());

            std::shared_ptr<USparseMatrix> constraint_matrix;

            if (pourous_matrix_.has_mortar_constraints()) {
                auto &mortar = pourous_matrix_.mortar();
                mortar.compute_mortar_matrix_without_slave_dofs();
                constraint_matrix = mortar.mortar_matrix_without_slave_dofs();
            }

            for (std::size_t i = 0; i < n_dfn; ++i) {
                if (!lagrange_multiplier_[i]->init(
                        pourous_matrix_.space(), constraint_matrix, fracture_network_[i]->space())) {
                    return false;
                }
            }

            return true;
        }

        inline bool compute_flow() {
            if (adaptivity_) {
                return compute_flow_with_refinement();
            } else {
                return compute_flow_no_refinement();
            }
        }

        inline bool compute_flow_with_refinement() {
            int i = 0;
            do {
                std::cout << "compute_flow_with_refinement: loop " << i++ << std::endl;
                bool ok = compute_flow_no_refinement();
                if (!ok) return ok;
                // TODO fracture networks
            } while (pourous_matrix_.refine(*x_m_));

            return true;
        }

        inline bool compute_flow_no_refinement() {
            Matrix A;
            Vector rhs;

            if (!init()) {
                std::cerr << "[Error] failed to init Lagrange multiplier" << std::endl;
                return false;
            }

            if (!assemble_flow(A, rhs)) {
                std::cerr << "[Error] failed to assemble" << std::endl;
                return false;
            }

            // if(use_mg_) {
            //      //TODO
            // }

            Vector x = local_zeros(local_size(rhs));

            // solver_ptr->describe(std::cout);
            if (!solver_ptr->solve(A, rhs, x)) {
                std::cerr << "[Error] failed to solve" << std::endl;
                // return false;
            }

            disassemble_flow(x);

            pourous_matrix_.post_process_flow(*x_m_);

            export_flow();
            return true;
        }

        inline bool assemble_flow(Matrix &A, Vector &rhs) {
            A_m_ = std::make_shared<Matrix>();
            rhs_m_ = std::make_shared<Vector>();
            x_m_ = std::make_shared<Vector>();

            auto &dof_map_m = pourous_matrix_.space().dof_map();

            *x_m_ = local_zeros(dof_map_m.n_local_dofs());
            pourous_matrix_.assemble_flow(*x_m_, *A_m_, *rhs_m_);

            // apply_boundary_conditions(dof_map_m, *A_m_, *rhs_m_);
            apply_boundary_conditions(pourous_matrix_.space(), *A_m_, *rhs_m_);

            if (report_) {
                report_->add_stat(pourous_matrix_.space(), *A_m_);
            }

            if (assembly_strategy_ == "static-condensation") {
                A = *A_m_;
                rhs = *rhs_m_;

                const std::size_t n_dfn = fracture_network_.size();

                assert(n_dfn == lagrange_multiplier_.size());

                A_f_.resize(n_dfn);
                rhs_f_.resize(n_dfn);
                x_f_.resize(n_dfn);

                for (std::size_t i = 0; i < n_dfn; ++i) {
                    auto dfn = fracture_network_[i];
                    auto &dof_map_f = dfn->space().dof_map();

                    A_f_[i] = std::make_shared<Matrix>();
                    x_f_[i] = std::make_shared<Vector>();
                    rhs_f_[i] = std::make_shared<Vector>();

                    *x_f_[i] = local_zeros(dof_map_f.n_local_dofs());
                    dfn->assemble_flow(*x_f_[i], *A_f_[i], *rhs_f_[i]);

                    apply_boundary_conditions(dof_map_f, *A_f_[i], *rhs_f_[i]);
                    // apply_boundary_conditions(dfn->space(), *A_f_[i], *rhs_f_[i]);

                    if (report_) {
                        report_->add_stat(dfn->space(), *A_f_[i]);
                    }

                    if (remove_constrained_dofs_) {
                        remove_constrained_dofs(dof_map_f, *A_f_[i], *rhs_f_[i]);
                    }

                    auto T = lagrange_multiplier_[i]->transfer_matrix();

                    Matrix A_temp = transpose(*T) * (*A_f_[i]) * (*T);

                    // utopia::write("A.m", A);
                    A += A_temp;

                    // utopia::write("A_temp.m", *A_f_[i]);

                    // utopia::write("T.m", *T);

                    rhs += transpose(*T) * (*rhs_f_[i]);
                }

                // apply_boundary_conditions(dof_map_m, *A_m_, *rhs_m_);
                // apply_boundary_conditions(pourous_matrix_.space(), *A_m_, *rhs_m_);
                apply_boundary_conditions(pourous_matrix_.space(), A, rhs);

            } else  // if(assembly_strategy_ == "monolithic")
            {
                std::cerr << "[Error] in assemble_flow; assembly-strategy = " << assembly_strategy_
                          << " not supported yet" << std::endl;
                assert(false && "TODO");
                return false;
            }

            rename("A", A);
            for (auto p : matrix_processors_) {
                p->post_process(A);
            }

            if (report_) {
                report_->post_process(A);
                report_->add_stat(pourous_matrix_.space(), A, "condensed");
            }

            return true;
        }

        inline void disassemble_flow(const Vector &x) {
            if (assembly_strategy_ == "static-condensation") {
                *x_m_ = x;
                pourous_matrix_.disassemble_flow(*x_m_);

                const std::size_t n_dfn = fracture_network_.size();

                for (std::size_t i = 0; i < n_dfn; ++i) {
                    auto T = lagrange_multiplier_[i]->transfer_matrix();
                    *x_f_[i] = (*T) * x;

                    static const bool plot_residual = false;
                    if (plot_residual) {
                        Vector r;
                        fracture_network_[i]->residual(*A_f_[i], *rhs_f_[i], *x_f_[i], r);
                        lagrange_multiplier_[i]->remove_mass(r);
                        write("multiplier_" + std::to_string(i) + ".e", fracture_network_[i]->space(), r);
                    }

                    fracture_network_[i]->disassemble_flow(*x_f_[i]);
                }
            } else  // if(assembly_strategy_ == "monolithic")
            {
                std::cerr << "[Error] in disassemble_flow; assembly-strategy = " << assembly_strategy_
                          << " not supported yet" << std::endl;
                assert(false && "TODO");
            }
        }

        inline bool export_flow() {
            // write(pourous_matrix_.space().equation_system().name() + ".e",
            // pourous_matrix_.space(), *x_m_);

            pourous_matrix_.write(*x_m_);

            const std::size_t n_dfn = fracture_network_.size();

            for (std::size_t i = 0; i < n_dfn; ++i) {
                write(fracture_network_[i]->space().equation_system().name() + ".e",
                      fracture_network_[i]->space(),
                      *x_f_[i]);
            }

            // FIXME (check for rank inside save)
            if (report_ && x_m_->comm().rank() == 0) {
                report_->save();
            }

            return true;
        }

        FracturedPourousMedia(libMesh::Parallel::Communicator &comm)
            : comm_(comm),
              pourous_matrix_(comm),
              assembly_strategy_("static-condensation"),
              use_mg_(false),
              remove_constrained_dofs_(false),
              rescale_(1.0),
              adaptivity_(false),
              solver_ptr(std::make_shared<Factorization<Matrix, Vector>>()) {}

    private:
        libMesh::Parallel::Communicator &comm_;
        PourousMatrix pourous_matrix_;
        std::vector<std::shared_ptr<DiscreteFractureNetwork>> fracture_network_;
        std::vector<std::shared_ptr<LagrangeMultiplier>> lagrange_multiplier_;

        std::shared_ptr<Matrix> A_m_;
        std::shared_ptr<Vector> rhs_m_;
        std::shared_ptr<Vector> x_m_;

        std::vector<std::shared_ptr<Matrix>> A_f_;
        std::vector<std::shared_ptr<Vector>> rhs_f_;
        std::vector<std::shared_ptr<Vector>> x_f_;

        std::string assembly_strategy_;
        bool use_mg_;
        bool remove_constrained_dofs_;
        Scalar rescale_;

        std::vector<std::shared_ptr<MatrixPostProcessor<Matrix>>> matrix_processors_;
        std::shared_ptr<DFMReport<Matrix, Vector>> report_;

        bool adaptivity_;
        std::shared_ptr<LinearSolver<Matrix, Vector>> solver_ptr;
    };

}  // namespace utopia

#endif  // UTOPIA_FRACTURED_POUROUS_MEDIA_HPP
