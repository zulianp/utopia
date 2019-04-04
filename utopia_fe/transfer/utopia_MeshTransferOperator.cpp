#include "utopia_MeshTransferOperator.hpp"

#include "utopia_L2LocalAssembler.hpp"
#include "utopia_ApproxL2LocalAssembler.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_Local2Global.hpp"
#include "utopia_BidirectionalL2LocalAssembler.hpp"
#include "utopia_MixedBidirectionalLocalAssembler.hpp"
#include "utopia_make_unique.hpp"

#include "libmesh/boundary_mesh.h"

#include <map>
#include <cmath>

namespace utopia {

    class MeshTransferOperator::Params : public Configurable {
    public:
        Params()
        : normalize_rows(true),
        tol(1e-14),
        bi_operator_mass_mat_outside(true),
        operator_type("l2-projection"),
        clamped(false),
        quad_order(-1),
        force_shell(false),
        write_operators_to_disk(false),
        force_zero_extension(false),
        from_boundary(false),
        to_boundary(false),
        normalize(false),
        use_composite_bidirectional(false),
        use_interpolation(false),
        nnz_x_row(0),
        output_path("./")
        {}

        void describe(std::ostream &os = std::cout) const
        {
            os << "normalize_rows " << normalize_rows << std::endl;
            os << "tol " << tol << std::endl;
            os << "bi_operator_mass_mat_outside " << bi_operator_mass_mat_outside << std::endl;
            os << "operator_type " << operator_type << std::endl;
            os << "clamped " << clamped << std::endl;
            os << "quad_order " << quad_order << std::endl;
            os << "force_shell " << force_shell << std::endl;
            os << "write_operators_to_disk " << write_operators_to_disk << std::endl;
            os << "force_zero_extension " << force_zero_extension << std::endl;
            os << "from_boundary " << from_boundary << std::endl;
            os << "to_boundary " << to_boundary << std::endl;
            os << "normalize " << normalize << std::endl;
            os << "use_composite_bidirectional " << use_composite_bidirectional << std::endl;
            os << "use_interpolation " << use_interpolation << std::endl;
            os << "nnz_x_row " << nnz_x_row << std::endl;
            os << "output_path " << output_path << std::endl;
        }

        inline void read(Input &is) override
        {
            is.get("type", operator_type);
            is.get("bi-op-mass-outside", bi_operator_mass_mat_outside);
            is.get("clamped", clamped);
            is.get("quad-order-approx", quad_order);

            is.get("force-shell", force_shell);
            is.get("write-operators-to-disk", write_operators_to_disk);
            is.get("force-zero-extension", force_zero_extension);
            // is.get("from-boundary", from_boundary);
            // is.get("to-boundary",   to_boundary);
            is.get("normalize", normalize);
            is.get("use-composite-bidirectional", use_composite_bidirectional);
            is.get("use-interpolation", use_interpolation);
            is.get("nnz-x-row", nnz_x_row);
            is.get("output-path", output_path);
        }

        bool normalize_rows;
        double tol;
        bool bi_operator_mass_mat_outside;
        std::string operator_type;
        bool clamped;
        int quad_order;
        bool force_shell;
        bool write_operators_to_disk;
        bool force_zero_extension;
        bool from_boundary;
        bool to_boundary;
        bool normalize;
        bool use_composite_bidirectional;
        bool use_interpolation;
        long nnz_x_row;
        std::string output_path;
    };

    static const std::map<std::string, TransferOperatorType> &get_str_to_type()
    {
        static std::map<std::string, TransferOperatorType> types;

        if(types.empty()) {
            types["INTERPOLATION"] 		  = INTERPOLATION;
            types["L2_PROJECTION"] 		  = L2_PROJECTION;
            types["PSEUDO_L2_PROJECTION"] = PSEUDO_L2_PROJECTION;
            types["APPROX_L2_PROJECTION"] = APPROX_L2_PROJECTION;
            types["BIDIRECTIONAL_L2_PROJECTION"] = BIDIRECTIONAL_L2_PROJECTION;
            types["BIDIRECTIONAL_PSEUDO_L2_PROJECTION"] = BIDIRECTIONAL_PSEUDO_L2_PROJECTION;

            //other way of writing them
            types["interpolation"] 		  = INTERPOLATION;
            types["l2-projection"] 		  = L2_PROJECTION;
            types["pseudo-l2-projection"] = PSEUDO_L2_PROJECTION;
            types["approx-l2-projection"] = APPROX_L2_PROJECTION;
            types["bidirectional-l2-projection"] = BIDIRECTIONAL_L2_PROJECTION;
            types["bidirectional-pseudo-l2-projection"] = BIDIRECTIONAL_PSEUDO_L2_PROJECTION;
        }

        return types;
    }

    MeshTransferOperator::MeshTransferOperator(
        const std::shared_ptr<MeshBase> &from_mesh,
        const std::shared_ptr<DofMap>   &from_dofs,
        const std::shared_ptr<MeshBase> &to_mesh,
        const std::shared_ptr<DofMap>   &to_dofs,
        const TransferOptions &opts
        ) :
    from_mesh(from_mesh),
    from_dofs(from_dofs),
    to_mesh(to_mesh),
    to_dofs(to_dofs),
    opts(opts),
    params_(utopia::make_unique<Params>())
    {
        assembly_strategies_[INTERPOLATION] =
        [this]() -> bool { return this->set_up_interpolation(); };

        assembly_strategies_[L2_PROJECTION] =
        [this]() -> bool { return this->set_up_l2_projection(); };

        assembly_strategies_[PSEUDO_L2_PROJECTION] =
        [this]() -> bool { return this->set_up_pseudo_l2_projection(); };

        assembly_strategies_[APPROX_L2_PROJECTION] =
        [this]() -> bool { return this->set_up_approx_l2_projection(); };

        assembly_strategies_[BIDIRECTIONAL_L2_PROJECTION] =
        [this]() -> bool { return this->set_up_bidirectional_transfer(); };

        assembly_strategies_[BIDIRECTIONAL_PSEUDO_L2_PROJECTION] =
        [this]() -> bool { return this->set_up_bidirectional_pseudo_transfer(); };
    }

    MeshTransferOperator::~MeshTransferOperator() {}

    void MeshTransferOperator::read(Input &is)
    {
        params_->read(is);
    }

    void MeshTransferOperator::assemble_mass_matrix(
        const libMesh::MeshBase &mesh,
        const libMesh::DofMap &dof_map,
        const int var,
        const int n_tensor,
        USparseMatrix &mat
        ) const
    {
        Chrono c;
        c.start();

        auto dim = mesh.mesh_dimension();
        auto fe_type = dof_map.variable_type(var);
        std::vector<libMesh::dof_id_type> indices;

        SizeType nnz_x_row = 0;

        if(params_->nnz_x_row <= 0 &&
            !dof_map.get_n_nz().empty()) {
            nnz_x_row =
            *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end());

            if(!dof_map.get_n_oz().empty()) {
                nnz_x_row += *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end());
            }
        } else {
            m_utopia_warning_once("MeshTransferOperator::assemble_mass_matrix(...) using user nnz_x_row");
            nnz_x_row = params_->nnz_x_row;
        }

        // std::cout << "nnz_x_row " << nnz_x_row << std::endl;

        assert(nnz_x_row != 0 && "super innefficient");
        mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);

        {
            Write<USparseMatrix> w_m(mat, utopia::GLOBAL_ADD);
            libMesh::DenseMatrix<libMesh::Real> el_mat;

            auto fe = libMesh::FEBase::build(dim, fe_type);
            libMesh::QGauss qrule(dim, fe_type.default_quadrature_order());
            fe->attach_quadrature_rule(&qrule);

            auto &phi = fe->get_phi();
            auto &JxW = fe->get_JxW();

            for(auto it = elements_begin(mesh); it != elements_end(mesh); ++it) {
                fe->reinit(*it);
                dof_map.dof_indices(*it, indices);

                auto n_shape_functions = phi.size();
                auto n_qp = qrule.n_points();
                auto n_indices = indices.size();

                assert(n_qp > 0);
                assert(n_shape_functions > 0);

                if(n_shape_functions * n_tensor != n_indices) {
                    m_utopia_warning_once("dof_map not consistent with tensorization of variable")
                }

                el_mat.resize(
                    n_indices,
                    n_indices
                );

                el_mat.zero();


                for(unsigned int i = 0; i < n_shape_functions; i++) {
                    for(unsigned int j = 0; j < n_shape_functions; j++) {
                        for(unsigned int qp = 0; qp < n_qp; qp++) {
                            auto value = JxW[qp] * phi[i][qp] * phi[j][qp];
                            assert(!std::isnan(value) && !std::isinf(value));

                            for(unsigned int k = 0; k < n_tensor; k++) {
                                el_mat(i + k * phi.size(), j + k * phi.size()) += value;
                            }
                        }
                    }
                }

                add_matrix(el_mat, indices, indices, mat);
            }
        }

        c.stop();
        std::cout << "assemble_mass_matrix (time)\n" << c << std::endl;
    }

    bool MeshTransferOperator::set_up_l2_projection()
    {
        std::cout << "[Status] using l2 projection" << std::endl;
        auto from = get_filtered_from_mesh();
        bool is_shell = params_->force_shell || from->mesh_dimension() < from->spatial_dimension();
        auto assembler = std::make_shared<L2LocalAssembler>(from->mesh_dimension(), false, true, is_shell);
        auto local2global = std::make_shared<Local2Global>(false);

        std::vector< std::shared_ptr<SparseMatrix> > mats;
        TransferAssembler transfer_assembler(assembler, local2global);
        if(!transfer_assembler.assemble(from, from_dofs, get_filtered_to_mesh(), to_dofs, mats, opts)) {
            return false;
        }

        auto l2_operator = std::make_shared<L2TransferOperator>(mats[0], mats[1], new_solver());
        l2_operator->fix_mass_matrix_operator(params_->tol);
        l2_operator->init();
        operator_ = l2_operator;
        return true;
    }

    bool MeshTransferOperator::set_up_interpolation()
    {
        std::cout << "[Status] using interpolation" << std::endl;
        auto assembler = std::make_shared<InterpolationLocalAssembler>(get_filtered_from_mesh()->mesh_dimension());
        auto local2global = std::make_shared<Local2Global>(true);

        std::vector< std::shared_ptr<SparseMatrix> > mats;
        TransferAssembler transfer_assembler(assembler, local2global);
        if(!transfer_assembler.assemble(get_filtered_from_mesh(), from_dofs, get_filtered_to_mesh(), to_dofs, mats, opts)) {
            return false;
        }

        auto interpolation_operator = std::make_shared<Interpolator>(mats[0]);
        //only necessary for non-conforming mesh in parallel
        interpolation_operator->normalize_rows();
        operator_ = interpolation_operator;
        return true;
    }

    bool MeshTransferOperator::set_up_approx_l2_projection()
    {
        std::cout << "[Status] using approx l2 projection" << std::endl;
        auto assembler = std::make_shared<ApproxL2LocalAssembler>(get_filtered_from_mesh()->mesh_dimension());
        assembler->set_quadrature_order(params_->quad_order);
        auto local2global = std::make_shared<Local2Global>(false);

        std::vector< std::shared_ptr<SparseMatrix> > mats;
        TransferAssembler transfer_assembler(assembler, local2global);
        if(!transfer_assembler.assemble(get_filtered_from_mesh(), from_dofs, get_filtered_to_mesh(), to_dofs, mats, opts)) {
            return false;
        }

        auto l2_operator = std::make_shared<L2TransferOperator>(mats[0], mats[1], new_solver());
        l2_operator->fix_mass_matrix_operator(params_->tol);
        l2_operator->init();
        operator_ = l2_operator;
        return true;
    }

    bool MeshTransferOperator::set_up_pseudo_l2_projection()
    {
        std::cout << "[Status] using pseudo l2 projection" << std::endl;

        auto from = get_filtered_from_mesh();
        bool is_shell = params_->force_shell || from->mesh_dimension() < from->spatial_dimension();

        auto assembler = std::make_shared<L2LocalAssembler>(from->mesh_dimension(), true, false);
        auto local2global = std::make_shared<Local2Global>(false);

        std::vector< std::shared_ptr<SparseMatrix> > mats;
        TransferAssembler transfer_assembler(assembler, local2global);
        if(!transfer_assembler.assemble(from, from_dofs, get_filtered_to_mesh(), to_dofs, mats, opts)) {
            return false;
        }

        if(params_->normalize_rows) {
            auto pseudo_l2_operator = std::make_shared<PseudoL2TransferOperator>();
            pseudo_l2_operator->init_from_coupling_operator(*mats[0]);
            operator_ = pseudo_l2_operator;
        } else {
            operator_ = std::make_shared<PseudoL2TransferOperator>(mats[0]);
        }

        return true;
    }

    bool MeshTransferOperator::set_up_bidirectional_transfer()
    {
        std::cout << "[Status] using bi l2 projection" << std::endl;
        std::shared_ptr<LocalAssembler> assembler;

        bool use_interpolation = false;
        const bool assemble_D =
            !params_->bi_operator_mass_mat_outside ||
             params_->use_composite_bidirectional;


        // params_->print();

        if(params_->use_composite_bidirectional) {
            auto from = get_filtered_from_mesh();

            if(params_->use_interpolation) {
                auto ass_1 = std::make_shared<InterpolationLocalAssembler>(from->mesh_dimension());
                auto ass_2 = std::make_shared<InterpolationLocalAssembler>(from->mesh_dimension());
                use_interpolation = true;
                auto ass = std::make_shared<MixedBidirectionalLocalAssembler>(ass_1, ass_2);
                ass->set_return_false_only_if_both_fail(true);
                assembler = ass;
            } else {
                auto ass_1 = std::make_shared<L2LocalAssembler>(from->mesh_dimension(), false, true);
                auto ass_2 = std::make_shared<ApproxL2LocalAssembler>(from->mesh_dimension());
                assembler = std::make_shared<MixedBidirectionalLocalAssembler>(ass_1, ass_2);
            }

        } else {

            std::cout << "params_->bi_operator_mass_mat_outside " << (params_->bi_operator_mass_mat_outside) << std::endl;
            assembler = std::make_shared<BidirectionalL2LocalAssembler>(
                get_filtered_from_mesh()->mesh_dimension(),
                false,
                !params_->bi_operator_mass_mat_outside);

            std::cout << "n_forms: " << assembler->n_forms() << std::endl;
        }

        auto local2global = std::make_shared<Local2Global>(use_interpolation);

        std::vector< std::shared_ptr<SparseMatrix> > mats;
        TransferAssembler transfer_assembler(assembler, local2global);
        if(!transfer_assembler.assemble(
            get_filtered_from_mesh(),
            from_dofs,
            get_filtered_to_mesh(),
            to_dofs,
            mats,
            opts)) {
            return false;
        }

        if(!assemble_D && !use_interpolation) {
            auto mass_mat_from = std::make_shared<USparseMatrix>();
            auto mass_mat_to   = std::make_shared<USparseMatrix>();

            assemble_mass_matrix(*get_filtered_from_mesh(), *from_dofs, opts.from_var_num, opts.n_var, *mass_mat_from);
            assemble_mass_matrix(*get_filtered_to_mesh(),   *to_dofs,   opts.to_var_num,   opts.n_var, *mass_mat_to);

            const bool restrict_mass_matrix = !params_->use_composite_bidirectional;
            auto forward = std::make_shared<L2TransferOperator>(mats[0], mass_mat_to, new_solver());

            if(!restrict_mass_matrix) {
                forward->fix_mass_matrix_operator(params_->tol);
            } else {
                // forward->restrict_mass_matrix(params_->tol);
                forward->restrict_mass_matrix_old(params_->tol);
            }

            auto backward = std::make_shared<L2TransferOperator>(mats[1], mass_mat_from, new_solver());

            if(!restrict_mass_matrix) {
                backward->fix_mass_matrix_operator(params_->tol);
            } else {
                // backward->restrict_mass_matrix(params_->tol);
                backward->restrict_mass_matrix_old(params_->tol);
            }

            forward->init();
            backward->init();

            if(params_->normalize) {
                operator_ = std::make_shared<BidirectionalOperator>(
                    std::make_shared<NormalizedOperator>(
                        local_size(*forward),
                        forward
                    ),
                    std::make_shared<NormalizedOperator>(
                        local_size(*backward),
                        backward
                    )
                );
            } else {
                operator_ = std::make_shared<BidirectionalOperator>(forward, backward);
            }

        } else {
            if(use_interpolation) {
                operator_ = std::make_shared<BidirectionalOperator>(
                    std::make_shared<Interpolator>(mats[0]),
                    std::make_shared<Interpolator>(mats[1])
                );
            } else {
                auto forward = std::make_shared<L2TransferOperator>(mats[0], mats[1], new_solver());
                forward->fix_mass_matrix_operator(params_->tol);

                auto backward = std::make_shared<L2TransferOperator>(mats[2], mats[3], new_solver());
                backward->fix_mass_matrix_operator(params_->tol);

                forward->init();
                backward->init();
                operator_ = std::make_shared<BidirectionalOperator>(forward, backward);

                if(params_->normalize) {
                    operator_ = std::make_shared<BidirectionalOperator>(
                        std::make_shared<NormalizedOperator>(
                            local_size(*forward),
                            forward
                        ),
                        std::make_shared<NormalizedOperator>(
                            local_size(*backward),
                            backward
                        )
                    );
                } else {
                    operator_ = std::make_shared<BidirectionalOperator>(forward, backward);
                }
            }
        }

        return true;
    }

    bool MeshTransferOperator::set_up_bidirectional_pseudo_transfer()
    {
        std::cout << "[Status] using bi pseudo l2 projection" << std::endl;
        auto assembler = std::make_shared<BidirectionalL2LocalAssembler>(
            get_filtered_from_mesh()->mesh_dimension(),
            true,
            false
        );

        auto local2global = std::make_shared<Local2Global>(false);

        std::vector< std::shared_ptr<SparseMatrix> > mats;
        TransferAssembler transfer_assembler(assembler, local2global);
        if(!transfer_assembler.assemble(get_filtered_from_mesh(), from_dofs, get_filtered_to_mesh(), to_dofs, mats, opts)) {
            return false;
        }

        auto mass_mat_from = std::make_shared<USparseMatrix>();
        auto mass_mat_to   = std::make_shared<USparseMatrix>();

        assemble_mass_matrix(*get_filtered_from_mesh(), *from_dofs, opts.from_var_num, opts.n_var, *mass_mat_from);
        assemble_mass_matrix(*get_filtered_to_mesh(),   *to_dofs,   opts.to_var_num,   opts.n_var, *mass_mat_to);

        auto forward = std::make_shared<PseudoL2TransferOperator>();
        forward->init_from_coupling_and_mass_operator(*mats[0], *mass_mat_to);

        auto backward = std::make_shared<PseudoL2TransferOperator>();
        backward->init_from_coupling_and_mass_operator(*mats[1], *mass_mat_from);
        operator_ = std::make_shared<BidirectionalOperator>(forward, backward);
        return true;
    }

    bool MeshTransferOperator::assemble()
    {
        return initialize(params_->operator_type);
    }

    bool MeshTransferOperator::initialize(const std::string operator_type)
    {
        const auto &m = get_str_to_type();
        params_->operator_type = operator_type;

        auto it = m.find(operator_type);

        if(it == m.end()) {
            return initialize(PSEUDO_L2_PROJECTION);
        } else {
            return initialize(it->second);
        }
    }

    static bool is_bidirectional(const TransferOperatorType operator_type) {
        return operator_type == BIDIRECTIONAL_L2_PROJECTION ||
               operator_type == BIDIRECTIONAL_PSEUDO_L2_PROJECTION;
    }

    bool MeshTransferOperator::initialize(const TransferOperatorType operator_type)
    {
        //apply boundary filter
        if(params_->from_boundary) {
            auto b_from_mesh = std::make_shared<libMesh::BoundaryMesh>(from_mesh->comm(), from_mesh->mesh_dimension()-1);
            from_mesh->boundary_info->sync(*b_from_mesh);
            filtered_from_mesh = b_from_mesh;
        }

        if(params_->to_boundary) {
            auto b_to_mesh = std::make_shared<libMesh::BoundaryMesh>(to_mesh->comm(), to_mesh->mesh_dimension()-1);
            to_mesh->boundary_info->sync(*b_to_mesh);
            filtered_to_mesh = b_to_mesh;
        }

        ///////////////////////////////////////////////////
        auto it = assembly_strategies_.find(operator_type);

        bool ok = false;
        if(it == assembly_strategies_.end()) {
            ok = assembly_strategies_.begin()->second();
        } else {
            ok = it->second();
        }
        ////////////////////////////////////////////////////

        if(!ok) {
            std::cerr << "[Error] assembly failed" << std::endl;
            return false;
        }

        if(params_->write_operators_to_disk) {
            operator_->write(params_->output_path);
        }

        if(params_->normalize && !is_bidirectional(operator_type)) {
            operator_ = std::make_shared<NormalizedOperator>(
                Size({
                        to_dofs->n_local_dofs(),
                        from_dofs->n_local_dofs()
                    }),
                operator_
            );
        }

        if(params_->force_zero_extension) {
            operator_ = std::make_shared<ForceZeroExtension>(operator_, params_->tol);
        }

        if(params_->clamped) {
            operator_ = std::make_shared<ClampedOperator>(operator_);
        }

        operator_->describe(std::cout);
        return true;
    }

    void MeshTransferOperator::set_tol(const double val)
    {
        params_->tol = val;
    }

    std::shared_ptr<libMesh::MeshBase> MeshTransferOperator::get_filtered_from_mesh()
    {
        if(filtered_from_mesh) {
            return filtered_from_mesh;
        }

        return from_mesh;
    }

    std::shared_ptr<libMesh::MeshBase> MeshTransferOperator::get_filtered_to_mesh()
    {
        if(filtered_to_mesh) {
            return filtered_to_mesh;
        }

        return to_mesh;
    }

    std::unique_ptr<LinearSolver<USparseMatrix, UVector> > MeshTransferOperator::new_solver() {
         return utopia::make_unique<ConjugateGradient<USparseMatrix, UVector>>();
        //return utopia::make_unique<GMRES<USparseMatrix, UVector>>("bjacobi");
        // return utopia::make_unique<Factorization<USparseMatrix, UVector>>("superlu_dist", "lu");
    }
}
