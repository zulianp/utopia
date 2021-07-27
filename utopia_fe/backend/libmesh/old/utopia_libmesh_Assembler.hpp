#ifndef UTOPIA_LIBMESH_ASSEMBLER_HPP
#define UTOPIA_LIBMESH_ASSEMBLER_HPP

#include "utopia_Adaptivity.hpp"
#include "utopia_libmesh_AssemblyContext.hpp"

namespace utopia {
    class LibMeshAssembler {
    public:
        typedef utopia::USparseMatrix GlobalMatrix;
        typedef utopia::UVector GlobalVector;
        typedef UTOPIA_SCALAR(GlobalVector) Scalar;

        LibMeshAssembler() : verbose_(Utopia::instance().verbose()) {}

        // FIXME put in utopia
        template <class T>
        static bool is_ghosted(const Tensor<T, 1> &vec) {
            return vec.derived().has_ghosts();
        }

        template <class Expr>
        bool assemble(/*const*/ Expr &expr, Scalar &val) {
            // perf
            Chrono c;
            c.start();

            typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
            typedef typename TraitsT::Matrix ElementMatrix;
            typedef typename TraitsT::Vector ElementVector;

            static const int Backend = TraitsT::Backend;

            const auto &space = find_space<LibMeshFunctionSpace>(expr);
            const auto &dof_map = space.dof_map();
            auto &m = space.mesh();

            val = 0.;

            auto e_begin = elements_begin(m);
            auto e_end = elements_end(m);

            if (e_begin != e_end) {
                init_context_on(expr, (*e_begin)->id());
                for (auto it = e_begin; it != e_end; ++it) {
                    if (it != e_begin) {
                        reinit_context_on(expr, (*it)->id());
                    }

                    Number<Scalar> el_val = 0.;

                    FormEvaluator<LIBMESH_TAG> eval;
                    eval.eval(expr, el_val, ctx_);

                    if (ctx_.has_assembled()) {
                        val += el_val;
                    }
                }
            }

            // FIXME (once libmesh is fixed go back to previous version)
            UVector dump;
            val = dump.comm().sum(val);
            // m.comm().sum(val);

            // perf
            c.stop();

            if (verbose_) {
                std::cout << "assemble: value" << std::endl;
                std::cout << c << std::endl;
            }

            return false;
        }

        template <class Expr>
        bool assemble(/*const*/ Expr &expr,
                      GlobalMatrix &mat,
                      GlobalVector &vec,
                      const bool apply_constraints = false) {
            // perf
            Chrono c;
            c.start();

            const bool disable_adaptivity = utopia::Utopia::instance().get("disable-adaptivity") == "true";

            typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
            typedef typename TraitsT::Matrix ElementMatrix;
            typedef typename TraitsT::Vector ElementVector;

            static const int Backend = TraitsT::Backend;

            const auto &space = find_space<LibMeshFunctionSpace>(expr);
            const auto &dof_map = space.dof_map();
            auto &m = space.mesh();

            auto s_m = size(mat);

            // FIXME trilinos backend is buggy
            if (Traits<GlobalMatrix>::Backend == utopia::TRILINOS || empty(mat) || s_m.get(0) != dof_map.n_dofs() ||
                s_m.get(1) != dof_map.n_dofs()) {
                auto nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
                                          *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));

                mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
            } else {
                if (mat.is_assembled()) {
                    mat *= 0.;
                }
            }

            typename Traits<GlobalVector>::IndexSet ghost_nodes;
            convert(dof_map.get_send_list(), ghost_nodes);
            GlobalVector temp_vec = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);
            // if(empty(vec) || size(vec).get(0) != dof_map.n_dofs() || !is_ghosted(vec)) {
            // vec = local_zeros(dof_map.n_local_dofs());

            // }
            // else {
            // 	vec.set(0.);
            // }

            libMesh::DofConstraints constraints;

            if (!disable_adaptivity) {
                Adaptivity::compute_all_constraints(m, dof_map, constraints);
            }

            {
                Write<GlobalMatrix> w_m(mat, utopia::GLOBAL_ADD);
                Write<GlobalVector> w_v(temp_vec, utopia::GLOBAL_ADD);

                ElementMatrix el_mat;
                ElementVector el_vec;

                auto e_begin = elements_begin(m);
                auto e_end = elements_end(m);

                if (e_begin != e_end) {
                    init_context_on(expr, (*e_begin)->id());
                }

                std::vector<libMesh::dof_id_type> dof_indices;
                for (auto it = e_begin; it != e_end; ++it) {
                    if (it != e_begin) {
                        reinit_context_on(expr, (*it)->id());
                    }

                    el_mat.set(0.0);
                    el_vec.set(0.0);

                    FormEvaluator<LIBMESH_TAG> eval;
                    eval.eval(expr, el_mat, el_vec, ctx_);

                    dof_map.dof_indices(*it, dof_indices);

                    if (ctx_.has_assembled()) {
                        if (apply_constraints) {
                            // std::cout<<"I am here heterogenously_constrain_element_matrix_and_vector"<<std::endl;
                            // dof_map.heterogenously_constrain_element_matrix_and_vector(
                            //     el_mat,
                            //     el_vec,
                            //     dof_indices
                            // );

                            assert(false);

                        } else {
                            // std::cout<<"Adaptivity::constrain_matrix_and_vector"<<std::endl;

                            if (!disable_adaptivity) {
                                Adaptivity::constrain_matrix_and_vector(
                                    *it, dof_map, constraints, el_mat, el_vec, dof_indices);
                            }
                        }

                        libMesh::Elem *ele = *it;
                        // std::cout<<"current_elem: "<<ele[0]<<std::endl;
                        // utopia::disp("el_mat");
                        // utopia::disp(el_mat);
                        // utopia::disp("el_vec");
                        // utopia::disp(el_vec);

                        add_matrix(el_mat, dof_indices, dof_indices, mat);
                        add_vector(el_vec, dof_indices, temp_vec);
                    }
                }
            }

            if (Traits<GlobalVector>::Backend == utopia::TRILINOS) {
                vec = 1. * temp_vec;  // avoid copying
            } else {
                vec = std::move(temp_vec);
            }

            // perf
            c.stop();
            if (verbose_) {
                std::cout << "assemble: lhs == rhs" << std::endl;
                std::cout << c << std::endl;
            }

            return true;
        }

        static void allocate_matrix(const libMesh::DofMap &dof_map, GlobalMatrix &mat) {
            auto s_m = size(mat);

            if (Traits<GlobalMatrix>::Backend == utopia::TRILINOS || empty(mat) || s_m.get(0) != dof_map.n_dofs() ||
                s_m.get(1) != dof_map.n_dofs()) {
                SizeType nnz_x_row = 0;
                if (!dof_map.get_n_nz().empty()) {
                    // nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
                    //  *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));

                    nnz_x_row = *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()) +
                                *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end());
                }

                mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
            } else {
                if (mat.is_assembled()) {
                    mat *= 0.;
                }
            }
        }

        template <class Expr>
        bool assemble(/*const*/ Expr &expr, GlobalMatrix &mat) {
            const bool disable_adaptivity = utopia::Utopia::instance().get("disable-adaptivity") == "true";

            // perf
            Chrono c;
            c.start();

            typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
            typedef typename TraitsT::Matrix ElementMatrix;

            static const int Backend = TraitsT::Backend;

            const auto &space = find_space<LibMeshFunctionSpace>(expr);
            const auto &dof_map = space.dof_map();
            auto &m = space.mesh();

            auto s_m = size(mat);

            libMesh::DofConstraints constraints;

            if (!disable_adaptivity) {
                Adaptivity::compute_all_constraints(m, dof_map, constraints);
            }

            // FIXME trilinos backend is buggy
            // if(Traits<GlobalMatrix>::Backend == utopia::TRILINOS || empty(mat) || s_m.get(0) != dof_map.n_dofs() ||
            // s_m.get(1) != dof_map.n_dofs()) {
            //     SizeType nnz_x_row = 0;
            //     if(!dof_map.get_n_nz().empty()) {
            //         // nnz_x_row = std::max(*std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()),
            //         // 	*std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end()));

            //         nnz_x_row =
            //             *std::max_element(dof_map.get_n_nz().begin(), dof_map.get_n_nz().end()) +
            //             *std::max_element(dof_map.get_n_oz().begin(), dof_map.get_n_oz().end());
            //     }

            //     mat = local_sparse(dof_map.n_local_dofs(), dof_map.n_local_dofs(), nnz_x_row);
            // } else {
            //     mat *= 0.;
            // }

            allocate_matrix(dof_map, mat);

            {
                Write<GlobalMatrix> w_m(mat, utopia::GLOBAL_ADD);

                auto e_begin = elements_begin(m);
                auto e_end = elements_end(m);

                if (e_begin != e_end) {
                    ElementMatrix el_mat;
                    init_context_on(expr, (*e_begin)->id());

                    for (auto it = e_begin; it != e_end; ++it) {
                        if (it != e_begin) {
                            reinit_context_on(expr, (*it)->id());
                        }

                        el_mat.set(0.0);

                        FormEvaluator<LIBMESH_TAG> eval;
                        eval.eval(expr, el_mat, ctx_, true);

                        std::vector<libMesh::dof_id_type> dof_indices;
                        dof_map.dof_indices(*it, dof_indices);

                        if (!disable_adaptivity) {
                            Adaptivity::constrain_matrix(*it, dof_map, constraints, el_mat, dof_indices);
                        }

                        if (ctx_.has_assembled()) {
                            add_matrix(el_mat, dof_indices, dof_indices, mat);
                        }
                    }
                }
            }

            // perf
            c.stop();

            if (verbose_) {
                std::cout << "assemble: lhs" << std::endl;
                std::cout << c << std::endl;
            }

            return true;
        }

        template <class Expr>
        bool assemble(/*const*/ Expr &expr, GlobalVector &vec, const bool apply_constraints = false) {
            const bool disable_adaptivity = utopia::Utopia::instance().get("disable-adaptivity") == "true";

            // perf
            Chrono c;
            c.start();

            typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
            typedef typename TraitsT::Vector ElementVector;

            static const int Backend = TraitsT::Backend;

            const auto &space = find_space<LibMeshFunctionSpace>(expr);
            const auto &dof_map = space.dof_map();
            auto &m = space.mesh();

            typename Traits<GlobalVector>::IndexSet ghost_nodes;

            convert(dof_map.get_send_list(), ghost_nodes);
            GlobalVector temp_vec = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);

            // if(empty(vec) || size(vec).get(0) != dof_map.n_dofs() || !is_ghosted(vec)) {
            // 	// vec = local_zeros(dof_map.n_local_dofs());
            // 	vec = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), dof_map.get_send_list());
            // } else {
            // 	vec *= 0.;
            // }

            libMesh::DofConstraints constraints;

            if (!disable_adaptivity) {
                Adaptivity::compute_all_constraints(m, dof_map, constraints);
            }

            {
                Write<GlobalVector> w_v(temp_vec, utopia::GLOBAL_ADD);
                ElementVector el_vec;

                auto e_begin = elements_begin(m);
                auto e_end = elements_end(m);

                if (e_begin != e_end) {
                    init_context_on(expr, (*e_begin)->id());

                    for (auto it = e_begin; it != e_end; ++it) {
                        if (it != e_begin) {
                            reinit_context_on(expr, (*it)->id());
                        }

                        el_vec.set(0.0);

                        FormEvaluator<LIBMESH_TAG> eval;
                        eval.eval(expr, el_vec, ctx_, true);

                        std::vector<libMesh::dof_id_type> dof_indices;
                        dof_map.dof_indices(*it, dof_indices);

                        if (!disable_adaptivity) {
                            Adaptivity::constrain_vector(*it, dof_map, constraints, el_vec, dof_indices);
                        }

                        if (ctx_.has_assembled()) {
                            add_vector(el_vec, dof_indices, temp_vec);
                        }
                    }
                }
            }

            if (Traits<GlobalVector>::Backend == utopia::TRILINOS) {
                vec = 1. * temp_vec;
            } else {
                vec = std::move(temp_vec);
            }

            // perf
            c.stop();

            if (verbose_) {
                std::cout << "assemble: rhs" << std::endl;
                std::cout << c << std::endl;
            }

            return true;
        }

        template <class Expr>
        void init_context_on(const Expr &expr, const int element_id) {
            ctx_.set_current_element(element_id);
            ctx_.set_has_assembled(false);
            ctx_.init(expr);
        }

        template <class Expr>
        void reinit_context_on(const Expr &expr, const int element_id) {
            ctx_.set_current_element(element_id);
            ctx_.set_has_assembled(false);
            ctx_.reinit(expr);
        }

    private:
        AssemblyContext<LIBMESH_TAG> ctx_;
        bool verbose_;
    };
}  // namespace utopia

#endif  // UTOPIA_LIBMESH_ASSEMBLER_HPP
