#ifndef UTOPIA_M3ELINSOL_IMPL_HPP
#define UTOPIA_M3ELINSOL_IMPL_HPP

#include <vector>
#include "M3Elinsol_CXX.hpp"
#include "utopia_m3elinsol.hpp"

namespace utopia {

    class M3ELinSol::Impl {
    public:
        M3Elinsol_Slv solver;
        M3Elinsol_Mat mat;
        M3Elinsol_Vec rhs, sol;
        M3Elinsol_Int ierr;

        // params
        M3Elinsol_Int plevel;
        M3Elinsol_Int iverbosity;
        M3Elinsol_Str logfile;

        // mat buffs
        std::vector<M3Elinsol_Int> row_ptr;
        std::vector<M3Elinsol_Int> col_ind;
        std::vector<M3Elinsol_Real> values;

        // vec buffs
        std::vector<M3Elinsol_Real> sol_buff;
        std::vector<M3Elinsol_Real> rhs_buff;

        Impl(const M3Elinsol_Int plevel = 1,
             const M3Elinsol_Int iverbosity = 1,
             const M3Elinsol_Str logfile = "m3elinsol_output.log")
            : plevel(plevel), iverbosity(iverbosity), logfile(logfile) {
            ierr = M3Elinsol_Init(plevel, iverbosity, logfile, &solver);
            M3Elinsol_Errchk(&solver, ierr);
        }

        template <class T>
        void update_vec(const Wrapper<T, 1> &vec_in, M3Elinsol_Vec &vec_out, std::vector<M3Elinsol_Real> &buff_out) {
            auto ls = local_size(vec_in).get(0);
            auto r = range(vec_in);
            buff_out.resize(ls);

            each_read(vec_in, [&buff_out, &r](const M3Elinsol_Int i, const M3Elinsol_Real val) {
                buff_out[i - r.begin()] = val;
            });

            vec_out.values = &buff_out[0];
            vec_out.nrows = ls;
        }

        template <class T>
        void update_vecs(const Wrapper<T, 1> &rhs, const Wrapper<T, 1> &sol) {
            update_vec(rhs, this->rhs, this->rhs_buff);
            update_vec(sol, this->sol, this->sol_buff);
        }

        template <class T>
        void copy_buff(const std::vector<M3Elinsol_Real> &buff, Wrapper<T, 1> &vec) {
            auto r = range(vec);
            each_write(vec, [&buff, &r](const M3Elinsol_Int i) -> M3Elinsol_Real { return buff[i - r.begin()]; });
        }

        template <class T>
        bool apply(const Wrapper<T, 1> &rhs, Wrapper<T, 1> &sol) {
            update_vecs(rhs, sol);

            ierr = M3Elinsol_Solve(&this->solver, &this->mat, &this->rhs, &this->sol);
            M3Elinsol_Errchk(&solver, ierr);

            copy_buff(this->sol_buff, sol);
            return ierr == 0;
        }

        void printSystem(const bool binwrite, const M3Elinsol_Str systemfile) {
            ierr = M3Elinsol_DumpLinSys(&this->solver, binwrite, systemfile);
            M3Elinsol_Errchk(&solver, ierr);
        }

        template <class Matrix, class Vector, int Backend>
        void set_options_from_interface(const ASPAMG<Matrix, Vector, Backend> &interface) {
            // example ... fill the rest
            // M3Elinsol_ASPAMG_SetTspMaxit(&solver.handle, interface.max_it());
        }

        // default options
        void set_options() {
            // Set Problem level
            // M3Elinsol_ASPAMG_SetProblemLevel(&solver.handle, 1);

            // Set Preconditioner Options
            // M3Elinsol_ASPAMG_SetSmtNstep(&solver.handle, 5);
            // M3Elinsol_ASPAMG_SetSmtNu(&solver.handle, 1);
            // M3Elinsol_ASPAMG_SetTspMethod(&solver.handle, 1);
            // M3Elinsol_ASPAMG_SetTspNtvecs(&solver.handle, 5);
            // M3Elinsol_ASPAMG_SetTspMaxit(&solver.handle, 50);
            // M3Elinsol_ASPAMG_SetCsnTheta(&solver.handle, 0.5);
            // M3Elinsol_ASPAMG_SetPrlNnzrmax(&solver.handle, 4);
            // M3Elinsol_ASPAMG_SetAmgMaxlvls(&solver.handle, 10);
            // M3Elinsol_ASPAMG_SetAmgMaxcsize(&solver.handle, 20);

            // Set Solver Options
            // M3Elinsol_ASPAMG_SetPcgRelConvTol(&solver.handle, 1e-9);
            // M3Elinsol_ASPAMG_SetPcgMaxIter(&solver.handle, 800);
        }

        template <class T>
        void update(const Wrapper<T, 2> &mat_in) {
            auto ls = local_size(mat_in);
            auto n_row_local = ls.get(0);
            row_ptr.resize(n_row_local + 1);
            std::fill(row_ptr.begin(), row_ptr.end(), 0);
            row_ptr[0] = 1;

            auto r = row_range(mat_in);

            each_read(mat_in, [this, &r](const M3Elinsol_Int i, const M3Elinsol_Int, const M3Elinsol_Real) {
                ++row_ptr[i - r.begin() + 1];
            });

            for (std::size_t i = 1; i < row_ptr.size(); ++i) {
                row_ptr[i] += row_ptr[i - 1];
            }

            col_ind.clear();
            col_ind.reserve(row_ptr.back());

            values.clear();
            values.reserve(row_ptr.back());

            each_read(mat_in, [this](const M3Elinsol_Int i, const M3Elinsol_Int j, const M3Elinsol_Real val) {
                // Fortran offset +1
                col_ind.push_back(j + 1);
                values.push_back(val);
            });

            mat.rowptr = &row_ptr[0];
            mat.colind = &col_ind[0];
            mat.values = &values[0];
            mat.nrows = n_row_local;
            mat.nterm = col_ind.size();

            ierr = M3Elinsol_Set(&this->solver, &this->mat);
            M3Elinsol_Errchk(&solver, ierr);
        }

        ~Impl() {
            ierr = M3Elinsol_Destroy(&solver);
            M3Elinsol_Errchk(&solver, ierr);
        }
    };

    M3ELinSol::~M3ELinSol() {}

    M3ELinSol::M3ELinSol() : impl(std::make_shared<Impl>()) {}

    template <class Matrix, class Vector, int Backend>
    bool ASPAMG<Matrix, Vector, Backend>::apply(const Vector &b, Vector &x) {
        return solver.impl->apply(b, x);
    }

    template <class Matrix, class Vector, int Backend>
    void ASPAMG<Matrix, Vector, Backend>::update(const std::shared_ptr<const Matrix> &op) {
        IterativeSolver<Matrix, Vector>::update(op);
        // solver.impl->set_options_from_interface(*this);
        solver.impl->update(*op);
    }

    template <class Matrix, class Vector, int Backend>
    void ASPAMG<Matrix, Vector, Backend>::print_system(const bool binwrite, const M3Elinsol_Str systemfile) {
        return solver.impl->printSystem(binwrite, systemfile);
    }

    template <class Matrix, class Vector, int Backend>
    void ASPAMG<Matrix, Vector, Backend>::read(Input &in) {
        assert((solver.impl));
        auto &handle = solver.impl->solver.handle;

        // Set Problem level and Preconditioner Options
        // FIXME is there a way to read from handle what was the setting?

        int problemlevel = 1;
        in.get("ProblemLevel", problemlevel);
        M3Elinsol_ASPAMG_SetProblemLevel(&handle, problemlevel);

        int smtnstep = 5;
        in.get("SmtNstep", smtnstep);
        M3Elinsol_ASPAMG_SetSmtNstep(&handle, smtnstep);

        int smtnu = 1;
        in.get("SmtNu", smtnu);
        M3Elinsol_ASPAMG_SetSmtNu(&handle, smtnu);

        int tspmethod = 1;
        in.get("TspMethod", tspmethod);
        M3Elinsol_ASPAMG_SetTspMethod(&handle, tspmethod);

        int tspntvecs = 5;
        in.get("TspNtvecs", tspntvecs);
        M3Elinsol_ASPAMG_SetTspNtvecs(&handle, tspntvecs);

        int tspmaxit = 50;
        in.get("TspMaxit", tspmaxit);
        M3Elinsol_ASPAMG_SetTspMaxit(&handle, tspmaxit);

        double csntheta = 0.5;
        in.get("CsnTheta", csntheta);
        M3Elinsol_ASPAMG_SetCsnTheta(&handle, csntheta);

        int prlnnzrmax = 4;
        in.get("PrlNnzrmax", prlnnzrmax);
        M3Elinsol_ASPAMG_SetPrlNnzrmax(&handle, prlnnzrmax);

        int amgmaxlvls = 10;
        in.get("AmgMaxlvls", amgmaxlvls);
        M3Elinsol_ASPAMG_SetAmgMaxlvls(&handle, amgmaxlvls);

        int amgmaxcsize = 20;
        in.get("AmgMaxcsize", amgmaxcsize);
        M3Elinsol_ASPAMG_SetAmgMaxcsize(&handle, amgmaxcsize);

        // Set Solver Options
        double pcgrelconvtol = 1e-9;
        in.get("PcgRelConvTol", pcgrelconvtol);
        M3Elinsol_ASPAMG_SetPcgRelConvTol(&handle, pcgrelconvtol);

        int pcgmaxiter = 800;
        in.get("PcgMaxIter", pcgmaxiter);
        M3Elinsol_ASPAMG_SetPcgMaxIter(&handle, pcgmaxiter);
    }
}  // namespace utopia

#endif  // UTOPIA_M3ELINSOL_IMPL_HPP
