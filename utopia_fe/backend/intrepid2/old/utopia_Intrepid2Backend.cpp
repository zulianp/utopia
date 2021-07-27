#include "utopia_Intrepid2Backend.hpp"

namespace utopia {

    void apply_boundary_conditions(libMesh::DofMap &dof_map, TSUSerialMatrix &mat, TUSerialVector &vec) {
        if (utopia::Utopia::instance().verbose()) {
            std::cout << "apply_boundary_conditions begin: " << std::endl;
        }

        Chrono c;
        c.start();

        assert(!empty(mat));
        assert(!empty(vec));

        using SizeType = Traits<TUSerialVector>::SizeType;

        const bool has_constraints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

        Size ls = local_size(mat);
        Size s = size(mat);

        std::vector<SizeType> index;

        Range rr = range(vec);

        if (has_constraints) {
            for (SizeType i = rr.begin(); i < rr.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    index.push_back(i);
                }
            }
        }

        set_zero_rows(mat, index, 1.);

        Write<TUSerialVector> w_v(vec);

        if (has_constraints) {
            libMesh::DofConstraintValueMap &rhs_values = dof_map.get_primal_constraint_values();

            Range r = range(vec);
            for (SizeType i = r.begin(); i < r.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    auto valpos = rhs_values.find(i);
                    vec.set(i, (valpos == rhs_values.end()) ? 0 : valpos->second);
                }
            }
        }

        c.stop();

        if (utopia::Utopia::instance().verbose()) {
            std::cout << "apply_boundary_conditions end: " << c << std::endl;
        }
    }

    void apply_boundary_conditions(libMesh::DofMap &dof_map, TUSerialVector &vec) {
        const bool has_constraints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

        Write<TUSerialVector> w_v(vec);

        if (has_constraints) {
            libMesh::DofConstraintValueMap &rhs_values = dof_map.get_primal_constraint_values();

            Range r = range(vec);
            for (SizeType i = r.begin(); i < r.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    auto valpos = rhs_values.find(i);
                    vec.set(i, (valpos == rhs_values.end()) ? 0 : valpos->second);
                }
            }
        }
    }

}  // namespace utopia
