#ifndef UTOPIA_MULTILEVEL_TEST_PROBLEM_HPP
#define UTOPIA_MULTILEVEL_TEST_PROBLEM_HPP

#include "utopia_assemble_laplacian_1D.hpp"

namespace utopia {
    template<class Matrix, class Vector>
    class MultiLevelTestProblem {
    public:
        typedef UTOPIA_SCALAR(Vector) Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        MultiLevelTestProblem(
            const SizeType n_coarse_elements,
            const SizeType n_levels,
            const bool remove_bc_dofs_from_interp = false)
        {

            assert(n_coarse_elements > 0);
            assert(n_levels > 1);

            n_dofs.resize(n_levels);
            n_dofs[0] = n_coarse_elements + 1;

            for(SizeType i = 1; i < n_levels; ++i) {
                n_dofs[i] = (n_dofs[i-1] - 1) * 2 + 1;
            }

            // Scalar h = 1.;//1./(n_dofs[n_levels -1]);
            Scalar h = 1./(n_dofs[n_levels -1] - 1);

            interpolators.resize(n_levels - 1);

            for(SizeType i = 0; i < n_levels - 1; ++i) {
                const auto n_coarse = n_dofs[i];
                const auto n_fine   = n_dofs[i + 1];
                interpolators[i] = std::make_shared<Matrix>(sparse(n_fine, n_coarse, 2));
                auto &I = *interpolators[i];

                Write<Matrix> w_(I, utopia::GLOBAL_INSERT);
                auto r = row_range(I);

                SizeType j = std::floor(r.begin()/2.);

                SizeType reminder = r.begin() % 2;
                SizeType r_begin  = r.begin() + reminder;

                if(reminder) {

                    if(j + 1 < n_coarse) {
                        I.set(r.begin(), j, 0.5/h);
                        I.set(r.begin(), j + 1, 0.5/h);
                    }

                    ++j;
                }

                for(auto k = r_begin; k < r.end(); k += 2, ++j) {
                    I.set(k, j, 1./h);

                    if(j + 1 < n_coarse) {
                        auto kp1 = k + 1;

                        if(r.inside(kp1)) {
                            I.set(kp1, j, 0.5/h);
                            I.set(kp1, j + 1, 0.5/h);
                        }
                    }
                }
            }

            SizeType n_finest = n_dofs.back();
            matrix = std::make_shared<Matrix>(sparse(n_finest, n_finest, 3));
            assemble_laplacian_1D(*matrix, true);

            rhs = std::make_shared<Vector>(values(n_finest, h*10.));


            Write<Vector> w_(*rhs);

            auto r = range(*rhs);
            if(r.begin() == 0) {
                rhs->set(0, -1.);
            }

            if(r.end() == n_finest) {
                rhs->set(n_finest - 1, -1.);
            }

            if(remove_bc_dofs_from_interp) {
                auto &I = *interpolators.back();

                Write<Matrix> w_(I);
                auto rr = row_range(I);

                if(rr.inside(0)) {
                    I.set(0, 0, 0.);
                }

                auto last_node_h = size(I).get(0) - 1;
                auto last_node_H = size(I).get(1) - 1;
                if(rr.inside(last_node_h)) {
                    I.set(last_node_h, last_node_H, 0.);
                }
            }
        }

        void describe() const
        {
            disp("----------------------------------");
            for(auto I_ptr : interpolators) {
                disp(*I_ptr);
            }
            disp("----------------------------------");
        }

        void write_matlab(const std::string &folder)
        {
            SizeType i = 0;
            for(auto I_ptr : interpolators) {
                //REMOVE ME
                // I_ptr->implementation().set_name("I_" + std::to_string(i));

                write("mat_I_" + std::to_string(i) + ".m", *I_ptr);
            }

            //REMOVE ME
            // rhs->implementation().set_name("r");

            //REMOVE ME
            // matrix->implementation().set_name("A");

            write("vec_r.m", *rhs);
            write("mat_A.m", *matrix);
        }

        std::vector<SizeType> n_dofs;
        std::vector<std::shared_ptr<Matrix>> interpolators;
        std::shared_ptr<Matrix> matrix;
        std::shared_ptr<Vector> rhs;
    };
}

#endif //UTOPIA_MULTILEVEL_TEST_PROBLEM_HPP
