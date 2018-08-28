#ifndef UTOPIA_BRATU_MULTILEVEL_TEST_PROBLEM_HPP
#define UTOPIA_BRATU_MULTILEVEL_TEST_PROBLEM_HPP

#include "utopia_Traits.hpp"
#include <memory>


namespace utopia {

		template<class Matrix, class Vector>
		class BratuMultilevelTestProblem {
		public:

			typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
			typedef UTOPIA_SCALAR(Vector) Scalar;

			BratuMultilevelTestProblem(const SizeType &n_levels = 2, bool remove_BC_contributions = false, bool verbose = false):
						n_levels(n_levels),
						n_coarse(10),
						remove_BC_contributions(remove_BC_contributions),
						verbose(verbose)
			{

				assert(n_coarse > 0);
				assert(n_levels > 1);

				n_dofs.resize(n_levels);
				n_dofs[0] = n_coarse;

				for(SizeType i = 1; i < n_levels; ++i)
					n_dofs[i] = (n_dofs[i-1] - 1) * 2 + 1;


				prolongations.resize(n_levels - 1);
				restrictions.resize(n_levels - 1);

				for(SizeType i = 0; i < n_levels - 1; ++i) {
					const auto n_coarse = n_dofs[i];
					const auto n_fine   = n_dofs[i + 1];
					prolongations[i] = std::make_shared<Matrix>(sparse(n_fine, n_coarse, 2));
					auto &I = *prolongations[i];

					Write<Matrix> w_(I);
					auto r = row_range(I);

					SizeType j = r.begin()/2;
					for(auto k = r.begin(); k < r.end(); k += 2, ++j) {
						I.set(k, j, 1.);

						if(j + 1 < n_coarse) {
							I.set(k + 1, j, 0.5);
							I.set(k + 1, j + 1, 0.5);
						}
					}
				}

				if(remove_BC_contributions)
				{
					auto &I = *prolongations.back();

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

				// restrictions, but let's use them as projections...
				// not very nice solution, but I am lazy to do something more sophisticated just for testing purposes...
				for(SizeType i = 0; i < prolongations.size(); ++i)
				{
					auto &I = *prolongations[i];
					Matrix R =  0.5 * transpose(I);
					restrictions[i] = std::make_shared<Matrix>(R);
				}
			}

			SizeType n_levels;
			SizeType n_coarse;
			bool remove_BC_contributions;

			std::vector<SizeType> n_dofs;

			std::vector<std::shared_ptr<Matrix>> prolongations;
			std::vector<std::shared_ptr<Matrix>> restrictions;

			bool verbose;
		};
}

#endif //UTOPIA_BRATU_MULTILEVEL_TEST_PROBLEM_HPP

