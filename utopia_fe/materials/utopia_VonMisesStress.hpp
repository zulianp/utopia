#ifndef UTOPIA_VON_MISES_STRESS_HPP
#define UTOPIA_VON_MISES_STRESS_HPP

#include "utopia_libmesh_Types.hpp"
#include <cassert>

namespace utopia {
    class VonMisesStress {
    public:

        static double apply(const LMDenseMatrix &stress)
        {
            auto s = size(stress);

            if(s.get(0) == 2) {
                return apply_2(stress);
            } else {
                assert(s.get(0) == 3);
                return apply_3(stress);
            }
        }

        static double apply_2(const LMDenseMatrix &stress)
        {
            using std::sqrt;
            Read<LMDenseMatrix> r(stress);

            double result =  0.5 * ( stress.get(0,0) - stress.get(1, 1) ) *
            ( stress.get(0,0) - stress.get(1, 1) ) +
            3.0  *  stress.get(0, 1) * stress.get(0, 1);
            
            result = sqrt( fabs(result) );
            assert(result == result && "apply_2: result is nan");
            return result;
        }

        static double apply_3(const LMDenseMatrix &stress)
        {
            using std::sqrt;
            Read<LMDenseMatrix> r(stress);

            double result = 0.5 * ( stress.get(0, 0) - stress.get(0, 1) ) *
            ( stress.get(0, 0) - stress.get(0, 1) ) +
            3.0  *  stress.get(0, 1) * stress.get(0, 1);
            
            result += 0.5 * (stress.get(2, 2) - stress.get(0, 1)) * (stress.get(2, 2) - stress.get(0, 1)) + 3.0  * stress.get(2, 1) * stress.get(2, 1);
            result += 0.5 * (stress.get(2, 2) - stress.get(0, 0)) * (stress.get(2, 2) - stress.get(0, 0)) + 3.0  * stress.get(2, 0) * stress.get(2, 0);
            
            result = sqrt( fabs(result) );
            
            assert(result == result && "apply_3: result is nan");
            return result;
        }
    };
}

#endif //UTOPIA_VON_MISES_STRESS_HPP
