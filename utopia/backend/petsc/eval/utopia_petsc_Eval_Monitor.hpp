#ifndef UTOPIA_PETSC_EVAL_MONITOR_HPP
#define UTOPIA_PETSC_EVAL_MONITOR_HPP

#include "utopia_Monitoring.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_Traits.hpp"

namespace utopia {

    template<>
    class EvalMonitor<PetscMatrix, PETSC> {
    public:
        using SizeType = Traits<PetscMatrix>::SizeType;

        inline static void apply(
            const SizeType &it,
            const PetscMatrix &t
        )
        {   
            apply(it, t, "H_log.m", "H");
        }

        static void apply(
            const SizeType &it,
            const PetscMatrix &t,
            const std::string &name_of_file,
            const std::string &name_of_instance
        );
    };


    template<>
    class EvalMonitor<PetscVector, PETSC> {
    public:
        using SizeType = Traits<PetscVector>::SizeType;

        inline static void apply(
            const SizeType &it,
            const PetscVector &t
        )
        {   
            apply(it, t, "g_log.m", "g");
        }

        static void apply(
            const SizeType &it,
            const PetscVector &t,
            const std::string &name_of_file,
            const std::string &name_of_instance
        );
    };

}

#endif //UTOPIA_PETSC_EVAL_MONITOR_HPP