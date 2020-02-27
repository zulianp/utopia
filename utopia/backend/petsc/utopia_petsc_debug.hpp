#ifndef UTOPIA_PETSC_DEBUG_HPP
#define UTOPIA_PETSC_DEBUG_HPP

#include <memory>
#include <iostream>

namespace utopia {

    class PetscDebugger {
    public:
        ~PetscDebugger();
        static PetscDebugger &instance();

        void memory_check_begin();
        void memory_check_end();
        void print_current_usage(std::ostream &os = std::cout) const;
        void print_current_collective_usage(const std::string &marker, std::ostream &os = std::cout) const;
        void describe(std::ostream &os) const;

    private:
        class Impl;
        PetscDebugger();
        std::unique_ptr<Impl> impl_;
    };
}

#define UTOPIA_PETSC_COLLECTIVE_MEMUSAGE() {  PetscDebugger::instance().print_current_collective_usage(__FILE__ + std::string(":") + std::to_string(__LINE__)); }
#define UTOPIA_PETSC_COLLECTIVE_MEMUSAGE(macro_marker_) {  PetscDebugger::instance().print_current_collective_usage(macro_marker_); }
#define UTOPIA_PETSC_MEMUSAGE() {  PetscDebugger::instance().print_current_usage(); }
#define UTOPIA_PETSC_MEMCHECK_BEGIN() { PetscDebugger::instance().memory_check_begin(); }
#define UTOPIA_PETSC_MEMCHECK_END() { PetscDebugger::instance().memory_check_end(); }

#endif //UTOPIA_PETSC_DEBUG_HPP
