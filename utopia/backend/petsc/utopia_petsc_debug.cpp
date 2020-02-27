#include "utopia_petsc_debug.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_petsc_Communicator.hpp"

#include "petscsys.h"
#include <stack>
#include <iostream>
#include <cassert>

namespace utopia {

    class PetscDebugger::Impl {
    public:
        //number of bytes currently allocated
        std::stack<PetscLogDouble> alloc_stack;
    };

    PetscDebugger::PetscDebugger()
    : impl_(utopia::make_unique<Impl>())
    {}

    PetscDebugger::~PetscDebugger()
    {

    }

    PetscDebugger &PetscDebugger::instance()
    {
        static PetscDebugger instance_;
        return instance_;
    }

    void PetscDebugger::memory_check_begin()
    {
        PetscLogDouble space;
        PetscErrorCode ierr = PetscMallocGetCurrentUsage(&space); assert(ierr == 0); (void) ierr;
        impl_->alloc_stack.push(space);
    }

    void PetscDebugger::memory_check_end()
    {
        PetscErrorCode ierr;
        PetscLogDouble space;
        ierr = PetscMallocGetCurrentUsage(&space); assert(ierr == 0); (void) ierr;

        PetscLogDouble previous_space = impl_->alloc_stack.top();

        if(space > previous_space) {
            std::cerr << "[Memory check] expected "
                      << (previous_space/1024) << "KB have "
                      << (space/1024) << "KB diff "
                      << ((space - previous_space)/1024) << "KB"
                      << std::endl;

            // ierr = PetscMallocDump(PETSC_STDOUT); assert(ierr == 0);
            // assert(space <= previous_space);
        }

        impl_->alloc_stack.pop();
    }

    void PetscDebugger::print_current_usage(std::ostream &os) const
    {
        PetscErrorCode ierr;
        PetscLogDouble space;
        ierr = PetscMemoryGetCurrentUsage(&space); assert(ierr == 0); (void) ierr;
        os << "[Memory Usage]  " << (space/1024) << "KB" << "\n";
    }

    void PetscDebugger::print_current_collective_usage(const std::string &marker, std::ostream &os) const
    {
        PetscErrorCode ierr;
        PetscLogDouble space;
        ierr = PetscMemoryGetCurrentUsage(&space); assert(ierr == 0); (void) ierr;

        
        space /= 1024; //to KB
        space /= 1024; //to MB
        space /= 1000; //to GB

        PetscCommunicator comm(PETSC_COMM_WORLD);
        space = comm.sum(space);
        comm.root_print("[Memory Usage]  (" + marker + ") " + std::to_string(space) + "GB");
    }

    void PetscDebugger::describe(std::ostream &os) const
    {
        //TODO
    }

}