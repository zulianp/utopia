#ifndef UTOPIA_PETSC_DM_HPP
#define UTOPIA_PETSC_DM_HPP

#include <memory>
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"

namespace utopia {

    class PetscDM final {
    public:
        using SizeType = PetscInt;
        using Scalar   = PetscScalar;

        class Impl;
        class Elements;
        class Nodes;

        void create_matrix(PetscMatrix &mat);

        PetscDM(
            const PetscCommunicator     &comm,
            const std::vector<SizeType> &dims,
            const std::vector<Scalar>   &box_min,
            const std::vector<Scalar>   &box_max
        );

        PetscDM();
        ~PetscDM();

        void describe() const;

        class Elem {
        public:

            Elem(const Elements &e, const SizeType &idx) : e_(e), idx_(idx) {}
            SizeType n_nodes() const;
            SizeType node_id(const SizeType k) const;
            SizeType idx() const { return idx_; }

        private:
            const Elements &e_;
            SizeType idx_;
        };

        class Node {
        public:
            Node(const Nodes &nodes, const SizeType &idx) : nodes_(nodes), idx_(idx) {}
            SizeType idx() const { return idx_; }
            bool is_ghost() const;

            private:
                const Nodes &nodes_;
                SizeType idx_;
        };

        void each_element(const std::function<void(const Elem &)> &f);
        void each_node(const std::function<void(const Node &)> &f);
        void each_node_with_ghosts(const std::function<void(const Node &)> &f);

    private:
        std::unique_ptr<Impl> impl_;
    };

}

#endif //UTOPIA_PETSC_DM_HPP
