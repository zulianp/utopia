#ifndef UTOPIA_PETSC_DM_HPP
#define UTOPIA_PETSC_DM_HPP

#include <memory>
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_Box.hpp"

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

        SizeType dim() const;
        Range local_node_range() const;
        SizeType n_local_nodes_with_ghosts() const;

        void dims(SizeType *arr) const;
        void box(Scalar *min, Scalar *max) const;

        void local_node_ranges(SizeType *begin, SizeType *end) const;
        void local_element_ranges(SizeType *begin, SizeType *end) const;

        template<class Array>
        void local_element_ranges(Array &begin, Array &end) const
        {
            SizeType temp_begin[3], temp_end[3];
            local_element_ranges(temp_begin, temp_end);

            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                begin[i] = temp_begin[i];
                end[i] = temp_end[i];
            }
        }

        template<class Array>
        void local_node_ranges(Array &begin, Array &end) const
        {
            SizeType temp_begin[3], temp_end[3];
            local_node_ranges(temp_begin, temp_end);

            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                begin[i] = temp_begin[i];
                end[i] = temp_end[i];
            }
        }

        template<class Array>
        void dims(Array &arr) const
        {
            SizeType temp[3];
            dims(temp);
            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                arr[i] = temp[i];
            }
        }

        template<class Array>
        void box(Array &min, Array &max) const
        {
            Scalar temp_min[3], temp_max[3];
            box(temp_min, temp_max);

            SizeType d = dim();
            for(SizeType i = 0; i < d; ++i) {
                min[i] = temp_min[i];
                max[i] = temp_max[i];
            }
        }

    private:
        std::unique_ptr<Impl> impl_;
    };

}

#endif //UTOPIA_PETSC_DM_HPP
