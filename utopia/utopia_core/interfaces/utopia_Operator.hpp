/**
 * @file utopia_Operator.hpp
 * @brief Defines the utopia::Operator class and related functions.
 */

#ifndef UTOPIA_OPERATOR_HPP
#define UTOPIA_OPERATOR_HPP

#include "utopia_DistributedObject.hpp"
#include "utopia_Size.hpp"
#include "utopia_Traits.hpp"

namespace utopia {

    /**
     * @brief A pure virtual base class for defining linear operators on vectors.
     *
     * This class represents a linear operator that can be applied to vectors. It defines
     * common interface functions for applying the operator, retrieving its size, and
     * accessing its communicator.
     *
     * @tparam Vector The vector type the operator operates on.
     */
    template <class Vector>
    class Operator : public virtual DistributedObject {
    public:
        /**
         * @brief Alias for the communicator type associated with the vector.
         */
        using Communicator = typename Traits<Vector>::Communicator;

        /**
         * @brief Destructor for the Operator class.
         */
        ~Operator() override = default;

        /**
         * @brief Applies the operator to a given vector.
         *
         * This pure virtual function applies the operator to the provided 'rhs' vector
         * and stores the result in the 'sol' vector.
         *
         * @param rhs The right-hand-side vector.
         * @param sol The solution vector.
         * @return True if the operation was successful, false otherwise.
         */
        virtual bool apply(const Vector &rhs, Vector &sol) const = 0;

        /**
         * @brief Gets the size of the operator.
         *
         * This pure virtual function returns the size (dimension) of the operator.
         *
         * @return The size of the operator.
         */
        virtual Size size() const = 0;

        /**
         * @brief Gets the local size of the operator.
         *
         * This pure virtual function returns the local size (dimension) of the operator.
         *
         * @return The local size of the operator.
         */
        virtual Size local_size() const = 0;

        /**
         * @brief Gets the communicator associated with the operator.
         *
         * This pure virtual function returns a reference to the communicator associated
         * with the operator.
         *
         * @return A const reference to the communicator.
         */
        const Communicator &comm() const override = 0;
    };

    /**
     * @brief Casts an operator to a specific vector type.
     *
     * This function is used to cast an operator of a different vector type to an
     * Operator of the specified 'Vector' type.
     *
     * @tparam Vector The desired vector type.
     * @tparam T The type of the operator to be casted.
     * @param op The operator to be casted.
     * @return A reference to the casted operator of type 'Vector'.
     */
    template <class Vector, class T>
    const Operator<Vector> &operator_cast(const T &op) {
        return static_cast<const Operator<Vector> &>(op);
    }

}  // namespace utopia

#endif  // UTOPIA_OPERATOR_HPP
