#ifndef UTOPIA_ONE_LEVEL_HPP
#define UTOPIA_ONE_LEVEL_HPP

#include "utopia_Core.hpp"
#include "utopia_Path.hpp"

#include <cassert>
#include <memory>
#include <utility>

namespace utopia {
    /**
     * @brief      Class keeps track on operators.
     *
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template <class Matrix, class Vector>
    class Level {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

    public:
        Level() : _A(std::make_shared<const Matrix>()) {}
        Level(std::shared_ptr<const Matrix> A) : _A(std::move(A)) { assert(_A); }

        virtual ~Level() = default;

        /**
         * @brief      Setter of stifness matrix.
         *
         * @param[in]  A     The stifness matrix.
         *
         */
        bool A(const std::shared_ptr<const Matrix> &A) {
            assert(A);
            _A = A;
            return true;
        }

        /**
         * @brief      Getter for stifness matrix.
         *
         * @return     The stifness on given level.
         */
        const Matrix &A() const {
            assert(_A);
            return *_A;
        }

        /**
         * @brief      Getter for stifness matrix.
         *
         * @return     The stifness on given level.
         */
        std::shared_ptr<const Matrix> A_ptr() {
            assert(_A);
            return _A;
        }

        /**
         * @brief      Enforce active set to the system on given level.
         *
         * @param[in]  active_set  The active set (indices).
         * @param      x           The solution vector.
         * @param      b           The right hand side.
         */
        bool enforce_active_set(const std::vector<SizeType> &active_set, Vector &x, Vector &b) {
            assert(_A);
            return apply_BC_to_system(_A, x, b, active_set);
        }

        virtual bool write(const Path &path) const {
            if (_A) {
                utopia::write(path.c_str(), *_A);
                return true;
            }

            return false;
        }

    protected:
        std::shared_ptr<const Matrix> _A;
    };

}  // namespace utopia

#endif  // UTOPIA_ONE_LEVEL_HPP
