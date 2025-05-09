#ifndef UTOPIA_FUNCTION_SPACE_BASE_HPP
#define UTOPIA_FUNCTION_SPACE_BASE_HPP

#include "utopia_Describable.hpp"
#include "utopia_FECoreForwardDeclarations.hpp"
#include "utopia_Field.hpp"
#include "utopia_SimulationTime.hpp"
#include "utopia_Traits.hpp"

#include <memory>
#include <string>

namespace utopia {

    template <class Mesh_>
    class Traits<FunctionSpaceBase<Mesh_>> : public Traits<Mesh_> {
    public:
        using Mesh = Mesh_;
    };

    template <class Mesh>
    class FunctionSpaceBase : public Configurable, public Describable {
    public:
        using Vector = typename Traits<FunctionSpaceBase>::Vector;
        using Matrix = typename Traits<FunctionSpaceBase>::Matrix;
        using Scalar = typename Traits<FunctionSpaceBase>::Scalar;
        using SizeType = typename Traits<FunctionSpaceBase>::SizeType;
        using IndexSet = typename Traits<FunctionSpaceBase>::IndexSet;
        using Communicator = typename Traits<FunctionSpaceBase>::Communicator;

        virtual ~FunctionSpaceBase() = default;

        virtual void init(const std::shared_ptr<Mesh> &mesh) = 0;
        virtual std::shared_ptr<Mesh> mesh_ptr() const = 0;
        virtual const Mesh &mesh() const = 0;
        virtual Mesh &mesh() = 0;
        virtual int n_var() const = 0;

        virtual void update(const SimulationTime<Scalar> &) {
            assert(false && "IMPLEMENT ME!");
            Utopia::Abort("FunctionSpaceBase::update(time) must be implemented in subclass!");
        }

        // Below could be made into a cross-backend interface (given that the same algebra is used)
        virtual bool write(const Path &path, const Vector &x) = 0;

        virtual const Communicator &comm() const = 0;

        virtual SizeType n_dofs() const = 0;
        virtual SizeType n_local_dofs() const = 0;

        virtual void create_vector(Vector &v) const = 0;
        virtual void create_matrix(Matrix &m) const = 0;

        virtual void apply_constraints(Matrix &m, const Scalar diag_value = 1.0) const = 0;
        virtual void apply_constraints(Vector &v) const = 0;
        virtual void apply_constraints_update(Vector &v) const { apply_constraints(v); }
        virtual void apply_constraints(Matrix &m, Vector &v) const = 0;
        virtual void apply_constraints_time_derivative(Vector &) const { assert(false); }
        virtual void apply_zero_constraints(Vector &vec) const = 0;

        virtual void copy_at_constrained_nodes(const Vector &, Vector &) const {
            assert(false && "IMPLEMENT ME!");
            Utopia::Abort("FunctionSpaceBase::copy_at_constrained_nodes(in, out) must be implemented in subclass!");
        }

        virtual void add_dirichlet_boundary_condition(const std::string &name,
                                                      const Scalar &value,
                                                      const int component = 0) = 0;

        virtual bool empty() const = 0;

        virtual void displace(const Vector &displacement) = 0;

        virtual const std::string &name() const = 0;
        virtual void initialize() = 0;

        virtual bool is_non_conforming() const { return false; }
        virtual std::shared_ptr<Matrix> constraint_matrix() const { return nullptr; }
    };

}  // namespace utopia

#endif  // UTOPIA_FUNCTION_SPACE_BASE_HPP