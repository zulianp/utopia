#ifndef UTOPIA_MARS_FUNCTION_SPACE_HPP
#define UTOPIA_MARS_FUNCTION_SPACE_HPP

// #include "utopia_FunctionSpace.hpp"
#include "utopia_Traits.hpp"

#include "mars_base.hpp"
#include "mars_mesh.hpp"

#include <algorithm>
#include <cassert>
#include <memory>
#include <string>

namespace utopia {
    template <mars::Integer Dim, mars::Integer ManifoldDim>
    class MarsFunctionSpace;

    template <mars::Integer Dim, mars::Integer ManifoldDim>
    class Traits<MarsFunctionSpace<Dim, ManifoldDim>> {
    public:
        using SizeType = mars::Integer;
    };

    // P0/P1 only
    template <mars::Integer Dim, mars::Integer ManifoldDim>
    class MarsFunctionSpace {
    public:
        using Integer = mars::Integer;

        class Variable {
        public:
            Variable(const std::string &name = "u", const Integer &order = 1) : name_(name), order_(order) {}

            const std::string &name() const { return name_; }

            inline Integer order() const { return order_; }

        private:
            std::string name_;
            Integer order_;
        };

        MarsFunctionSpace(const std::shared_ptr<mars::Mesh<Dim, ManifoldDim>> &mesh, const Variable &var = Variable())
            : mesh_(mesh) {}

        void dofs(const Integer elem_id, std::vector<Integer> &indices) {
            const auto &e = mesh_->elem(elem_id);

            if (var_.order() == 1) {
                // P1
                indices.resize(e.nodes.size());
                std::copy(std::begin(e.nodes), std::end(e.nodes), std::begin(indices));
            } else {
                // P0
                assert(var_.order() == 0);
                indices.resize(1);
                indices[0] = elem_id;
            }
        }

        inline const Variable &var() const { return var_; }

    private:
        Variable var_;
        std::shared_ptr<mars::Mesh<Dim, ManifoldDim>> mesh_;
    };

}  // namespace utopia

#endif  // UTOPIA_MARS_FUNCTION_SPACE_HPP
