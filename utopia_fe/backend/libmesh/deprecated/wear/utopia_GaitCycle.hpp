#ifndef UTOPIA_GAIT_CYCLE_HPP
#define UTOPIA_GAIT_CYCLE_HPP

#include <array>
#include "utopia.hpp"
#include "utopia_libmesh_Types.hpp"

#include "libmesh/dof_map.h"
#include "libmesh/mesh.h"

#include "utopia_AffineTransform.hpp"
#include "utopia_Input.hpp"
#include "utopia_ProductFunctionSpace.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

namespace utopia {

    class InputStream;

    class GaitCycle final : public Configurable {
    public:
        ////////////////////////////////////////////////////////////////////////////////////////

        class Configuration : public Configurable {
        public:
            Configuration() {}
            virtual ~Configuration() {}

            virtual void update(const int time_step) = 0;
            virtual void init(ProductFunctionSpace<LibMeshFunctionSpace> &) {}
            virtual void describe(std::ostream &) const {}

            virtual void displacement_and_forces(ProductFunctionSpace<LibMeshFunctionSpace> &space,
                                                 UVector &displacement,
                                                 UVector &forces) const = 0;

            virtual int n_steps() const = 0;
            virtual double dt() const = 0;
        };

        ////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////

        GaitCycle();
        void read(Input &is) override;
        void update(const std::size_t time_step);

        inline int n_steps() const { return conf_->n_steps(); }

        inline void init(ProductFunctionSpace<LibMeshFunctionSpace> &space) { conf_->init(space); }

        void displacement_and_forces(ProductFunctionSpace<LibMeshFunctionSpace> &space,
                                     UVector &displacement,
                                     UVector &forces) const;

        virtual void describe(std::ostream &os) const {
            if (conf_) conf_->describe(os);
        }

        Configuration &conf() {
            assert(conf_);
            return *conf_;
        }

        const Configuration &conf() const {
            assert(conf_);
            return *conf_;
        }

        std::shared_ptr<Configuration> conf_;
        std::map<std::string, std::shared_ptr<Configuration>> conf_types_;
    };
}  // namespace utopia

#endif  // UTOPIA_GAIT_CYCLE_HPP
