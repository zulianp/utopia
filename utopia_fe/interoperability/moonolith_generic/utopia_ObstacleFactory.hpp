#ifndef UTOPIA_OBSTACLE_FACTORY_HPP
#define UTOPIA_OBSTACLE_FACTORY_HPP

#include "utopia_fe_Core.hpp"

#include "utopia_IObstacle.hpp"

#include "utopia_AnalyticObstacle_impl.hpp"
#include "utopia_ImplicitObstacle_impl.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ObstacleFactory {
    public:
        using IObstacle = utopia::IObstacle<FunctionSpace>;
        using IObstaclePtr = std::unique_ptr<IObstacle>;

        // Use specialized components for function space
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;
        using Obstacle_t = utopia::Obstacle<FunctionSpace>;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using ImplicitObstacle_t = utopia::ImplicitObstacle<FunctionSpace>;
        using AnalyticObstacle_t = utopia::AnalyticObstacle<FunctionSpace>;

        static IObstaclePtr new_obstacle(Input &in) {
            std::string type;
            in.get("type", type);

            IObstaclePtr obstacle;
            if (type == "implicit") {
#ifdef UTOPIA_WITH_LIBMESH
                // FIXME Once we go to c++17
                // if constexpr (IsNotSupported<ImplicitObstacle_t>::value) {
                Utopia::Abort("ImplicitObstacle not supported for this backend!");
#else
                // } else {
                obstacle = utopia::make_unique<ImplicitObstacle_t>();
                obstacle->read(in);
#endif
                // }
            } else if (type == "analytic") {
                // FIXME Once we go to c++17
                // if constexpr (IsNotSupported<AnalyticObstacle_t>::value) {
#ifdef UTOPIA_WITH_LIBMESH
                Utopia::Abort("AnalyticObstacle not supported for this backend!");
                // } else {
#else
                obstacle = utopia::make_unique<AnalyticObstacle_t>();
                obstacle->read(in);
                // }
#endif

            } else {
                auto obs = utopia::make_unique<Obstacle_t>();
                typename Obstacle_t::Params params;
                params.read(in);
                // Must be created for every process independently and the same
                Mesh_t obstacle_mesh(Communicator_t::self());
                obstacle_mesh.read(in);

                obs->set_params(params);
                obs->init_obstacle(obstacle_mesh);

                // bool export_obstacle = false;
                // in.get("export_obstacle", export_obstacle);

                // if (export_obstacle) {
                //     if (this->space()->comm().rank() == 0) {
                //         obstacle_mesh.write("obstacle.e");
                //     }
                // }

                obstacle = std::move(obs);
            }

            return obstacle;
        }
    };

}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_FACTORY_HPP
