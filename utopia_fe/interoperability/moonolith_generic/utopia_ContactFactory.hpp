#ifndef UTOPIA_OBSTACLE_FACTORY_HPP
#define UTOPIA_OBSTACLE_FACTORY_HPP

#include "utopia_fe_Core.hpp"

#include "utopia_FECoreForwardDeclarations.hpp"

#include "utopia_ContactInterface.hpp"

#include "utopia_AnalyticObstacle_impl.hpp"
#include "utopia_ImplicitObstacle_impl.hpp"

#include "utopia_FactoryMethod.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ContactFactory {
    public:
        using ContactInterface = utopia::ContactInterface<FunctionSpace>;
        using ContactInterfacePtr = std::unique_ptr<ContactInterface>;

        // Use specialized components for function space
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;
        using Obstacle_t = utopia::Obstacle<FunctionSpace>;
        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using ImplicitObstacle_t = utopia::ImplicitObstacle<FunctionSpace>;
        using AnalyticObstacle_t = utopia::AnalyticObstacle<FunctionSpace>;
        using Contact_t = utopia::Contact<FunctionSpace>;

        static ContactInterfacePtr new_contact(Input &in) {
            auto contact = utopia::make_unique<Contact_t>();
            contact->read(in);
            return std::move(contact);
        }

        inline static ContactFactory &instance() {
            static ContactFactory instance_;
            return instance_;
        }

        template <class Type>
        static int register_contact_with_name(const std::string &name) {
            instance().factory_[name] = utopia::make_unique<FactoryMethod<ContactInterface, Type>>();
            return 0;
        }

        static ContactInterfacePtr new_contact_from_name(const std::string &name) {
            auto it = instance().factory_.find(name);

            if (it == instance().factory_.end()) {
                return nullptr;
            } else {
                return it->second->make();
            }
        }

        static ContactInterfacePtr new_obstacle(Input &in) {
            std::string type;
            in.get("type", type);

            ContactInterfacePtr obstacle;

            auto contact_ptr = instance().new_contact_from_name(type);

            if (contact_ptr) {
                contact_ptr->read(in);
                return contact_ptr;
            }

            if (type == "implicit") {
#ifdef UTOPIA_WITH_LIBMESH
                Utopia::Abort("ImplicitObstacle not supported for this backend!");
#else
                obstacle = utopia::make_unique<ImplicitObstacle_t>();
                obstacle->read(in);
#endif
            } else if (type == "analytic") {
#ifdef UTOPIA_WITH_LIBMESH
                Utopia::Abort("AnalyticObstacle not supported for this backend!");
#else
                obstacle = utopia::make_unique<AnalyticObstacle_t>();
                obstacle->read(in);
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

                obstacle = std::move(obs);
            }

            return obstacle;
        }

    private:
        std::map<std::string, std::unique_ptr<IFactoryMethod<ContactInterface>>> factory_;
    };

#define UTOPIA_CONTACT_DEFINE_VAR(macro_in) dummy_contact_variable_##macro_in##__LINE__

#define UTOPIA_CONTACT_REGISTER(FunctionSpace_, type_, name_) \
    static char UTOPIA_CONTACT_DEFINE_VAR(type_) =            \
        utopia::ContactFactory<FunctionSpace_>::register_contact_with_name<type_>(name_)

}  // namespace utopia

#endif  // UTOPIA_OBSTACLE_FACTORY_HPP
