#ifndef UTOPIA_FACTORY_METHOD_HPP
#define UTOPIA_FACTORY_METHOD_HPP

#include <memory>
#include "utopia_make_unique.hpp"

namespace utopia {

    template<class OutInterface>
    class IFactoryMethod {
    public:
        virtual ~IFactoryMethod()
        {}

        virtual std::unique_ptr<OutInterface> make() const = 0;
    };


    template<class OutInterface, class OutObject = OutInterface>
    class FactoryMethod final : public IFactoryMethod<OutInterface> {
    public:
        FactoryMethod()
        {}

        std::unique_ptr<OutInterface> make() const final override
        {
            return utopia::make_unique<OutObject>();
        }
    };

    template<class OutInterface, class ArgIn, class OutObject = OutInterface>
    class UnaryFactoryMethod final : public IFactoryMethod<OutInterface> {
    public:
        UnaryFactoryMethod(const ArgIn &arg)
        : arg_(arg)
        {}

        std::unique_ptr<OutInterface> make() const final override
        {
            return utopia::make_unique<OutObject>(arg_);
        }

    private:
        ArgIn arg_;
    };


    template<class OutInterface, class OutObject>
    inline std::unique_ptr<IFactoryMethod<OutInterface>> make_factory()
    {
        return utopia::make_unique<FactoryMethod<OutInterface, OutObject>>();
    }

    template<class OutInterface, class OutObject, class Args>
    inline std::unique_ptr<IFactoryMethod<OutInterface>> make_factory(const Args &args)
    {
        return utopia::make_unique<UnaryFactoryMethod<OutInterface, Args, OutObject>>(args);
    }
}

#endif
