#ifndef UTOPIA_OPTIONS_HPP
#define UTOPIA_OPTIONS_HPP

#include "utopia_Describable.hpp"
#include "utopia_Input.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_MPI.hpp"
#include "utopia_TypeToString.hpp"

#include <iostream>

namespace utopia {

    class Options final : public Configurable, public Describable {
    public:
        class OptBase : public Configurable, public Describable {
        public:
            virtual ~OptBase() = default;
        };

        template <typename T>
        class Opt : public OptBase {
        public:
            void read(Input &in) override { in.get(key, ref); }

            void describe(std::ostream &os) const override {
                os << '-' << key << " <" << TypeToString<T>::get() << ">\tdefault: " << ref << "\t\t";
                os << explanation << "\n";
            }

            std::string key;
            T &ref;
            std::string explanation;

            Opt(std::string key, T &ref, std::string explanation)
                : key(std::move(key)), ref(ref), explanation(std::move(explanation)) {}
        };

        template <typename T>
        Options &add_option(std::string key, T &ref, std::string explanation) {
            args_.push_back(utopia::make_unique<Opt<T>>(key, ref, explanation));
            return *this;
        }

        bool parse(Input &in, std::ostream &os = std::cout) {
            bool help = false;
            in.get("help", help);

            if (!help) {
                in.get("show", help);
            }

            if (help) {
                describe(os);
                return false;
            } else {
                read(in);
            }

            return true;
        }

        void describe(std::ostream &os) const override {
            if (mpi_world_rank() == 0) {
                os << "Options are:\n";
                for (auto &a : args_) {
                    a->describe(os);
                }
            }
        }

        void read(Input &in) override {
            for (auto &a : args_) {
                a->read(in);
            }
        }

    private:
        std::vector<std::unique_ptr<OptBase>> args_;
    };
}  // namespace utopia

#endif  // UTOPIA_OPTIONS_HPP
