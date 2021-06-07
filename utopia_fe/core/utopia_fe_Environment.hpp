#ifndef UTOPIA_FE_ENVIRONMENT_HPP
#define UTOPIA_FE_ENVIRONMENT_HPP

#include "utopia_IOStream.hpp"
#include "utopia_Instance.hpp"

#include "utopia_Field.hpp"

#include <cassert>
#include <map>
#include <memory>
#include <string>

namespace utopia {

    template <class FunctionSpace>
    class Environment : public Describable {
    public:
        // using Vector_t = typename Traits<FunctionSpace>::Vector;
        // using Matrix_t = typename Traits<FunctionSpace>::Matrix;
        // using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;

        bool add_field(const std::shared_ptr<Field<FunctionSpace>> &field) {
            if (field->empty()) {
                utopia::err() << "Environment: Cannot add empty field to environment!\n";
                assert(false);
                return false;
            }

            if (field->name().empty()) {
                utopia::err() << "Environment: Cannot add field without a proper name to environment!\n";
                assert(false);
                return false;
            }

            const Scalar_t norm_field = norm2(field->data());
            utopia::out() << "Environment: read field with norm: " << norm_field << '\n';

            auto space = field->space();
            auto ret = space_to_fields_[space->name()].fields.insert(std::make_pair(field->name(), field));

            if (!ret.second) {
                utopia::err() << "Environment: Field with name " + field->name() +
                                     " already exists, name of field must be unique\n";
                Utopia::Abort();
            }

            return true;
        }

        std::shared_ptr<Field<FunctionSpace>> find_field(const FunctionSpace &space, const std::string &name) const {
            auto fields_it = space_to_fields_.find(space.name());

            if (fields_it == space_to_fields_.end()) {
                utopia::err() << "Environment: No fields for space with name " << space.name() << "!\n";
                return nullptr;
            }

            auto f = fields_it->second.fields.find(name);
            if (f == fields_it->second.fields.end()) {
                utopia::err() << "Environment: field with name " << name << " does not exists, returning null\n";
                return nullptr;
            }

            return f->second;
        }

        bool add_space(const std::shared_ptr<FunctionSpace> &space) {
            if (space->empty()) {
                utopia::err() << "Environment: Cannot add empty space to environment!\n";
                assert(false);
                return false;
            }

            if (space->name().empty()) {
                utopia::err() << "Environment: Cannot add space without a proper name to environment!\n";
                assert(false);
                return false;
            }

            auto ret = spaces_.insert(std::make_pair(space->name(), space));

            if (!ret.second) {
                utopia::err() << "Environment: Space with name " + space->name() +
                                     " already exists, name of space must be unique\n";
                Utopia::Abort();
            }

            return true;
        }

        std::shared_ptr<FunctionSpace> find_space(const std::string &name) const {
            auto f = spaces_.find(name);
            if (f == spaces_.end()) {
                utopia::err() << "Environment: space with name " << name << " does not exists, returning null\n";
                return nullptr;
            }

            return f->second;
        }

        void describe(std::ostream &os) const override {
            for (auto &stf : space_to_fields_) {
                os << "Space: " << stf.first << " fields:\n";

                for (auto &f : stf.second.fields) {
                    os << "\t- " << f.first << "\n";
                }
            }
        }

    private:
        class Fields {
        public:
            std::map<std::string, std::shared_ptr<Field<FunctionSpace>>> fields;
        };

        std::map<std::string, Fields> space_to_fields_;
        std::map<std::string, std::shared_ptr<FunctionSpace>> spaces_;
    };
}  // namespace utopia

#endif