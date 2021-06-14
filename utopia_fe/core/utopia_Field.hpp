#ifndef UTOPIA_FIELD_HPP
#define UTOPIA_FIELD_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_Instance.hpp"

#include <memory>
#include <string>

namespace utopia {

    template <class FunctionSpace>
    class Field : public Configurable {
    public:
        using Vector = typename Traits<FunctionSpace>::Vector;

        Field(const std::string &name = "",
              const std::shared_ptr<FunctionSpace> &space = nullptr,
              const std::shared_ptr<Vector> &data = nullptr)
            : name_(name), space_(space), data_(data) {}

        void read(Input &in) override {
            in.get("name", name_);
            in.get("tensor_size", tensor_size_);
            in.get("offset", offset_);
        }

        inline const std::shared_ptr<Vector> data_ptr() const { return data_; }

        void zero() {
            if (!space_) {
                Utopia::Abort("Cannot generate zero field from nullptr space");
            }

            if (!data_) {
                data_ = std::make_shared<Vector>();
                space_->create_vector(*data_);
            }

            data_->set(0.0);
        }

        inline Vector &data() {
            assert(data_);
            return *data_;
        }

        inline const Vector &data() const {
            assert(data_);
            return *data_;
        }

        inline void set_data(const std::shared_ptr<Vector> &data) { data_ = data; }
        inline void set_name(const std::string &name) { name_ = name; }

        inline const std::string &name() const { return name_; }

        inline int tensor_size() const { return tensor_size_; }
        inline void set_tensor_size(const int tensor_size) { tensor_size_ = tensor_size; }
        inline void set_offset(const int offset) { offset_ = offset; }
        inline int offset() const { return offset_; }

        /// @return true if anything goes with this field
        inline bool is_blank() const { return name_.empty(); }
        inline bool empty() const { return !data_; }

        void set_space(const std::shared_ptr<FunctionSpace> &space) { space_ = space; }
        const std::shared_ptr<FunctionSpace> &space() const { return space_; }

    private:
        std::string name_;
        std::shared_ptr<FunctionSpace> space_;
        std::shared_ptr<Vector> data_;
        int offset_{0};
        int tensor_size_{1};
    };

}  // namespace utopia

#endif  // UTOPIA_FIELD_HPP