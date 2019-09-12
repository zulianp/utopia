#ifndef UTOPIA_BOX_CONSTRAINTS_HPP
#define UTOPIA_BOX_CONSTRAINTS_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_Factory.hpp"

#include <memory>
#include <limits>

namespace utopia {

    template<class Vector>
    class BoxConstraints {
    public:
        DEF_UTOPIA_SCALAR(Vector); 
        typedef UTOPIA_SIZE_TYPE(Vector)  SizeType;

        BoxConstraints(const std::shared_ptr<Vector> &lower_bound,
                       const std::shared_ptr<Vector> &upper_bound): 
        lower_bound_(lower_bound),
        upper_bound_(upper_bound),
        min_val_(-std::numeric_limits<Scalar>::max()),
        max_val_(std::numeric_limits<Scalar>::max()), 
        uniform_(false)
        {

        }

        BoxConstraints() {}

        inline std::shared_ptr<Vector> &upper_bound()
        {
            return upper_bound_;
        }

        inline std::shared_ptr<const Vector> upper_bound() const
        {
            return upper_bound_;
        }

        inline std::shared_ptr<Vector> &lower_bound()
        {
            return lower_bound_;
        }

        inline std::shared_ptr<const Vector> lower_bound() const
        {
            return lower_bound_;
        }

        inline bool has_lower_bound() const
        {
            return static_cast<bool>(lower_bound_);
        }

        inline bool has_upper_bound() const
        {
            return static_cast<bool>(upper_bound_);
        }

        inline bool has_bound() const
        {
            return has_lower_bound() || has_upper_bound();
        }

        inline void fill_empty_bounds()
        {
            if(lower_bound_ == nullptr && upper_bound_ == nullptr) {
                return;
            }

            if(!lower_bound_) {
                lower_bound_ = std::make_shared<Vector>(local_values(local_size(*upper_bound_), min_val_));
            }

            if(!upper_bound_) {
                upper_bound_ = std::make_shared<Vector>(local_values(local_size(*lower_bound_), max_val_));
            }
        }  

        inline bool uniform() const 
        {
            return uniform_; 
        }

        inline void uniform(const bool & flg)
        {
            uniform_ = flg; 
        }

    private:
        std::shared_ptr<Vector> lower_bound_;
        std::shared_ptr<Vector> upper_bound_;
        Scalar min_val_;
        Scalar max_val_;
        bool uniform_; 
    };


    template<class Vector>
    inline BoxConstraints<Vector> make_box_constaints(const std::shared_ptr<Vector> &lower_bound,
                                                      const std::shared_ptr<Vector> &upper_bound)
    {
        return BoxConstraints<Vector>(lower_bound, upper_bound);
    }

    template<class Vector>
    inline BoxConstraints<Vector> make_lower_bound_constraints(const std::shared_ptr<Vector> &lower_bound)
    {
        return BoxConstraints<Vector>(lower_bound, nullptr);
    }

    template<class Vector>
    inline BoxConstraints<Vector> make_upper_bound_constraints(const std::shared_ptr<Vector> &upper_bound)
    {
        return BoxConstraints<Vector>(nullptr, upper_bound);
    }



}

#endif //UTOPIA_BOX_CONSTRAINTS_HPP
