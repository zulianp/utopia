#ifndef UTOPIA_FE_FORM_EVAL_TEST_HPP
#define UTOPIA_FE_FORM_EVAL_TEST_HPP

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

    class FormEvalTest final : public FETest {
    public:
        void run(Input &in) override;

        inline static std::string command()
        {
            return "fet";
        }
    };

}


#endif //UTOPIA_FE_FORM_EVAL_TEST_HPP
