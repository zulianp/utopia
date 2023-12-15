
#include "utopia_Testing.hpp"
#include "utopia_ui.hpp"

#include "utopia.hpp"
#include "utopia_InputParameters.hpp"
#include "utopia_Instance.hpp"
#include "utopia_SymbolicFunction.hpp"
#include "utopia_ui.hpp"

#ifdef UTOPIA_WITH_YAML_CPP
#include "utopia_YAMLInput.hpp"
#endif

namespace utopia {

    void generic_stream(Input &is) {
        utopia_test_assert(is.good());
        // TODO(zulianp):
    }

    void xml_stream() {
        Path path = Utopia::instance().get("data_path") + "/xmlsamples/xml_test.xml";

        auto is_ptr = open_istream(path);
        utopia_test_assert(bool(is_ptr));

        if (!is_ptr) {
            return;
        }
        generic_stream(*is_ptr);
    }

#ifdef UTOPIA_WITH_TINY_EXPR
    void symbolic_expr() {
        {
            SymbolicFunction f("x + y + z");
            double w = f.eval(1, 2, 3);
            utopia_test_assert(f.valid());
            utopia_test_assert(approxeq(w, 6.));

            w = f.eval(1, 2);
            utopia_test_assert(approxeq(w, 3.));
        }

        {
            SymbolicFunction f("x*y");
            double w = f.eval({2, 2, 3});
            utopia_test_assert(f.valid());
            utopia_test_assert(approxeq(w, 4.));
        }
    }

#endif  // UTOPIA_WITH_TINY_EXPR

#ifdef UTOPIA_WITH_YAML_CPP

    void yaml_vector_test() {
        std::string yaml = "vec: [1, 2, 3]";
        auto in = YAMLInput::FromString(yaml);

        {
            std::vector<std::string> v;
            in->get("vec", v);
            utopia_test_assert(v.size() == 3);
        }

        {
            std::vector<double> v;
            in->get("vec", v);
            utopia_test_assert(v.size() == 3);
        }

        {
            std::vector<float> v;
            in->get("vec", v);
            utopia_test_assert(v.size() == 3);
        }
    }
#endif

    void input_parameters() {
        InputParameters in;
        in.set("string-key", std::string("value"));
        in.set("double-key", 1.);
        in.set("int-key", 20);

        double d_value = 2.;
        int i_value = 10;
        std::string s_value = "blah";

        in.get("string-key", s_value);
        in.get("double-key", d_value);
        in.get("int-key", i_value);

        utopia_test_assert(s_value == "value");
        utopia_test_assert(d_value == 1.);
        utopia_test_assert(i_value == 20);

        int inexistent_key = -6;
        in.get("inexistent-key", inexistent_key);
        utopia_test_assert(inexistent_key == -6);
    }

    void newton_ui() {
        const std::string data_path = Utopia::instance().get("data_path");

#ifdef UTOPIA_WITH_PETSC

        // auto cg = std::make_shared<ConjugateGradient<PetscMatrix, PetscVector, HOMEMADE>>();
        auto cg = std::make_shared<ConjugateGradient<PetscMatrix, PetscVector>>();
        Newton<PetscMatrix, PetscVector> newton(cg);

#ifdef UTOPIA_WITH_JSON
        newton.import("Newton", data_path + "/json/default.json");
#endif  // UTOPIA_WITH_JSON

        newton.import("Newton", data_path + "/xml/default.xml");
#endif  // UTOPIA_WITH_PETSC
    }

    static void ui() {
        UTOPIA_RUN_TEST(xml_stream);
        UTOPIA_RUN_TEST(input_parameters);
        UTOPIA_RUN_TEST(newton_ui);
#ifdef UTOPIA_WITH_TINY_EXPR
        UTOPIA_RUN_TEST(symbolic_expr);
#endif  // UTOPIA_WITH_TINY_EXPR

#ifdef UTOPIA_WITH_YAML_CPP
        UTOPIA_RUN_TEST(yaml_vector_test);
#endif
    }

    UTOPIA_REGISTER_TEST_FUNCTION(ui);

}  // namespace utopia
