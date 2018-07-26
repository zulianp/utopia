
#include "utopia_UITest.hpp"
#include "utopia_ui.hpp"

#include "utopia.hpp"
#include "utopia_ui.hpp"
#include "utopia_Instance.hpp"
#include "utopia_SymbolicFunction.hpp"

namespace utopia {
	void generic_stream(InputStream &is)
	{
		utopia_test_assert(is.good());

		//set-default values here
		std::string type = "undefined", op = "undefined", solver = "undefined", algorithm = "undefined";
		SizeType max_iter = -1;

		//this should not be overriden
		double tol = 1e-16;

		utopia_test_assert( is.object_begin("solve") );

		{
			is.read("type", type);
			is.read("operator", op);

			utopia_test_assert( is.object_begin("solver") );

			{
				is.read("algorithm", algorithm);
				is.read("max_iter", max_iter);
				is.read("tol", tol);
			}

			utopia_test_assert( is.object_end() ); //close solver

		}

		utopia_test_assert( is.object_end() ); //close solve

		utopia_test_assert( type == "linear" );
		utopia_test_assert( op == "../data/mg/A.bin" );
		utopia_test_assert( algorithm == "CG" );
		utopia_test_assert( max_iter == 10 );
		utopia_test_assert( tol == 1e-16 );
	}

	void xml_stream()
	{
		/** XML file
		<solve>
			<type>linear</type>
			<operator>../data/mg/A.bin</operator>
			<rhs>../data/mg/rhs.bin</rhs>
			<solver>
				<algorithm>CG</algorithm>
				<max_iter>10</max_iter>
			</solver>
		</solve>
		*/

		Path path = Utopia::instance().get("data_path") + "/xmlsamples/xml_test.xml";

		auto is_ptr = open_istream(path);
		utopia_test_assert(is_ptr);

		if(!is_ptr) return;
		generic_stream(*is_ptr);
	}

#ifdef WITH_TINY_EXPR
	void symbolic_expr()
	{
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

#endif //WITH_TINY_EXPR

	void run_ui_test()
	{
		UTOPIA_UNIT_TEST_BEGIN("UITest");
		UTOPIA_RUN_TEST(xml_stream);
#ifdef WITH_TINY_EXPR
		UTOPIA_RUN_TEST(symbolic_expr);
#endif //WITH_TINY_EXPR
		UTOPIA_UNIT_TEST_END("UITest");
	}

}
