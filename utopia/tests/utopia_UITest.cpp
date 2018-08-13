
#include "utopia_UITest.hpp"
#include "utopia_ui.hpp"

#include "utopia.hpp"
#include "utopia_ui.hpp"
#include "utopia_Instance.hpp"
#include "utopia_SymbolicFunction.hpp"

namespace utopia {

	class SolveDesc : public Serializable {
	public:

		void read(InputStream &is) override {
			is.read("type", type);
			is.read("operator", op);
			is.read("rhs", rhs);

			//instead of creating another serializable use lambdas
			is.read("solver", [this](InputStream &sub_is) {
				sub_is.read("algorithm", algorithm);
				sub_is.read("max_iter", max_iter);
			});



			is.read("array", [this](InputStream &sub_is) {
				sub_is.read_all([this](InputStream &sub_is) {
					std::string v;
					sub_is.read(v);
					array.push_back(v);
				});
			});

			is.read("tol", tol);
			is.read("values", values);
		}

		std::string type, op, rhs;
		std::string algorithm;
		int max_iter;
		double tol = 1e-16;
		std::vector<double> values;
		std::vector<std::string> array;
	};

	void generic_stream(InputStream &is)
	{
		utopia_test_assert(is.good());

		SolveDesc desc;
		is.read("solve", desc);

		utopia_test_assert( desc.type == "linear" );
		utopia_test_assert( desc.op == "../data/mg/A.bin" );
		utopia_test_assert( desc.algorithm == "CG" );
		utopia_test_assert( desc.max_iter == 10 );
		utopia_test_assert( desc.tol == 1e-16 );

		utopia_test_assert( desc.values.size() == 3 );
		utopia_test_assert( desc.values[0] == 1. );
		utopia_test_assert( desc.values[1] == 2. );
		utopia_test_assert( desc.values[2] == 3. );
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
				<values>
					<value>1</value>
					<value>2</value>
					<value>3</value>
				</values>
				<array>
					<entry>first</entry>
					<entry>last</entry>
				</array>
			</solve>
		*/

		Path path = Utopia::instance().get("data_path") + "/xmlsamples/xml_test.xml";

		auto is_ptr = open_istream(path);
		utopia_test_assert(bool(is_ptr));

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
