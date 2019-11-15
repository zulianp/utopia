#include "utopia_FETests.hpp"


#include "utopia_AssemblyTest.hpp"
#include "utopia_BoundaryIntegralTest.hpp"
#include "utopia_BoundaryMeshTest.hpp"
#include "utopia_CoarsenerTest.hpp"
#include "utopia_EnergyAssemblyTest.hpp"
#include "utopia_FEEvalTest.hpp"
#include "utopia_FSITest.hpp"
#include "utopia_FormEvalTest.hpp"
#include "utopia_GeometryTest.hpp"
#include "utopia_IntersectTest.hpp"
#include "utopia_Intrepid2Test.hpp"
#include "utopia_MSHReaderTest.hpp"
#include "utopia_MechTest.hpp"
#include "utopia_NonLinearElasticityTest.hpp"
#include "utopia_SDCTest.hpp"
#include "utopia_SemigeometricMultigridTest.hpp"
#include "utopia_RefactoredContactTest.hpp"
#include "utopia_DualBasisTest.hpp"
#include "utopia_IntegratorTest.hpp"
#include "utopia_FETensorTest.hpp"
#include "utopia_NewNeohookeanTest.hpp"

namespace utopia {

    void FETests::print_usage(std::ostream &os) const
    {
        os << "----------------- test suite -------------------------\n";
        os << "-test <test_name>\n";
        os << "available tests:\n";

        for(const auto &t : tests_) {
            os << "\t" << t.first << "\n";
        }

        os << std::endl;

        os << "-test all\t for running all tests\n";
        os << "------------------------------------------------------\n";
    }

    void FETests::run(libMesh::Parallel::Communicator &comm, int argc, char * argv[])
    {
        bool all = false;

        for(int i = 1; i < argc; ++i) {
            const int ip1 = i+1;
            const int ip2 = i+2;

            if(argv[i] == std::string("-test")) {

                if(ip1 < argc) {

                    if(argv[ip1] == std::string("all")) {
                        all = true;
                        i = ip1;
                        continue;
                    }

                    auto it = tests_.find(argv[ip1]);
                    if(it == tests_.end()) {
                        std::cerr << "[Error] " << argv[ip1] << " not found" << std::endl;
                        print_usage(std::cout);
                    } else {
                        std::cout << "--------------------------------------------" << std::endl;
                        std::cout << "[Status] Running: " << argv[ip1] << std::endl;
                        std::cout << "--------------------------------------------" << std::endl;
                        it->second->init(comm);

                        if(ip2 < argc) {
                            auto in = open_istream(argv[ip2]);

                            if(in) {
                                it->second->run(*in);
                            } else {
                                InputParameters in;
                                it->second->run(in);
                            }

                        } else {
                            InputParameters in;
                            it->second->run(in);
                        }

                        std::cout << "--------------------------------------------" << std::endl;
                        std::cout << "--------------------------------------------" << std::endl;
                        std::cout << "Exiting runtime..." << std::endl;
                    }
                } else {
                    std::cerr << "[Error] run requires an input string" << std::endl;
                }
            }
        }

        if(all) {
            for(const auto &t : tests_) {
                std::cout << "running " << t.first << std::endl;
                t.second->init(comm);
                InputParameters in;
                t.second->run(in);
            }
        }
    }

    int FETests::add_test(const std::string &command, std::unique_ptr<FETest> &&app)
    {
        tests_[command] = std::move(app);
        return 0;
    }

    FETests::FETests()
    {
        add_test(AssemblyTest::command(), utopia::make_unique<AssemblyTest>());
        add_test(BoundaryIntegralTest::command(), utopia::make_unique<BoundaryIntegralTest>());
        add_test(BoundaryMeshTest::command(), utopia::make_unique<BoundaryMeshTest>());
        // add_test(CoarsenerTest::command(), utopia::make_unique<CoarsenerTest>()); //FIXME missing mesh
        add_test(EnergyTest::command(), utopia::make_unique<EnergyTest>());
        add_test(FEEvalTest::command(), utopia::make_unique<FEEvalTest>());
        // add_test(FSITest::command(), utopia::make_unique<FSITest>()); //FIXME too slow
        add_test(FormEvalTest::command(), utopia::make_unique<FormEvalTest>());
        // add_test(GeometryTest::command(), utopia::make_unique<GeometryTest>()); //FIXME missing mesh
        add_test(IntersectTest::command(), utopia::make_unique<IntersectTest>());
        // add_test(Intrepid2Test::command(), utopia::make_unique<Intrepid2Test>());
        add_test(MSHReaderTest::command(), utopia::make_unique<MSHReaderTest>());
        // add_test(MechTest::command(), utopia::make_unique<MechTest>()); //FIXME missing mesh
        add_test(NonLinearElasticityTest::command(), utopia::make_unique<NonLinearElasticityTest>());
        add_test(SMGTest::command(), utopia::make_unique<SMGTest>());
        add_test(RefactoredContactTest::command(), utopia::make_unique<RefactoredContactTest>());
        add_test(DualBasisTest::command(), utopia::make_unique<DualBasisTest>());
        add_test(IntegratorTest::command(), utopia::make_unique<IntegratorTest>());
        add_test(FETensorTest::command(), utopia::make_unique<FETensorTest>());
        add_test(NewNeohookeanTest::command(), utopia::make_unique<NewNeohookeanTest>());
    }

}