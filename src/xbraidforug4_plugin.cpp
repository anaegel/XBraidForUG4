//
// Created by parnet on 12.05.19.
//


#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"
#include "XCommunicator.h"
#include "PBraid.h"

using namespace std;
using namespace ug::bridge;

namespace ug{

    namespace XBraidForUG4 {

        struct Functionality {
            template<typename TDomain, typename TAlgebra>
            static void DomainAlgebra(Registry &reg, string grp) {
                string suffix = GetDomainAlgebraSuffix<TDomain, TAlgebra>();
                string tag = GetDomainAlgebraTag<TDomain, TAlgebra>();

                typedef PBraid<TDomain, TAlgebra> TBraid;
                string name = string("Braid").append(suffix);
                reg.add_class_<TBraid>(name, grp)
                        .add_constructor()
                        .add_method("run", &TBraid::run, "", "", "starts XBraid iteration")
                        .add_method("test", &TBraid::test, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setVerbose", &TBraid::setVerbose, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setStartVector", &TBraid::setStartVector, "", "", "set vector for time t = t0")
                        .add_method("setStartTime", &TBraid::setStartTime, "", "", "set vector for time t = t0")
                        .add_method("setEndTime", &TBraid::setEndTime, "", "", "set vector for time t = t0")
                        .add_method("setNumberOfTimesteps", &TBraid::setNumberOfTimesteps, "", "", "set vector for time t = t0")
                        .add_method("setRemainingVector", &TBraid::setRemainingVector, "", "", "set vector for time t != t0")
                        .add_method("setTimeDisc", &TBraid::setTimeDisc, "", "", "sets time discretization")
                        .add_method("setLinearSolver", &TBraid::setLinearSolver, "", "", "set linear Solver")
                        .add_method("setOutput", &TBraid::createAccess, "", "", "set a VTK Output")
                        .add_method("setFilename", &TBraid::setFilename, "", "", "set a filename for VTK Output")
                        .add_method("setCommunicator", &TBraid::setCommunicator, "", "", "Splits the MPI communicator for temporal and spatial parallization")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "XBraidForUG4", tag);

            }

            template<typename TDomain>
            static void Domain(Registry &reg, string grp) {
                string suffix = GetDomainSuffix<TDomain>();
                string tag = GetDomainTag<TDomain>();

            }

            template<int dim>
            static void Dimension(Registry &reg, string grp) {
                string suffix = GetDimensionSuffix<dim>();
                string tag = GetDimensionTag<dim>();

            }

            template<typename TAlgebra>
            static void Algebra(Registry &reg, string grp) {
                string suffix = GetAlgebraSuffix<TAlgebra>();
                string tag = GetAlgebraTag<TAlgebra>();


            }

            static void Common(Registry &reg, string grp) {
            }

        };
    } // end namespace XBraidForUG4

    extern "C" void
    InitUGPlugin_XBraidForUG4(Registry *reg, string grp) {
        using namespace XBraidForUG4;
        grp.append("XBraidForUG4");
        reg->add_class_<XCommunicator>("XCommunicator","XBraid","")
                .add_method("split",&XCommunicator::split)
                .add_method("setOutput",&XCommunicator::setOutput)

                .add_method("getGlobalRank",&XCommunicator::getGlobalRank)
                .add_method("getSpatialRank",&XCommunicator::getSpatialRank)
                .add_method("getTemporalRank",&XCommunicator::getTemporalRank)

                .add_method("getGlobalSize",&XCommunicator::getGlobalSize)
                .add_method("getSpatialSize",&XCommunicator::getSpatialSize)
                .add_method("getTemporalSize",&XCommunicator::getTemporalSize)

                .add_method("processPrint",&XCommunicator::processPrint)
                .add_method("print",&XCommunicator::print)
                .add_constructor()
                .set_construct_as_smart_pointer(true);
        try {
            RegisterCommon<Functionality>(*reg, grp);
            RegisterDimensionDependent<Functionality>(*reg, grp);
            RegisterDomainDependent<Functionality>(*reg, grp);
            RegisterAlgebraDependent<Functionality>(*reg, grp);
            RegisterDomainAlgebraDependent<Functionality>(*reg, grp);
        }
        UG_REGISTRY_CATCH_THROW(grp);

    }
}