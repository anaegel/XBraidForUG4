//
// Created by parnet on 12.05.19.
//


#include "bridge/util.h"
#include "bridge/util_domain_algebra_dependent.h"
#include "SpaceTimeCommunicator.h"
#include "Scriptor.h"
#include "Talasma.h"
#include "PBraid.h"

#include "fieldgenerator.h"

using namespace std;
using namespace ug::bridge;

namespace ug{

    namespace XBraidForUG4 {

        struct Functionality {


            template<typename TDomain, typename TAlgebra>
            static void DomainAlgebra(Registry &reg, string grp) {
                string suffix = GetDomainAlgebraSuffix<TDomain, TAlgebra>();
                string tag = GetDomainAlgebraTag<TDomain, TAlgebra>();

                typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
                typedef SmartPtr<ug::GridFunction<TDomain, TAlgebra> > SPGridFunction;
                typedef SmartPtr<ug::UserData<double, TGridFunction::dim> > SPData;
                typedef SmartPtr<XBraidForUG4::SpaceTimeCommunicator> SPXCommunicator;
                typedef SmartPtr<ug::VTKOutput<TDomain::dim> > SPVTKOutput;



                typedef Scriptor<TDomain, TAlgebra> TScriptor;
                typedef SmartPtr<Scriptor<TDomain, TAlgebra> > SPScriptor;
                string name_scr = string("Scriptor").append(suffix);
                reg.add_class_<TScriptor>(name_scr, grp);
                reg.add_class_to_group(name_scr, "Scriptor", tag);


                typedef MultiScriptor<TDomain, TAlgebra> TMultiScriptor;
                string name_multi = string("MultiScriptor").append(suffix);
                reg.add_class_<TMultiScriptor,TScriptor>(name_multi, grp)
                        .add_constructor()
                        .add_method("addScriptor", &TMultiScriptor::addScriptor, "None", "a derived class of Scriptor", "adds a Scriptor to the storage")
                        .add_method("print", &TMultiScriptor::print, "None", "Gridfunction u#index#timestamp", "accesses a given gridfunction and delegate it to the stored Scriptors")
                        .add_method("write_time_pvd", &TMultiScriptor::write_time_pvd, "None", "Gridfunction u", "finishes the access and write a pvd")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_multi, "MultiScriptor", tag);

                typedef EvalScriptor<TDomain, TAlgebra> TEvalScriptor;
                string name_eval = string("EvalScriptor").append(suffix);
                reg.add_class_<TEvalScriptor,TScriptor>(name_eval, grp)
                        .add_constructor()
                        .add_method("setFile", static_cast<void (TEvalScriptor::*)(const char *)>(&TEvalScriptor::setFile), "void","filename","Creates a file")
                        .add_method("setGeneratorComponent", &TEvalScriptor::setGeneratorComponent, "None", "function component", "Set the function component for interpolation")
                        .add_method("setVectorGenerator", static_cast<void (TEvalScriptor::*)(const char*)>(&TEvalScriptor::setVectorGenerator), "None", "lua function name", "Set a LUA function to generate Gridfunctions to different timestamps")
                        .add_method("setVectorGenerator", static_cast<void (TEvalScriptor::*)(SPData)>(&TEvalScriptor::setVectorGenerator), "None", "userdata", "Set user data to generate Gridfunctions")
                        .add_method("setDomain", &TEvalScriptor::setDomain, "None", "domain discretization", "Set the Domain in which the Gridfunctions are generated")
                        .add_method("setRelative", &TEvalScriptor::setRelative, "None", "display norm of u and norm of norm(u-v)/nrom(u)", "")
                        .add_method("print", &TEvalScriptor::print, "None", "Gridfunction u# index # timestamp", "generates a Gridfunction and compare this with the given Gridfunction")
                        .add_method("write_time_pvd", &TEvalScriptor::write_time_pvd, "None", "", "finishes the access")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_eval, "EvalScriptor", tag);




                typedef VTKModScriptor<TDomain, TAlgebra> TVTKModScriptor;
                string name_vtkmod = string("VTKModScriptor").append(suffix);
                reg.add_class_<TVTKModScriptor,TScriptor>(name_vtkmod, grp)
                        .template add_constructor<void (*)(SPVTKOutput, const char *)>("")
                                .add_method("print", &TVTKModScriptor::print, "None","Gridfunction u# index # timestamp","writes a VTK file for the grid function")
                                .add_method("setModal", &TVTKModScriptor::setModal, "None","n","write only every n'th value")
                                .add_method("write_time_pvd", &TVTKModScriptor::write_time_pvd, "None","Gridfunction u","creates the PVD for the series")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_vtkmod, "VTKModScriptor", tag);
                {
                    typedef VTKScriptor<TDomain, TAlgebra> TVTKScriptor;
                    string name_vtk = string("VTKScriptor").append(suffix);
                    reg.add_class_<TVTKScriptor, TScriptor>(name_vtk, grp)
                            .template add_constructor<void (*)(SPVTKOutput, const char *)>("")
                            .add_method("print", &TVTKScriptor::print, "None", "Gridfunction u# index # timestamp",
                                        "writes a VTK file for the grid function")
                            .add_method("write_time_pvd", &TVTKScriptor::write_time_pvd, "None", "Gridfunction u",
                                        "creates the PVD for the series")
                            .set_construct_as_smart_pointer(true);;
                    reg.add_class_to_group(name_vtk, "VTKScriptor", tag);
                }

                {
                    typedef VTKIterationScriptor<TDomain, TAlgebra> TVTKScriptor;
                    string name_vtk = string("VTKIterationScriptor").append(suffix);
                    reg.add_class_<TVTKScriptor,TScriptor>(name_vtk, grp)
                            .template add_constructor<void (*)(SPVTKOutput, const char *)>("")
                            .add_method("print", &TVTKScriptor::print, "None","Gridfunction u# index # timestamp","writes a VTK file for the grid function")
                            .add_method("write_time_pvd", &TVTKScriptor::write_time_pvd, "None","Gridfunction u","creates the PVD for the series")
                            .set_construct_as_smart_pointer(true);;
                    reg.add_class_to_group(name_vtk, "VTKIterationScriptor", tag);
                }



                typedef UG4BraidApp<TDomain, TAlgebra> TBraidApp;
                string name_gf = string("UG4BraidApp").append(suffix);
                reg.add_class_<TBraidApp>(name_gf,grp)
                        .add_method("setVerbose", &TBraidApp::setVerbose, "None","verbose","set the level of verbose (true / false)")
                        .add_method("setStartTime", &TBraidApp::setStartTime, "None", "initial time", "set t0 as initial time")
                        .add_method("setEndTime", &TBraidApp::setEndTime, "None", "end time", "set tN as endtime")
                        .add_method("setNumberOfTimesteps", &TBraidApp::setNumberOfTimesteps, "", "number of timesteps", "set N as number of timesteps")
                        .add_method("setTimeValues", static_cast<void (TBraidApp::*)(double,double,int)>(&TBraidApp::setTimeValues), "None", "t0#tN#N", "sets tstart, tstop, number of timesteps")
                        .add_method("setDomainDisc", &TBraidApp::setDomainDisc, "None", "domain discretization", "set the domain")
                        .add_method("setStartVector", &TBraidApp::setStartVector, "None", "Gridfunction u0", "set the vector for t=t0")
                        .add_method("setMaxLevels", &TBraidApp::setMaxLevels, "None", "number", "set maximum number of level")
                        .add_method("setGeneratorComponent", &TBraidApp::setGeneratorComponent, "None", "function component", "set the function component for interpolation")
                        .add_method("setVectorGenerator", static_cast<void (TBraidApp::*)(const char*)>(&TBraidApp::setVectorGenerator), "None", "LuaFunction name", "set vector for time t != t0 as guess to speed up iterations")
                        .add_method("setVectorGenerator", static_cast<void (TBraidApp::*)(SPData)>(&TBraidApp::setVectorGenerator), "None", "UserData", "set vector for time t != t0 as guess to speed up iterations")
                        .add_method("setVectorGenerator", static_cast<void (TBraidApp::*)(ug::LuaFunctionHandle)>(&TBraidApp::setVectorGenerator), "None", "LuaFunctionHandle", "set vector for time t != t0 as guess to speed up iterations");
                reg.add_class_to_group(name_gf,"UG4BraidApp",tag);



                typedef RGFBraidApp<TDomain, TAlgebra> TRGFBraidApp;


                // SingleStageUniform
                string name_rgf = string("RGFBraidApp").append(suffix);
                reg.add_class_<TRGFBraidApp,TBraidApp>(name_rgf, grp)
                        .add_constructor()
                        .add_method("setTightTol", &TRGFBraidApp::setTightTol, "None", "tight tol", "sets the tight tolerance for adaptive convergence check")
                        .add_method("setLooseTol", &TRGFBraidApp::setLooseTol, "None", "loose tol", "sets the loose tolerance for adaptive convergence check")
                        .add_method("setAdaptConv", &TRGFBraidApp::setAdaptConv, "None", "adaptivity", "set adaptive convergence check")
                        .add_method("setTimeDisc", &TRGFBraidApp::setTimeDisc, "None", "single stage time discretization", "sets a single stage time discretization")
                        .add_method("setLinearSolver", &TRGFBraidApp::setLinearSolver, "None", "linear solver", "sets the linear Solver")
                        .add_method("setForceConv", &TRGFBraidApp::setForceConv, "None", "force convergence", "set cancel on convergence check fail")
                        .add_method("setStrongFirstIteration", &TRGFBraidApp::setStrongFirstIteration, "None", "value", "apply the first iteration with tight tolerance")
                        .add_method("setRecurringStrongIteration", &TRGFBraidApp::setRecurringStrongIteration, "None", "interval", "apply every n'th iteration with tight tolerance")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_rgf, "RGFBraidApp", tag);


                // SingeStageLevelDependend
                typedef MGFBraidApp<TDomain, TAlgebra> TMGFBraidApp;
                string name_mgf = string("MGFBraidApp").append(suffix);
                reg.add_class_<TMGFBraidApp,TBraidApp>(name_mgf, grp)
                        .add_constructor()
                        .add_method("setTimeDisc", &TMGFBraidApp::setTimeDisc, "None", "single stage time discretization", "sets a single stage time discretization")
                        .add_method("setLinearSolver", &TMGFBraidApp::setLinearSolver, "None", "linear solver", "sets the linear Solver")
                        .add_method("setForceConv", &TMGFBraidApp::setForceConv, "None", "force convergence", "cancle operation if convergence check fail")
                        .add_method("setStoreOperator", &TMGFBraidApp::setStoreOperator, "None", "level 0 is finest level", "store the operator for all level < x ")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_mgf, "MGFBraidApp", tag);


                //todo MultiStageUniform

                // MultiStageLevelDependend
                typedef ITSGFBraidApp<TDomain, TAlgebra> TITSGFBraidApp;
                string name_itsgf = string("ITSGFBraidApp").append(suffix);
                reg.add_class_<TITSGFBraidApp,TBraidApp>(name_itsgf, grp)
                        .add_constructor()
                        .add_method("setTimeDisc", &TITSGFBraidApp::setTimeDisc, "None", "timediscretization#level", "sets an arbitrary time discretization for level l ")
                        .add_method("setLinearSolver", &TITSGFBraidApp::setLinearSolver, "None", "linear solver", "set a linear Solver for level l")
                        .add_method("setSolver", &TITSGFBraidApp::setSolver, "None", "non linear solver", "set a non linear solver for level l")
                        .add_method("setForceConv", &TITSGFBraidApp::setForceConv, "None", "force", "cancel if convergence check fail")
                        .add_method("setStoreOperator", &TITSGFBraidApp::setStoreOperator, "None", "store", "store operator for level l")
                        .add_method("setResidual", &TITSGFBraidApp::setResidual, "None", "method", "set the residual method")
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name_itsgf, "ITSGFBraidApp", tag);

                // todo IterationDependendTimeDiscretization
                // todo completelua


                typedef SmartPtr<XBraidForUG4::SpaceTimeCommunicator> SPXCommunicator;
                typedef SmartPtr<UG4BraidApp<TDomain, TAlgebra> > SPBraidApp;
                typedef PBraid<TDomain, TAlgebra> TBraid;
                string name = string("Braid").append(suffix);
                reg.add_class_<TBraid>(name, grp)
                        .template add_constructor<void (*)(SPXCommunicator, SPBraidApp)>("Function(s)#Subset(s)") // todo parameter!
                        .add_method("run", static_cast<int (TBraid::*)(void)>(&TBraid::run), "", "", "starts XBraid iteration")
#ifdef UG_FOR_LUA
                        .add_method("run", static_cast<int (TBraid::*)(SPGridFunction, const char *, const char *, SPScriptor )>(&TBraid::run), "", "", "starts XBraid iteration")
                        .add_method("test", static_cast<int (TBraid::*)(SPGridFunction, const char *, const char *, SPScriptor )>(&TBraid::test), "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
#endif
                        .add_method("setApp", &TBraid::setApp, "", "", "set vector for time t = t0")
                        .add_method("getApp", &TBraid::getApp, "braid application", "", "set vector for time t = t0")

                        // .add_method("toggleResidual", &TBraid::toggleResidual, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setResidual", &TBraid::setResidual, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setMaxIterations", &TBraid::setMaxIterations, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setNRelax", &TBraid::setNRelax, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setCFactor", &TBraid::setCFactor, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setMaxLevels", &TBraid::setMaxLevels, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setSkip", &TBraid::setSkip, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setMinCoarse", &TBraid::setMinCoarse, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setAbsoluteTol", &TBraid::setAbsoluteTol, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setReduction", &TBraid::setReduction, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setTemporalNorm", &TBraid::setTemporalNorm, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setSequential", &TBraid::setSequential, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setStoreValues", &TBraid::setStoreValues, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setSpatialCoarsenAndRefine", &TBraid::setSpatialCoarsenAndRefine, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setRefine", &TBraid::setRefine, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setMaxRefinements", &TBraid::setMaxRefinements, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setAccessLevel", &TBraid::setAccessLevel, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setPrintLevel", &TBraid::setPrintLevel, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setPrintFile", &TBraid::setPrintFile, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setCycleFMG", &TBraid::setCycleFMG, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setCycleNFMG", &TBraid::setCycleNFMG, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setCycleNFMGV", &TBraid::setCycleNFMGV, "", "", "tests XBraid core functionality, grid has to be defined at time t=0.0")
                        .add_method("setOutput", &TBraid::createAccess, "", "", "set a VTK Output")
                        .add_method("setFilename", &TBraid::setFilename, "", "", "set a filename for VTK Output") // todo delete
                        .set_construct_as_smart_pointer(true);
                reg.add_class_to_group(name, "Braid", tag);


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
                {
                    typedef SinSourceOneCube<number, dim> T;
                    typedef CplUserData<number, dim> TBase;
                    string name = string("SinSourceOneCube").append(suffix);
                    reg.add_class_<T, TBase>(name, grp)
                            .add_constructor()
                            .add_method("setAlpha", &T::setAlpha, "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "SinSourceOneCube", tag);
                }
                {
                    typedef SinAnalyticSolutionOneCube<number, dim> T;
                    typedef CplUserData<number, dim> TBase;
                    string name = string("SinAnalyticSolutionOneCube").append(suffix);
                    reg.add_class_<T, TBase>(name, grp)
                            .add_constructor()
                            .add_method("setAlpha", &T::setAlpha, "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "SinAnalyticSolutionOneCube", tag);
                }
                {
                    typedef SinSourcePiCube<number, dim> T;
                    typedef CplUserData<number, dim> TBase;
                    string name = string("SinSourcePiCube").append(suffix);
                    reg.add_class_<T, TBase>(name, grp)
                            .add_constructor()
                            .add_method("setAlpha", &T::setAlpha, "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "SinSourcePiCube", tag);
                }
                {
                    typedef SinAnalyticSolutionPiCube<number, dim> T;
                    typedef CplUserData<number, dim> TBase;
                    string name = string("SinAnalyticSolutionPiCube").append(suffix);
                    reg.add_class_<T, TBase>(name, grp)
                            .add_constructor()
                            .add_method("setAlpha", &T::setAlpha, "", "")
                            .set_construct_as_smart_pointer(true);
                    reg.add_class_to_group(name, "SinAnalyticSolutionPiCube", tag);
                }
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
        using namespace TraceTools;
        using XBraidForUG4::SpaceTimeCommunicator;
        grp.append("XBraidForUG4");
        reg->add_class_<SpaceTimeCommunicator>("SpaceTimeCommunicator","XBraid","")
                .add_method("split",&SpaceTimeCommunicator::split)
                .add_method("unsplit",&SpaceTimeCommunicator::unsplit)

                .add_method("getGlobalRank",&SpaceTimeCommunicator::getGlobalRank)
                .add_method("getSpatialRank",&SpaceTimeCommunicator::getSpatialRank)
                .add_method("getTemporalRank",&SpaceTimeCommunicator::getTemporalRank)

                .add_method("getGlobalSize",&SpaceTimeCommunicator::getGlobalSize)
                .add_method("getSpatialSize",&SpaceTimeCommunicator::getSpatialSize)
                .add_method("getTemporalSize",&SpaceTimeCommunicator::getTemporalSize)

                .add_constructor()
                .set_construct_as_smart_pointer(true);

        reg->add_class_<Talasma>("Talasma", grp, "Class to measure time differences")
                .add_method("start", &Talasma::start)
                .add_method("stop", &Talasma::stop)
                .add_method("get", &Talasma::get)
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
