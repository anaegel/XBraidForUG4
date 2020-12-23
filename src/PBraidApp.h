//
// Created by parnet on 12.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_PBRAIDAPP_H
#define UG_PLUGIN_XBRAIDFORUG4_PBRAIDAPP_H



// UG4 lib.
#include "ugbase.h"
#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "common/serialization.h"
#include "lib_disc/io/vtkoutput.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_algebra/vector_interface/vec_functions.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "../../plugins/Limex/time_disc/time_integrator.hpp"

// This plugin.
#include "UG4BraidApp.h"

template<typename TDomain, typename TAlgebra>
class PBraidApp : public UG4BraidApp<TDomain, TAlgebra> {

public:
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;

    typedef ug::ThetaTimeStep<TAlgebra> TTimeStep;
    typedef SmartPtr<TTimeStep> SPTimeStep;

    typedef ug::VectorTimeSeries<typename TAlgebra::vector_type> TTimeSeries;
    typedef SmartPtr<TTimeSeries> SPTimeSeries;

    typedef SmartPtr<ug::ILinearOperatorInverse<typename TAlgebra::vector_type> > SPSolver;

    typedef SmartPtr<ug::AssembledLinearOperator<TAlgebra> > SPAssembledOperator;

    typedef SmartPtr<ug::ITimeIntegrator<TDomain, TAlgebra> > SPTimeIntegrator;

    typedef SmartPtr<ug::VTKOutput<TDomain::dim> > SPOutput;

    typedef typename ug::XBraidForUG4::SpaceTimeCommunicator TSpaceTimeCommunicator;

    SPTimeStep timeDisc;
    SPSolver linSolver;
    SPTimeIntegrator m_spIntegratorC;
    SPTimeIntegrator m_spIntegratorF;
    SPOutput out;


    SPGridFunction ux; // for t > tstart
    SmartPtr<TSpaceTimeCommunicator> comm;

    const char *filename = "FixedFilename";

    std::string pointerToString(void *c) {
        std::stringstream ss;
        ss << c;
        return ss.str();
    }

    bool write(TGridFunction *u, int index, double time) {
        //std::cout << "write:\t" << index << "\t @ " << time << std::endl;
        out->print(filename, *u, index, time);
        return true;
    }

    bool write(TGridFunction *u, int index, double time, const char *type) {
        std::stringstream ss;
        ss << filename;
        ss << "_T";
        ss << this->comm->getTemporalRank();
        ss << "_";
        ss << type;
        ss << "_";

        //std::cout << "write:\t" << index << "\t @ " << time << std::endl;
        out->print(ss.str().c_str(), *u, index, time);
        return true;
    }



    // bool verbose = true; todo

public:

    PBraidApp(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) :
            UG4BraidApp<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {}

    ~PBraidApp() override = default;


    void setStartVector(SPGridFunction p_u0) {
        this->u0 = p_u0;
    }

    void setRemainingVector(SPGridFunction p_ux) {
        this->ux = p_ux;
    }

    void setTimeDisc(SPTimeStep p_timeDisc) {
        this->timeDisc = p_timeDisc;
    }

    void setLinSolver(SPSolver p_linSolver) {
        this->linSolver = p_linSolver;
    }


    double assembled_dt = 0;
    SmartPtr<ug::AssembledLinearOperator<TAlgebra> > A;
    //SPGridFunction defect;

    void init() {
        const ug::GridLevel gridlevel = this->u0->grid_level();
        A = SmartPtr<ug::AssembledLinearOperator<TAlgebra> >(
                new ug::AssembledLinearOperator<TAlgebra>(timeDisc, gridlevel));
    }

    //todo change for user function
    braid_Int Init(braid_Real t, braid_Vector *u_ptr) override;


    /**
    input vector u corresponding to time tstart
    input ustop previous approximate solution at tstop
    input fstop Additional source at time tstop. (NULL?)
    output vector u the computed result for time tstop.

        virtual braid_Int Step(braid_Vector     u_,
                          braid_Vector     ustop_,
                          braid_Vector     fstop_,
                          BraidStepStatus &pstatus) = 0;

                          */
    virtual braid_Int Step(braid_Vector u, //
                           braid_Vector ustop, // estimated solution?
                           braid_Vector fstop,
                           BraidStepStatus &pstatus) override;
protected:
    int step_ignoring_residual(braid_Vector u, //
                               braid_Vector ustop, // estimated solution?
                               braid_Vector fstop,
                               BraidStepStatus &pstatus) {
        double tstart;
        double tstop;
        pstatus.GetTstartTstop(&tstart, &tstop);
        double current_dt = tstop - tstart;

        auto *sp_u_start = (SPGridFunction *) u->value;
        auto *sp_u_stop = (SPGridFunction *) u->value;
        auto *sp_u_stop_approx = (SPGridFunction *) ustop->value;

        auto *sp_rhs = new SPGridFunction(new TGridFunction(*this->u0));
        /*std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << pointerToString((void*)sp_rhs) << std::endl;
        std::cout << pointerToString((void*)sp_u_start) << std::endl;
        std::cout << pointerToString((void*)sp_u_stop) << std::endl;
        std::cout << pointerToString((void*)sp_u_stop_approx) << std::endl;*/



        //auto *rhs = (BraidVector *) malloc(sizeof(BraidVector));
        //auto *sp_rhs_val_ref = new SPGridFunction(new TGridFunction(*u0));
        //rhs->value = sp_rhs_val_ref;
        //fstop = rhs;

        const ug::GridLevel gridlevel = sp_u_stop->get()->grid_level();


        SPTimeSeries series = SPTimeSeries(new TTimeSeries());
        series->push(sp_u_start->get()->clone(), tstart);
        timeDisc->prepare_step(series, current_dt);

        if (current_dt != assembled_dt) { // todo is close?

            timeDisc->assemble_linear(*A, *sp_rhs->get(), gridlevel);
            linSolver->init(A, *sp_u_stop_approx->get());
            assembled_dt = current_dt;
        } else {
            timeDisc->assemble_rhs(*sp_rhs->get(), gridlevel);
        }

        /*if (fstop != nullptr) {
            auto *sp_rhs_stop = (SPGridFunction *) fstop->value;

            write(sp_u_stop_approx->get(),7055,0.01,"u_stop_approx");
            write(sp_u_start->get(),7055,0.01,"u_start");
            write(sp_rhs->get(),75055,0.01, "rhs");
            write(sp_rhs_stop->get(),75055,0.01,"rhs_stop");
            write(sp_u_stop->get(),75055,0.01,"u_stop");
            exit(80);
            VecAdd(1,*sp_u_stop_approx->get(),1,*sp_rhs_stop->get());
        }*/

        bool success = linSolver->apply(*sp_u_stop_approx->get(), *sp_rhs->get());

        if (!success) {
            std::cout << "Failure" << std::endl;
            exit(127);
        }


        /*std::cout << pointerToString((void*)sp_rhs) << std::endl;
        std::cout << pointerToString((void*)sp_u_start) << std::endl;
        std::cout << pointerToString((void*)sp_u_stop) << std::endl;
        std::cout << pointerToString((void*)sp_u_stop_approx) << std::endl;*/

        free(sp_rhs);

        if (sp_u_stop != sp_u_stop_approx) {
            free(sp_u_stop);
            u->value = new SPGridFunction(new TGridFunction(*sp_u_stop_approx->get()));
        }
        //sp_u_stop_approx; // solution

        return 0;
    }
public:

    braid_Int Residual(braid_Vector u, braid_Vector r,
                       BraidStepStatus &pstatus) override;

    // todo replace?
    braid_Int SpatialNorm(braid_Vector u, braid_Real *norm_ptr) override;
    braid_Int Access(braid_Vector u, BraidAccessStatus &astatus) override;


    //todo check
    braid_Int Coarsen(braid_Vector fu, braid_Vector *cu, BraidCoarsenRefStatus &status) override;

    //todo check
    braid_Int Refine(braid_Vector cu, braid_Vector *fu, BraidCoarsenRefStatus &status) override;


};


#include "PBraidApp_impl.hh"

#endif //UG_PLUGIN_XBRAIDFORUG4_PBRAIDAPP_H
