//
// Created by parnet on 29.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_MGFBRAIDAPP_H
#define UG_PLUGIN_XBRAIDFORUG4_MGFBRAIDAPP_H

#include <ugbase.h>

#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "common/serialization.h"
#include "lib_disc/io/vtkoutput.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_algebra/vector_interface/vec_functions.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "../../plugins/Limex/time_disc/time_integrator.hpp"

#include "GFBraidApp.h"

template<typename TDomain, typename TAlgebra>
class MRGFBraidApp : public GFBraidApp<TDomain, TAlgebra> {

public:


    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;

    typedef SmartPtr<TGridFunction> SPGridFunction;

    typedef ug::ThetaTimeStep<TAlgebra> TTimeStep;

    typedef SmartPtr<ug::ThetaTimeStep<TAlgebra>> SPTimeStep;

    typedef ug::VectorTimeSeries<typename TAlgebra::vector_type> TTimeSeries;

    typedef SmartPtr<ug::VectorTimeSeries<typename TAlgebra::vector_type>> SPTimeSeries;

    typedef SmartPtr<ug::ILinearOperatorInverse<typename TAlgebra::vector_type>> SPSolver;

    typedef SmartPtr<ug::AssembledLinearOperator<TAlgebra>> SPAssembledOperator;

    typedef SmartPtr<ug::ITimeIntegrator<TDomain, TAlgebra>> SPTimeIntegrator;

    typedef SmartPtr<ug::VTKOutput<TDomain::dim>> SPOutput;

    SPTimeIntegrator m_spIntegratorC;
    SPTimeIntegrator m_spIntegratorF;

    SPOutput out;
    SPGridFunction ux; // for t > tstart

    struct LevelSolver {
        double assembled_dt;
        SPAssembledOperator A;
        SPTimeStep timeDisc;
        SPSolver linSolver;
    };

    size_t levels = 1;
    LevelSolver *lvl = new LevelSolver[levels];

    const char *filename{};

    std::string pointerToString(void *c) {
        std::stringstream ss;
        ss << c;
        return ss.str();
    }

    bool write(TGridFunction *u, int index, double time) {
        std::cout << "write:\t" << index << "\t @ " << time << std::endl;
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

        std::cout << "write:\t" << index << "\t @ " << time << std::endl;
        out->print(ss.str().c_str(), *u, index, time);
        return true;
    }



    // bool verbose = true; todo

public:

    MRGFBraidApp(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) : GFBraidApp<TDomain, TAlgebra>(
            mpi_temporal, tstart, tstop,
            steps) {}

    ~MRGFBraidApp() override {
        delete[] this->lvl;
    };


    void setStartVector(SPGridFunction p_u0) {
        this->u0 = p_u0;
    }

    void setRemainingVector(SPGridFunction p_ux) {
        this->ux = p_ux;
    }

    void setTimeDisc(SPTimeStep p_timeDisc, size_t level = 0) {
        this->lvl[level].timeDisc = p_timeDisc;
    }

    void setLinSolver(SPSolver p_linSolver, size_t level = 0) {
        this->lvl[level].linSolver = p_linSolver;
    }


    //SPGridFunction defect;
    void setLevels(size_t levelcount) {
        delete[] this->lvl;
        this->levels = levelcount;
        this->lvl = new LevelSolver[this->levels];
    }


    void init() {
        const ug::GridLevel gridlevel = this->u0->grid_level();
        for (size_t i = 0; i < this->levels; i++) {
            this->lvl[i].A = SmartPtr<ug::AssembledLinearOperator<TAlgebra>>(
                    new ug::AssembledLinearOperator<TAlgebra>(this->lvl[i].timeDisc, gridlevel));
        }
    }

    //todo change for user function
    braid_Int Init(braid_Real t, braid_Vector *u_ptr)
    override {
        auto *u = (BraidVector *) malloc(sizeof(BraidVector));
        SPGridFunction *vec;
        if (t == this->tstart) {
            vec = new SPGridFunction(new TGridFunction(*this->u0));
        } else {
            vec = new SPGridFunction(new TGridFunction(*this->ux));
        }
        u->value = vec;
        *u_ptr = u;
        return 0;
    };


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
                           BraidStepStatus &pstatus) {
        double tstart;
        double tstop;
        int current_level; // current level

        pstatus.GetTstartTstop(&tstart, &tstop);
        pstatus.GetLevel(&current_level);

        double current_dt = tstop - tstart;

        auto *sp_u_start = (SPGridFunction *) u->value;
        auto *sp_u_stop = (SPGridFunction *) u->value;
        auto *sp_u_stop_approx = (SPGridFunction *) ustop->value;
        auto *sp_rhs = new SPGridFunction(new TGridFunction(*this->u0));

        /* Show Pointer Adresses
         std::cout << pointerToString((void*)sp_u_start) << std::endl;
         std::cout << pointerToString((void*)sp_u_stop) << std::endl;
         std::cout << pointerToString((void*)sp_u_stop_approx) << std::endl;
         std::cout << pointerToString((void*)sp_rhs) << std::endl;
         */



        //auto *rhs = (BraidVector *) malloc(sizeof(BraidVector));
        //auto *sp_rhs_val_ref = new SPGridFunction(new TGridFunction(*u0));
        //rhs->value = sp_rhs_val_ref;
        //fstop = rhs;

        const ug::GridLevel gridlevel = sp_u_stop->get()->grid_level();


        SPTimeSeries series = SPTimeSeries(new TTimeSeries());
        series->push(sp_u_start->get()->clone(), tstart);
        this->lvl[current_level].timeDisc->prepare_step(series, current_dt);

        if (current_dt != this->lvl[current_level].assembled_dt) { // todo is close?

            this->lvl[current_level].timeDisc->assemble_linear(*this->lvl[current_level].A, *sp_rhs->get(), gridlevel);
            this->lvl[current_level].linSolver->init(this->lvl[current_level].A, *sp_u_stop_approx->get());
            this->lvl[current_level].assembled_dt = current_dt;
        } else {
            this->lvl[current_level].timeDisc->assemble_rhs(*sp_rhs->get(), gridlevel);
        }

        if (fstop != nullptr) {
            auto *sp_rhs_stop = (SPGridFunction *) fstop->value;
            /*
            write(sp_u_stop_approx->get(),7055,0.01,"u_stop_approx");
            write(sp_u_start->get(),7055,0.01,"u_start");
            write(sp_rhs->get(),75055,0.01, "rhs");
            write(sp_rhs_stop->get(),75055,0.01,"rhs_stop");
            write(sp_u_stop->get(),75055,0.01,"u_stop");
            */
            //VecAdd(1,*sp_rhs->get(),1,*sp_rhs_stop->get());
            VecAdd(1, *sp_rhs->get(), -1, *sp_rhs_stop->get());
            //exit(80);
        }

        bool success = this->lvl[current_level].linSolver->apply(*sp_u_stop_approx->get(), *sp_rhs->get());

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
    };

    braid_Int Residual(braid_Vector u, braid_Vector r,
                       BraidStepStatus &pstatus) override {

        // std::cout << this->comm->getGlobalRank() << "\tResidual requested" << std::endl;
        auto *sp_u_start = (SPGridFunction *) u->value; // input
        auto *sp_u_stop = (SPGridFunction *) r->value; // input
        auto *sp_r_stop = (SPGridFunction *) r->value; // output == u_stop!

        double t_start, t_stop;
        int current_level;
        pstatus.GetTstartTstop(&t_start, &t_stop);
        pstatus.GetLevel(&current_level);
        double current_dt = t_stop - t_start;

        const ug::GridLevel gridlevel = sp_u_start->get()->grid_level();

        SmartPtr<TTimeSeries> series = SmartPtr<TTimeSeries>(new TTimeSeries());

        series->push(*sp_u_start, this->tstart); //todo clone sp_current_u_ref ?
        this->lvl[current_level].timeDisc->prepare_step(series, current_dt);
        auto *sp_rhs = new SPGridFunction(new TGridFunction(*this->u0));
        if (current_dt != this->lvl[current_level].assembled_dt) { // todo check for is close?
            this->lvl[current_level].timeDisc->assemble_linear(*this->lvl[current_level].A, *sp_rhs->get(), gridlevel);
            this->lvl[current_level].linSolver->init(this->lvl[current_level].A, *sp_u_stop->get()); // todo linearization point ?
            this->lvl[current_level].assembled_dt = current_dt;
        } else {
            this->lvl[current_level].timeDisc->assemble_rhs(*sp_rhs->get(), gridlevel);
        }

        this->lvl[current_level].linSolver->linear_operator()->apply_sub(
                *sp_rhs->get(), // f co domain function [in / out]
                *sp_u_stop->get() // u domain function [in]
        );

        free(r->value);
        r->value = sp_rhs;

        //exit(77);
        return 0;

    };

    // todo replace?
    braid_Int SpatialNorm(braid_Vector u, braid_Real *norm_ptr)
    override {
        //a->setSpatialNorm(&ug::VecTwoNormSq<typename TAlgebra::vector_type>);
        //a->setSpatialNorm(&VecNorm2<typename TAlgebra::vector_type>);
        auto *uref = (SPGridFunction *) u->value;
        *norm_ptr = (*uref)->norm();
        return 0;
    };

    braid_Int Access(braid_Vector u, BraidAccessStatus &astatus)
    override {
        int v = 0;

        int index;
        astatus.GetTIndex(&index);
        double timestamp;
        astatus.GetT(&timestamp);

        auto *ref = (SPGridFunction *) u->value;
        v = write(ref->get(), index, timestamp);
        return v;

    };


    //todo check
    braid_Int Coarsen(braid_Vector fu, braid_Vector *cu, BraidCoarsenRefStatus &status)
    override {
        this->Clone(fu, cu);

        auto *sp_fu = (SPGridFunction *) fu->value;
        auto *sp_cu = (SPGridFunction *) (*cu)->value;

        double t_upper;
        double t_lower;
        status.GetT(&t_lower);
        status.GetCTstop(&t_upper); // todo C or F?
        //status.GetFTstop(&t_upper);

        m_spIntegratorC->apply(*sp_fu, t_upper, sp_cu->cast_const(), t_lower); // todo check fu, cu order
        return 0;
    }

    //todo check
    braid_Int Refine(braid_Vector cu, braid_Vector *fu, BraidCoarsenRefStatus &status)
    override {
        this->Clone(cu, fu);

        auto *sp_fu = (SPGridFunction *) (*fu)->value;
        auto *sp_cu = (SPGridFunction *) cu->value;

        double t_upper;
        double t_lower;
        status.GetT(&t_lower);
        status.GetCTstop(&t_upper); // todo C or F?
        //status.GetFTstop(&t_upper);

        m_spIntegratorF->apply(*sp_cu, t_upper, sp_fu->cast_const(), t_lower); // todo check fu, cu order

        return 0;
    }
}


#endif //UG_PLUGIN_XBRAIDFORUG4_MGFBRAIDAPP_H
