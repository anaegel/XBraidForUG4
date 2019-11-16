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
class MGFBraidApp : public GFBraidApp<TDomain, TAlgebra>
{

public:

    /****************************************************************************
    * Typedefs
    ***************************************************************************/
    typedef GFBraidApp<TDomain, TAlgebra> base_type;
    typedef MGFBraidApp<TDomain, TAlgebra> this_type;

    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;

    typedef SmartPtr <TGridFunction> SPGridFunction;

    typedef ug::ThetaTimeStep <TAlgebra> TTimeStep;

    typedef SmartPtr<ug::ThetaTimeStep<TAlgebra> > SPTimeStep;

    typedef ug::VectorTimeSeries<typename TAlgebra::vector_type> TTimeSeries;

    typedef SmartPtr<ug::VectorTimeSeries<typename TAlgebra::vector_type> > SPTimeSeries;

    typedef SmartPtr<ug::ILinearOperatorInverse<typename TAlgebra::vector_type> > SPSolver;

    typedef SmartPtr<ug::AssembledLinearOperator<TAlgebra> > SPAssembledOperator;

    typedef ug::AssembledLinearOperator<TAlgebra> TAssembledOperator;

    typedef SmartPtr<ug::ITimeIntegrator<TDomain, TAlgebra> > SPTimeIntegrator;

    typedef SmartPtr<ug::VTKOutput<TDomain::dim> > SPOutput;


    struct LevelSolver {
        double assembled_dt = -1;
        bool storeOperator = false;
        bool forceConv = true;
        SPAssembledOperator A;
        SPTimeStep timeDisc;
        SPSolver linSolver;
    };


    SPAssembledOperator sharedA;
    double shared_assembled_dt;

    SPTimeSeries series = SPTimeSeries(new TTimeSeries());
    /****************************************************************************
    * Members
    ***************************************************************************/
    SPTimeIntegrator m_spIntegratorC;
    SPTimeIntegrator m_spIntegratorF;

    int m_progress = 0;

    /*
     * XBraid uses low index for fine grids
     */
    LevelSolver *m_lvl = new LevelSolver[this->m_levels];

public:
    /****************************************************************************
    * Constructor / Destructur
    ***************************************************************************/
    MGFBraidApp() : GFBraidApp<TDomain, TAlgebra>() {
        this->name = "leveldependend";
    }

    virtual ~MGFBraidApp() {
        delete[] this->m_lvl;
    };

    void setMaxLevels(size_t levelcount) {
        delete[] this->m_lvl;
        this->m_levels = levelcount;
        this->m_lvl = new LevelSolver[this->m_levels];
        //std::cout << "Max level is " << levelcount << std::endl;
    }


    /**
     * Assembled operator of coarse level ( level > x ) will not be stored.
     * To save all assembled operators use x == MaxLevel
     * @param level
     */
    void setStoreOperator(size_t x) {
        //std::cout << "store operator " << x<< std::endl;
        if (x > this->m_levels) {
            x = this->m_levels;
        }

        for (size_t i = 0; i < x; ++i) {
            m_lvl[i].storeOperator = true;
        }
        for (size_t i = x; i < this->m_levels; ++i) {
            m_lvl[i].storeOperator = false;
        }
    }

    void setTimeDisc(SPTimeStep p_timeDisc, size_t level = 0) {
        // std::cout << "set timeDisc for level " << level << std::endl;
        this->m_lvl[level].timeDisc = p_timeDisc;
    }

    void setLinearSolver(SPSolver p_linSolver, size_t level = 0) {
        // std::cout << "set linearSolver for level " << level << std::endl;
        this->m_lvl[level].linSolver = p_linSolver;
    }

    void setForceConv(bool pForceConvergence, size_t level = 0) {
        // std::cout << "set force conv for level " << level << "\t" << pForceConvergence << std::endl;
        this->m_lvl[level].forceConv = pForceConvergence;
    }

    void init() override {


        base_type::timer.start();
#if TRACE_TIMINGS == 1
        this->redoran = Redoran(this->m_levels);
#endif

#if TRACE_GRIDFUNCTION == 1
        std::stringstream ss;
        ss << "job_" << this->m_comm->getTemporalRank() << ".output";
        debugwriter.open(ss.str());
        this->matlab = SmartPtr<MATLABScriptor<TDomain, TAlgebra> >(new MATLABScriptor<TDomain, TAlgebra>(this->debugwriter));
#endif
        const ug::GridLevel gridlevel = this->m_u0->grid_level();
        for (size_t i = 0; i < this->m_levels; ++i) {
            if (this->m_lvl[i].storeOperator) {
                this->m_lvl[i].A = SPAssembledOperator(new TAssembledOperator(this->m_lvl[i].timeDisc, gridlevel));
            }
        }
    }


    braid_Int Step(braid_Vector u, //
                   braid_Vector ustop,
                   braid_Vector fstop,
                   BraidStepStatus &pstatus) override {
        int l; // level
        pstatus.GetLevel(&l);
        StartRedoranLevel(LevelObserver::T_STEP,l);


        //std::cout << "step" << std::endl;
        double t_start, t_stop;
        pstatus.GetTstartTstop(&t_start, &t_stop);

        double current_dt = t_stop - t_start;

        if (this->m_verbose) {
            int tindex;
            pstatus.GetTIndex(&tindex);
            int iteration;
            pstatus.GetIter(&iteration);
#if TRACE_INDEX == 1


            if (fstop == nullptr) {
                base_type::debugwriter << "u_" << u->index << " = step_"<<l<<"_n( u_" << u->index << ", u_" << ustop->index
                                  << ", null, " << t_start << ", " << current_dt << ", " << t_stop << ", " << l
                                  << ")"
                                  << "\t\t % " << tindex << std::endl;
            } else {
            	base_type::debugwriter << "u_" << u->index << " = step_"<<l<<"_r( u_" << u->index << ", u_" << ustop->index
                                  << ", u_"
                                  << fstop->index << ", " << t_start << ", " << current_dt << ", " << t_stop << ", "
                                  << l << ")"
                                  << " % " << tindex << std::endl;
            }

#else
            this->debugwriter << std::setw(13) << iteration << "step for level " << l << " at position " << tindex << " and iteration" << std::endl;
#endif
        }


        auto *sp_u_approx_tstart = (SPGridFunction *) u->value;
        auto *constsp_u_approx_tstop = (SPGridFunction *) ustop->value;
        SPGridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();
        SPGridFunction lp = constsp_u_approx_tstop->get()->clone();

        SPGridFunction sp_rhs = this->m_u0->clone_without_values(); // for rhs
        const ug::GridLevel gridlevel = sp_u_approx_tstart->get()->grid_level();

        //SPGridFunction last_time_step = sp_u_approx_tstart->get()->clone(); // todo delete?
        series->push(*sp_u_approx_tstart, t_start);
        this->m_lvl[l].timeDisc->prepare_step(series, current_dt);




        //this->debugwriter << "check operator" << std::endl;
        SPAssembledOperator *spA;
        bool reassemble = false;
        if (this->m_lvl[l].storeOperator) {
            spA = &this->m_lvl[l].A;
            if (this->m_lvl[l].assembled_dt == -1.0) {
                std::cout << "create operator for level " << l << std::endl;
                this->m_lvl[l].assembled_dt = current_dt;
                reassemble = true;
            }
        } else {
            spA = &this->sharedA;
            if (fabs(this->shared_assembled_dt - current_dt) > 1e-14) {
                std::cout << "create shared operator" << std::endl;
                shared_assembled_dt = current_dt;
                reassemble = true;
            }
        }
        if (reassemble) {

            //this->debugwriter << "assemble operator" << std::endl;
            this->m_lvl[l].timeDisc->assemble_linear(*spA->get(), *sp_rhs.get(), gridlevel);
            this->m_lvl[l].linSolver->init(*spA, *sp_u_tstop_approx.get()); // todo check if exists?
        } else {

            // this->debugwriter << "assemble vector" << std::endl;
            this->m_lvl[l].timeDisc->assemble_rhs(*sp_rhs.get(), gridlevel);
        }

        //this->debugwriter << "use fstop" << std::endl;
        if (fstop != nullptr) {
            auto *sp_fstop = (SPGridFunction *) fstop->value;
            //VecAdd(a,x,b,y) <=> x = a * x + b * y
            VecAdd(1, *sp_rhs.get(), 1, *sp_fstop->get());
        }


        //this->debugwriter << "apply operator" << std::endl;
        bool success = this->m_lvl[l].linSolver->apply(*sp_u_tstop_approx.get(), *sp_rhs.get());

        if (!success) {
            if (this->m_lvl[l].forceConv) {
                //std::cout << "Failure Convergence not reached" << std::endl;
                exit(127);
            }
        }

        *sp_u_approx_tstart = sp_u_tstop_approx;

        series->clear();
        //this->debugwriter << "finished step" << std::endl;
        return 0;
    };

    braid_Int Residual(braid_Vector u, braid_Vector r,
                       BraidStepStatus &pstatus) override {
#if TRACE_CONST == 1
        r->m_const = false;
        // u->m_const = false;
#endif
        int l; // level;
        pstatus.GetLevel(&l);
        StartRedoranLevel(LevelObserver::T_RESIDUAL,l);

        double t_start, t_stop;
        pstatus.GetTstartTstop(&t_start, &t_stop);
        double current_dt = t_stop - t_start;

        int tindex;
        pstatus.GetTIndex(&tindex);

#if TRACE_INDEX == 1
        if (this->m_verbose) {
            this->debugwriter << "u_" << r->index << " =  residual( u_" << u->index << " , u_" << r->index
                              << ", "
                              << t_start << ", " << current_dt << ", " << t_stop << ", " << l << ")"
                              << " % " << tindex << std::endl;
        }
#endif
        auto *const_u_approx_tstop = (SPGridFunction *) u->value;
        auto *u_approx_tstart = (SPGridFunction *) r->value;

        const ug::GridLevel gridlevel = const_u_approx_tstop->get()->grid_level();


        series->push(*u_approx_tstart, t_start);
        // r = utstart
        // addboundary(r)
        // addforce(tstop, r) [ r+= dt*F(x,y)
        this->m_lvl[l].timeDisc->prepare_step(series, current_dt);
        auto sp_rhs = this->m_u0->clone_without_values();

        SPAssembledOperator *spA;
        bool reassemble = false;
        if (this->m_lvl[l].storeOperator) {
            spA = &this->m_lvl[l].A;
            if (this->m_lvl[l].assembled_dt == -1.0) {
                std::cout << "create operator for level " << l << std::endl;
                this->m_lvl[l].assembled_dt = current_dt;
                reassemble = true;
            }
        } else {
            spA = &this->sharedA;
            if (fabs(shared_assembled_dt - current_dt) > 1e-14) {
                std::cout << "create shared operator" << std::endl;
                shared_assembled_dt = current_dt;
                reassemble = true;
            }
        }
        if (reassemble) {
            StartRedoranLevel(LevelObserver::T_ASSEMBLE_OP,l);
            this->m_lvl[l].timeDisc->assemble_linear(*spA->get(), *sp_rhs.get(), gridlevel);
            this->m_lvl[l].linSolver->init(*spA, *u_approx_tstart->get());
            StopRedoranLevel(LevelObserver::T_ASSEMBLE_OP,l);
        } else {
            StartRedoranLevel(LevelObserver::T_ASSEMBLE_RHS,l);
            this->m_lvl[l].timeDisc->assemble_rhs(*sp_rhs.get(), gridlevel);
            StopRedoranLevel(LevelObserver::T_ASSEMBLE_RHS,l);
        }


        this->m_lvl[l].linSolver->linear_operator()->apply_sub(
                *sp_rhs.get(), // f co domain function [in / out]
                *const_u_approx_tstop->get() // u domain function [in]
        ); // calculates r = r - A * u

        (*sp_rhs) *= -1; // r = -r + A*u
        *u_approx_tstart = sp_rhs;

        //SPGridFunction  def = sp_rhs->clone(); // todo alternative?
        //m_timeDisc->assemble_defect(*def.get(),*const_u_approx_tstop->get(),gridlevel);
        //VecAdd(1,*def.get(),-1,*u_approx_tstart->get());

#if TRACE_INDEX == 1
        MATLAB(sp_rhs.get(), u->index, t_stop);
#endif
        series->clear();
        StopRedoranLevel(LevelObserver::T_RESIDUAL,l);
        return 0;

    };


    braid_Int SpatialNorm(braid_Vector u, braid_Real *norm_ptr) override {
        auto *uref = (SPGridFunction *) u->value;
        // todo clone for non residual
        *norm_ptr = (*uref)->norm();
#if TRACE_INDEX == 1
        if (this->m_verbose) {
            this->debugwriter << "norm( u_" << u->index << ") % " << *norm_ptr << std::endl;
        }
#endif
        return 0;
    };


    braid_Int Coarsen(braid_Vector fu, braid_Vector *cu, BraidCoarsenRefStatus &status)
    override {
       /* Clone(fu, cu);

        auto *sp_fu = (SPGridFunction *) fu->value;
        auto *sp_cu = (SPGridFunction *) (*cu)->value;

        double t_upper;
        double t_lower;
        status.GetT(&t_lower);
        status.GetCTstop(&t_upper);
        //status.GetFTstop(&t_upper);

        m_spIntegratorC->apply(*sp_fu, t_upper, sp_cu->cast_const(), t_lower); // todo check fu, cu order*/
        return 0;
    }


    braid_Int Refine(braid_Vector cu, braid_Vector *fu, BraidCoarsenRefStatus &status)
    override {
       /* Clone(cu, fu);

        auto *sp_fu = (SPGridFunction *) (*fu)->value;
        auto *sp_cu = (SPGridFunction *) cu->value;

        double t_upper;
        double t_lower;
        status.GetT(&t_lower);
        status.GetCTstop(&t_upper);
        //status.GetFTstop(&t_upper);

        m_spIntegratorF->apply(*sp_cu, t_upper, sp_fu->cast_const(), t_lower); // todo check fu, cu order*/

        return 0;
    }
};


#endif //UG_PLUGIN_XBRAIDFORUG4_MGFBRAIDAPP_H
