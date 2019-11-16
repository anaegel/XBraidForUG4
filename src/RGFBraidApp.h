//
// Created by parnet on 24.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_RGFBRAIDAPP_H
#define UG_PLUGIN_XBRAIDFORUG4_RGFBRAIDAPP_H

#include <ugbase.h>

#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "common/serialization.h"

#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_algebra/vector_interface/vec_functions.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "../../plugins/Limex/time_disc/time_integrator.hpp"
#include "lib_disc/dof_manager/function_pattern.h"

#include "GFBraidApp.h"

template<typename TDomain, typename TAlgebra>
class RGFBraidApp : public GFBraidApp<TDomain, TAlgebra> {

public:

    /****************************************************************************
    * Typedefs
    ***************************************************************************/
    typedef GFBraidApp<TDomain, TAlgebra> base_type;

	typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;

    typedef SmartPtr<TGridFunction> SPGridFunction;

    typedef ug::ThetaTimeStep <TAlgebra> TTimeStep;

    typedef SmartPtr<ug::ThetaTimeStep<TAlgebra> > SPTimeStep;

    typedef ug::VectorTimeSeries<typename TAlgebra::vector_type> TTimeSeries;

    typedef SmartPtr<ug::VectorTimeSeries<typename TAlgebra::vector_type> > SPTimeSeries;

    typedef ug::ILinearOperatorInverse<typename TAlgebra::vector_type> TSolver;

    typedef SmartPtr<TSolver> SPSolver;

    typedef ug::StdConvCheck<typename TAlgebra::vector_type> TConv;

    typedef SmartPtr<TConv> SPConv;

    typedef ug::AssembledLinearOperator <TAlgebra> TAssembledOperator;

    typedef SmartPtr<ug::AssembledLinearOperator<TAlgebra> > SPAssembledOperator;

    typedef ug::ITimeIntegrator<TDomain, TAlgebra> TTimeIntegrator;

    typedef SmartPtr<TTimeIntegrator> SPTimeIntegrator;

    /****************************************************************************
    * Members
    ***************************************************************************/
    SPTimeStep m_timeDisc;
    SPSolver m_linSolver;
    double m_assembled_dt = 0;
    SPAssembledOperator m_A;

    SPTimeIntegrator m_spIntegratorC;
    SPTimeIntegrator m_spIntegratorF;


    int m_progress = 0;
    bool forceConv = false;
    bool adaptConv = false;

    bool strongFirstIteration = false;
    bool strongRecurring = false;
    int recurringInterval = 4;
    SPTimeSeries series = SPTimeSeries(new TTimeSeries());


public:

    void setAdaptConv(bool b_state){
        this->adaptConv = b_state;
    }

    void setStrongFirstIteration(bool b_state){
        this->strongFirstIteration = b_state;
    }

    void setRecurringStrongIteration(int interval){
        this->strongRecurring = true;
        this->strongFirstIteration = interval;
    }

    /****************************************************************************
    * Constructor / Destructur
    ***************************************************************************/
    RGFBraidApp() : GFBraidApp<TDomain, TAlgebra>() {
        this->name = "uniform";
    }

    ~RGFBraidApp() = default;

    /****************************************************************************
    * Functions to outsource // todo
    ***************************************************************************/

    void setTimeDisc(SPTimeStep p_timeDisc, size_t level = 0) {
        this->m_timeDisc = p_timeDisc;
    }

    void setLinearSolver(SPSolver p_linSolver, size_t level = 0) {
        this->m_linSolver = p_linSolver;
    }


    void init() override {

    	 std::stringstream ss;
    	 ss << "job_" << this->m_comm->getTemporalRank() << ".output";
    	base_type::debugwriter.open(ss.str());
    	base_type::timer.start();
#if TRACE_TIMINGS == 1
        this->redoran = Redoran(this->m_levels);
#endif



#if TRACE_GRIDFUNCTION == 1

        this->matlab = SmartPtr<MATLABScriptor<TDomain, TAlgebra> >(new MATLABScriptor<TDomain, TAlgebra>(this->debugwriter));
#endif
        const ug::GridLevel gridlevel = this->m_u0->grid_level();
        m_A = SPAssembledOperator(new TAssembledOperator(m_timeDisc, gridlevel));
        std::cout << "finished init" << std::endl;
    }

    void print_settings(){
        std::cout << "Settings: " << this->name << std::endl;
        std::cout << "Max number of level: " << this->m_levels << std::endl;
        std::cout << "Solver: " << this->m_linSolver->name() << std::endl;
        std::cout << "Force convergence: " <<  this->forceConv << std::endl;
        std::cout << "Use adaptive ConvCheck: " << this->adaptConv << std::endl;
        if(this->adaptConv){
            std::cout << "\t loose: " << this->loose_tol << std::endl;
            std::cout << "\t tight: " << this->tight_tol << std::endl;
            std::cout << "\t strong first iteration: " << this->strongFirstIteration << std::endl;
            std::cout << "\t strong recurring iteration: " << this->strongRecurring << std::endl;
            if(this->strongRecurring){
                std::cout << "\t \tinterval: " << this->recurringInterval << std::endl;
            }
        }
        //this->m_timeDisc;
        //this->m_out;

    }


    void setForceConv(bool pForceConvergence) {
        this->forceConv = pForceConvergence;
    }


    double tight_tol = 1e-9;
    void setTightTol(double pTightTol){
        this->tight_tol = pTightTol;
    }
    double loose_tol = 1e-2;

    void setLooseTol(double pLooseTol){
        this->loose_tol = pLooseTol;
    }
    int lastiter[15] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};


    braid_Int GetSpatialAccuracy(BraidStepStatus &pstatus,
                                 int level,
                                 double loose_tol,
                                 double tight_tol,
                                 double *tol_ptr) {
        int nrequest = 2;
        braid_Real stol, tol, rnorm, rnorm0, old_fine_tolx;
        braid_Real l_rnorm, l_ltol, l_ttol, l_tol;
        double *rnorms = (braid_Real *) malloc(2 * sizeof(braid_Real));
        pstatus.StepStatusGetTol(&tol);
        pstatus.SetOldFineTolx(old_fine_tolx);

        /* Get the first and then the current residual norms */
        rnorms[0] = -1.0;
        rnorms[1] = -1.0;
        pstatus.GetRNorms(&nrequest, rnorms);

        if ((rnorms[0] == -1.0) && (rnorms[1] != -1.0)) {
            rnorm0 = rnorms[1];
        } else {
            rnorm0 = rnorms[0];
        }
        nrequest = -2;
        pstatus.GetRNorms(&nrequest, rnorms);
        if ((rnorms[1] == -1.0) && (rnorms[0] != -1.0)) {
            rnorm = rnorms[0];
        } else {
            rnorm = rnorms[1];
        }


        if ((level > 0) || (nrequest == 0) || (rnorm0 == -1.0)) {
            /* Always return the loose tolerance, if
             * (1) On a coarse grid computation
             * (2) There is no residual history yet (this is the first Braid iteration with skip turned on) */
            *tol_ptr = loose_tol;
        } else {
            /* Else, do a variable tolerance for the fine grid */
            l_rnorm = -log10(rnorm / rnorm0);
            l_tol = -log10(tol / rnorm0);
            l_ltol = -log10(loose_tol);
            l_ttol = -log10(tight_tol);

            if (l_rnorm >= (7.0 / 8.0) * l_tol) {
                /* Close to convergence, return tight_tol */
                *tol_ptr = tight_tol;
            } else {
                /* Do linear interpolation between loose_tol and tight_tol (but with respect to log10) */
                stol = (l_rnorm / l_tol) * (l_ttol - l_ltol) + l_ltol;
                *tol_ptr = pow(10, -stol);

                /* The fine grid tolerance MUST never decrease */
                if (((*tol_ptr) > old_fine_tolx) && (old_fine_tolx > 0)) {
                    *tol_ptr = old_fine_tolx;
                }
            }
        }

        if (level == 0) {
            /* Store this fine grid tolerance */
            pstatus.SetOldFineTolx(*tol_ptr);

            /* If we've reached the "tight tolerance", then indicate to Braid that we can halt */
            if (*tol_ptr == tight_tol) {
                pstatus.SetTightFineTolx(1);
            } else {
                pstatus.SetTightFineTolx(0);
            }
        }

        free(rnorms);
        /* printf( "lev: %d, accuracy: %1.2e, nreq: %d, rnorm: %1.2e, rnorm0: %1.2e, loose: %1.2e, tight: %1.2e, old: %1.2e, braid_tol: %1.2e \n", level, *tol_ptr, nrequest, rnorm, rnorm0, loose_tol, tight_tol, old_fine_tolx, tol); */
        return _braid_error_flag;
    }



    /* @brief Apply the time stepping routine to the input vector @a u_
    corresponding to time @a tstart, and return in the same vector @a u_ the
    computed result for time @a tstop. The values of @a tstart and @a tstop
    can be obtained from @a pstatus.

    @param[in,out] u_ Input: approximate solution at time @a tstart.
                      Output: computed solution at time @a tstop.
    @param[in] ustop_ Previous approximate solution at @a tstop?
    @param[in] fstop_ Additional source at time @a tstop. May be set to NULL,
                      indicating no additional source.

    @see braid_PtFcnStep.
    virtual braid_Int Step(braid_Vector     u_,
                           braid_Vector     ustop_,
                           braid_Vector     fstop_,
                           BraidStepStatus &pstatus) = 0;*/

    /*
     * Defines the central time stepping function that the user must write.
     *
     * The user must advance the vector *u* from time *tstart* to *tstop*.  The time
     * step is taken assuming the right-hand-side vector *fstop* at time *tstop*.
     * The vector *ustop* may be the same vector as *u* (in the case where not all
     * unknowns are stored).  The vector *fstop* is set to NULL to indicate a zero
     * right-hand-side.
     *
     * Query the status structure with *braid_StepStatusGetTstart(status, &tstart)*
     * and *braid_StepStatusGetTstop(status, &tstop)* to get *tstart* and *tstop*.
     * The status structure also allows for steering.  For example,
     * *braid_StepStatusSetRFactor(...)* allows for setting a refinement factor,
     * which tells XBraid to refine this time interval.
     *
    typedef braid_Int
    (*braid_PtFcnStep)(braid_App        app,    / **< user-defined _braid_App structure /
                       braid_Vector     ustop,  / **< input, u vector at *tstop* /
                       braid_Vector     fstop,  / **< input, right-hand-side at *tstop* /
                       braid_Vector     u     , / **< input/output, initially u vector at *tstart*, upon exit, u vector at *tstop* /
                       braid_StepStatus status  / **< query this struct for info about u (e.g., tstart and tstop), allows for steering (e.g., set rfactor) /
    );*/


    /*
     *
     * @param u [in/out] input approximate solution at t=tstart, output computed solution
     * @param ustop previous approximate solution at t=tstop
     * @param fstop additional source at t=tstop or nullptr if not aviable
     * @param pstatus
     * @return
     */
    braid_Int Step(braid_Vector u, //
                   braid_Vector ustop,
                   braid_Vector fstop,
                   BraidStepStatus &pstatus) override {
#if TRACE_CONST == 1
        u->m_const = false;
        //ustop->m_const = false;
#endif
        int l; // level
        pstatus.GetLevel(&l);
        StartRedoranLevel(LevelObserver::T_STEP,l);
        double t_start, t_stop;
        pstatus.GetTstartTstop(&t_start, &t_stop);

        double current_dt = t_stop - t_start;

        // std::cout << "num stages " << this->m_timeDisc->num_stages() << std::endl;

        //this->debugwriter << "adaptConv" << std::endl;
        if (this->adaptConv) {
            StartRedoranLevel(LevelObserver::T_ADAPTIVE_TOL, l);
            int iteration;
            pstatus.GetIter(&iteration);

            int iter;
            double tol;
            if (iteration != this->lastiter[l]) {
                this->lastiter[l] = iteration;
                if (this->m_timing) {
                    double diff, total;
                    base_type::timer.now(total, diff);
                    base_type::debugwriter << std::setw(10) << "@time:"
                            << std::setw(12) << total << " ; "
                            << std::setw(12) << diff << " Begin iteration for level" << l << std::endl;
                }
                //if (this->m_verbose) {
                //    this->debugwriter << "========== ==========<< " << l << " >>========== ==========" << std::endl;
                //}
                if(this->strongFirstIteration && iteration == 0) {
                    tol = tight_tol;
                    iter = 10;
                    this->debugwriter << "strong first cycle    ";
                } else if (this->strongRecurring && iteration % this->recurringInterval == 0){
                    tol = tight_tol;
                    iter = 10;
                    this->debugwriter << "strong recurring    ";
                } else if (l == 0) {
                        this->GetSpatialAccuracy(pstatus, l, loose_tol, tight_tol, &tol);
                        iter = 100;
                    this->debugwriter << "level 0     ";
                } else {
                    tol = this->loose_tol;
                    iter = 2;
                    this->debugwriter << "level != 0    ";
                }


                SPConv conv = SPConv(new TConv(iter, 1e-32, tol, false));
                this->debugwriter << "% new Tolerance \t " << tol << std::endl;
                this->m_linSolver->set_convergence_check(conv);
            }
            StopRedoranLevel(LevelObserver::T_ADAPTIVE_TOL, l);
        }

        //this->debugwriter << "message" << std::endl;
        if (this->m_verbose) {
            int tindex;
            pstatus.GetTIndex(&tindex);
            int iteration;
            pstatus.GetIter(&iteration);
#if TRACE_INDEX == 1


            if (fstop == nullptr) {
                this->debugwriter << "u_" << u->index << " = step_"<<l<<"_n( u_" << u->index << ", u_" << ustop->index
                                  << ", null, " << t_start << ", " << current_dt << ", " << t_stop << ", " << l
                                  << ")"
                                  << "\t\t % " << tindex << std::endl;
            } else {
                this->debugwriter << "u_" << u->index << " = step_"<<l<<"_r( u_" << u->index << ", u_" << ustop->index
                                  << ", u_"
                                  << fstop->index << ", " << t_start << ", " << current_dt << ", " << t_stop << ", "
                                  << l << ")"
                                  << " % " << tindex << std::endl;
            }

#else
        this->debugwriter << std::setw(13) << iteration << "step for level " << l << " at position " << tindex << " and iteration" << std::endl;
#endif
        }

        //this->debugwriter << "preparation" << std::endl;
        auto *sp_u_approx_tstart = (SPGridFunction *) u->value;
        auto *constsp_u_approx_tstop = (SPGridFunction *) ustop->value;
        SPGridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();
        SPGridFunction lp = constsp_u_approx_tstop->get()->clone();

        SPGridFunction sp_rhs = this->m_u0->clone_without_values(); // for rhs
        const ug::GridLevel gridlevel = sp_u_approx_tstart->get()->grid_level();


        series->push(*sp_u_approx_tstart, t_start);
        m_timeDisc->prepare_step(series, current_dt);
        //this->debugwriter << "assemblation" << std::endl;
        if (fabs(current_dt - this->m_assembled_dt) > 1e-14) {
            StartRedoranLevel(LevelObserver::T_ASSEMBLE_OP,l);
            if(this->m_verbose){
                this->debugwriter << "% Assemble operator " << current_dt<< std::endl;
            }
            m_timeDisc->assemble_linear(*m_A, *sp_rhs.get(), gridlevel);
            m_linSolver->init(m_A, *lp.get());
            this->m_assembled_dt = current_dt;
            StopRedoranLevel(LevelObserver::T_ASSEMBLE_OP,l);
        } else {
            StartRedoranLevel(LevelObserver::T_ASSEMBLE_RHS,l);
            m_timeDisc->assemble_rhs(*sp_rhs.get(), gridlevel);
            StopRedoranLevel(LevelObserver::T_ASSEMBLE_RHS,l);
        }


        //this->debugwriter << "fstop" << std::endl;
        if (fstop != nullptr) {
            auto *sp_fstop = (SPGridFunction *) fstop->value;
            //VecAdd(a,x,b,y) <=> x = a * x + b * y
            VecAdd(1, *sp_rhs.get(), 1, *sp_fstop->get());
        }
#if TRACE_INDEX == 1
        MATLAB(sp_rhs->clone().get(), u->index, t_stop);
#endif
        //this->debugwriter << "solve" << std::endl;
        StartRedoranLevel(LevelObserver::T_SOLVE,l);
        bool success = m_linSolver->apply(*sp_u_tstop_approx.get(), *sp_rhs.get());
        StopRedoranLevel(LevelObserver::T_SOLVE,l);

#if TRACE_DEFECT == 1
        if(this->m_verbose){
        this->debugwriter << std::setw(20) << "@conv Iterations: " << std::setw(12) << m_linSolver->step()
                          << std::setw(20) << "Reduction: " << std::setw(12) << m_linSolver->reduction()
                          << std::setw(20) << "Defect: " << std::setw(12) << m_linSolver->defect() << std::endl;
    }
#endif

        if (!success && forceConv) {
            this->debugwriter << "!!! Failure convergence not reached" << std::endl;
            exit(127);
        }
        //this->debugwriter << "output" << std::endl;
#if TRACE_INDEX == 1
        MATLAB(sp_u_tstop_approx.get(), u->index, t_stop);
#endif
        *sp_u_approx_tstart = sp_u_tstop_approx;
        series->clear();
        //this->debugwriter << "end" << std::endl;
        StopRedoranLevel(LevelObserver::T_STEP,l);
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
        m_timeDisc->prepare_step(series, current_dt);
        auto sp_rhs = this->m_u0->clone_without_values();
        if (fabs(current_dt - this->m_assembled_dt) > 1e-14) {
            StartRedoranLevel(LevelObserver::T_ASSEMBLE_OP,l);
            m_timeDisc->assemble_linear(*m_A, *sp_rhs.get(), gridlevel);

            m_linSolver->init(m_A, *const_u_approx_tstop->get());

            this->m_assembled_dt = current_dt;
            StopRedoranLevel(LevelObserver::T_ASSEMBLE_OP,l);
        } else {
            StartRedoranLevel(LevelObserver::T_ASSEMBLE_RHS,l);
            m_timeDisc->assemble_rhs(*sp_rhs.get(), gridlevel);
            StopRedoranLevel(LevelObserver::T_ASSEMBLE_RHS,l);
        }


        m_linSolver->linear_operator()->apply_sub(
                *sp_rhs.get(), // f co domain function [in / out]
                *const_u_approx_tstop->get() // u domain function [in]
        ); // calculates r = r - A * u
        // r = rhs  - L*u_stop

        (*sp_rhs) *= -1; // r = -r + A*u
         *u_approx_tstart = sp_rhs;

         //SPGridFunction  def = sp_rhs->clone(); // todo alternative?
        //m_timeDisc->assemble_defect(*def.get(),*const_u_approx_tstop->get(),gridlevel);
        //VecAdd(1,*def.get(),-1,*u_approx_tstart->get());
#if TRACE_INDEX == 1
        MATLAB(sp_rhs.get(), u->index, t_stop);
#endif
        //*u_approx_tstart = *const_u_approx_tstop;
        series->clear();
        StopRedoranLevel(LevelObserver::T_RESIDUAL,l);
        return 0;
    };


    braid_Int SpatialNorm(braid_Vector u, braid_Real *norm_ptr) override {
        auto *uref = (SPGridFunction *) u->value;
        // todo clone for non residual
	SPGridFunction tempobject = uref->get()->clone();
        *norm_ptr = tempobject->norm();
#if TRACE_INDEX == 1
        if (this->m_verbose) {
            this->debugwriter << "norm( u_" << u->index << ") % " << *norm_ptr << std::endl;
        }
#endif
        return 0;
    };


    braid_Int Coarsen(braid_Vector fu, braid_Vector *cu, BraidCoarsenRefStatus &status)
    override {
     /*	Clone(fu, cu);

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
    	/*Clone(cu, fu);

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

#endif //UG_PLUGIN_XBRAIDFORUG4_RGFBRAIDAPP_H
