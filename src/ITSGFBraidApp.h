//
// Created by parnet on 29.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_ITSGFBRAIDAPP_H
#define UG_PLUGIN_XBRAIDFORUG4_ITSGFBRAIDAPP_H


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
class ITSGFBraidApp : public GFBraidApp<TDomain, TAlgebra> {

public:

    /****************************************************************************
    * Typedefs
    ***************************************************************************/
    typedef GFBraidApp<TDomain, TAlgebra> base_type;
    typedef MGFBraidApp<TDomain, TAlgebra> this_type;

    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;

    typedef SmartPtr<TGridFunction> SPGridFunction;

    typedef ug::ITimeDiscretization <TAlgebra> TTimeStep;

    typedef SmartPtr<TTimeStep> SPTimeStep;

    typedef ug::VectorTimeSeries<typename TAlgebra::vector_type> TTimeSeries;

    typedef SmartPtr<ug::VectorTimeSeries<typename TAlgebra::vector_type> > SPTimeSeries;

    typedef ug::ILinearOperatorInverse<typename TAlgebra::vector_type> TLinearSolver;

    typedef SmartPtr <TLinearSolver> SPLinearSolver;


    //typedef ug::IOperatorInverse<typename TAlgebra::vector_type> TSolver;
    typedef ug::NewtonSolver<TAlgebra> TSolver;

    typedef SmartPtr <TSolver> SPSolver;


    typedef ug::AssembledLinearOperator <TAlgebra> TAssembledLinearOperator;

    typedef SmartPtr <TAssembledLinearOperator> SPAssembledLinearOperator;

    typedef ug::AssembledOperator <TAlgebra> TAssembledOperator;

    typedef SmartPtr <TAssembledOperator> SPAssembledOperator;


    typedef SmartPtr<ug::ITimeIntegrator<TDomain, TAlgebra> > SPTimeIntegrator;

    typedef SmartPtr<ug::VTKOutput<TDomain::dim> > SPOutput;


    struct LevelSolver {
        /*
         * Linear only Theta, ImplEuler, ExplEuler, CrankNicholson and BDF
         * NonLinear SDIRK, Alexander, FracStep,
         */
        bool linear = false;
        bool forceConv = false;

        short storeOperator = 0;
        double assembled_dt = -1;

        SPAssembledLinearOperator A;
        SPAssembledOperator NLA;

        SPTimeStep timeDisc;

        SPLinearSolver linearSolver;
        SPSolver solver;
    };


    double shared_assembled_dt = -1;
    SPAssembledLinearOperator shared_A;
    SPAssembledOperator shared_NLA;


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
    ITSGFBraidApp() : GFBraidApp<TDomain, TAlgebra>()
	{
        this->name = "Level dependent Higher Order Time Discretization";
    }

    ~ITSGFBraidApp()
    {
        delete[] this->m_lvl;
    };

    void setMaxLevels(size_t levelcount) {
        delete[] this->m_lvl;
        this->m_levels = levelcount;
        this->m_lvl = new LevelSolver[this->m_levels];
    }


    /**
     * Assembled operator of coarse level ( level > x ) will not be stored.
     * To save all assembled operators use x == MaxLevel
     * @param level

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
    }*/

    bool m_residual = true;

    void setResidual(bool p_residual, size_t level = -1){
        m_residual = p_residual;
    }
    void setStoreOperator(bool p_value, size_t level = 0) {
        this->m_lvl[level].storeOperator = p_value;
    }

    void setTimeDisc(SPTimeStep p_timeDisc, size_t level = 0) {
        //std::cout << "set timeDisc for level " << level << std::endl;
        this->m_lvl[level].timeDisc = p_timeDisc;
    }

    void setSolver(SPSolver p_solver, size_t level = 0) {
        // std::cout << "set linearSolver for level " << level << std::endl;
        this->m_lvl[level].linear = false;
        this->m_lvl[level].solver = p_solver;

        //std::cout << "solver set for level " << level << std::endl;
    }

    void setLinearSolver(SPLinearSolver p_solver, size_t level = 0) {
        // std::cout << "set linearSolver for level " << level << std::endl;
        this->m_lvl[level].linear = true;
        this->m_lvl[level].linearSolver = p_solver;
        //std::cout << "linear solver set for level " << level << std::endl;
    }

    void setForceConv(bool pForceConvergence, size_t level = 0) {
        // std::cout << "set force conv for level " << level << "\t" << pForceConvergence << std::endl;
        this->m_lvl[level].forceConv = pForceConvergence;
    }

    void init() override {
    	std::stringstream ss;
    	ss << "job_" << this->m_comm->getTemporalRank() << ".output";
    	this->debugwriter.open(ss.str());

    	base_type::timer.start();
#if TRACE_TIMINGS == 1
        this->redoran = Redoran(this->m_levels);
#endif

#if TRACE_GRIDFUNCTION == 1

        this->matlab = SmartPtr<MATLABScriptor<TDomain, TAlgebra> >(new MATLABScriptor<TDomain, TAlgebra>(this->debugwriter));
#endif
        const ug::GridLevel gridlevel = this->m_u0->grid_level();

        this->shared_A = SPAssembledLinearOperator(new TAssembledLinearOperator());
        this->shared_NLA = SPAssembledOperator(new TAssembledOperator());
        this->shared_assembled_dt = -1;

        for (size_t i = 0; i < this->m_levels; ++i) {
            if (this->m_lvl[i].storeOperator) {
                if (this->m_lvl[i].linear) {
                    this->m_lvl[i].A = SPAssembledLinearOperator(
                            new TAssembledLinearOperator(this->m_lvl[i].timeDisc, gridlevel));
                } else {
                    //this->debugwriter << "create non linear operator for level " << i << std::endl;
                    this->m_lvl[i].NLA = SPAssembledOperator( new TAssembledOperator(this->m_lvl[i].timeDisc, gridlevel));
                }
            }
        }
    }


    braid_Int Step(braid_Vector u, //
                   braid_Vector ustop,
                   braid_Vector fstop,
                   BraidStepStatus &pstatus) override {
        int l; // level
        pstatus.GetLevel(&l);
        StartRedoranLevel(LevelObserver::T_STEP, l);


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
                debugwriter << "u_" << u->index << " = step_" << l << "_n( u_" << u->index << ", u_"
                                  << ustop->index
                                  << ", null, " << t_start << ", " << current_dt << ", " << t_stop << ", " << l
                                  << ")"
                                  << "\t\t % " << tindex << std::endl;
            } else {
                debugwriter << "u_" << u->index << " = step_" << l << "_r( u_" << u->index << ", u_"
                                  << ustop->index
                                  << ", u_"
                                  << fstop->index << ", " << t_start << ", " << current_dt << ", " << t_stop << ", "
                                  << l << ")"
                                  << " % " << tindex << std::endl;
            }
#else
            this->debugwriter << std::setw(13) << iteration << "step for level " << l << " at position " << tindex
                              << " and iteration" << std::endl;
#endif
        }

        auto *sp_u_approx_tstart = (SPGridFunction *) u->value;
        auto *constsp_u_approx_tstop = (SPGridFunction *) ustop->value;
        SPGridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();
        SPGridFunction lp = constsp_u_approx_tstop->get()->clone();

        SPGridFunction sp_rhs = this->m_u0->clone_without_values(); // for rhs
        const ug::GridLevel gridlevel = sp_u_approx_tstart->get()->grid_level();

        //SPGridFunction last_time_step = sp_u_approx_tstart->get()->clone(); // todo delete?
        bool success = false;
        if (this->m_lvl[l].linear) {
            //this->debugwriter << "linear case" << std::endl;

            series->push(*sp_u_approx_tstart, t_start);
            this->m_lvl[l].timeDisc->prepare_step(series, current_dt);
            bool reassemble = false;
            SPAssembledLinearOperator *spA;
            if (this->m_lvl[l].storeOperator) {
                //this->debugwriter << "store operator case" << std::endl;
                spA = &this->m_lvl[l].A;
                if (this->m_lvl[l].assembled_dt == -1.0) {
                    //this->debugwriter << "create operator for level " << l << std::endl;
                    this->m_lvl[l].assembled_dt = current_dt;
                    reassemble = true;
                }
            } else {
                //this->debugwriter << "shared operator case " << std::endl;
                spA = &this->shared_A;
                if (fabs(this->shared_assembled_dt - current_dt) > 1e-14) {
                    (*spA)->set_discretization(this->m_lvl[l].timeDisc);
                    (*spA)->set_level(gridlevel);
                    //this->debugwriter << "create shared operator" << std::endl;
                    this->shared_assembled_dt = current_dt;
                    reassemble = true;
                }
            }

            if (reassemble) {
                //this->debugwriter << "assemble operator" << std::endl;
                this->m_lvl[l].timeDisc->assemble_linear(*spA->get(), *sp_rhs.get(), gridlevel);
                //this->debugwriter << spA->valid() << std::endl;
                //this->debugwriter << sp_u_tstop_approx.valid() << std::endl;
                //this->debugwriter << this->m_lvl[l].linearSolver.valid() << std::endl;
                this->m_lvl[l].linearSolver->init(*spA, *sp_u_tstop_approx.get());


            } else {
                //this->debugwriter << "assemble vector" << std::endl;
                this->m_lvl[l].timeDisc->assemble_rhs(*sp_rhs.get(), gridlevel);
            }
            //this->debugwriter << "use fstop" << std::endl;
            if (fstop != nullptr) {
                auto *sp_fstop = (SPGridFunction *) fstop->value;
                //VecAdd(a,x,b,y) <=> x = a * x + b * y
                VecAdd(1, *sp_rhs.get(), 1, *sp_fstop->get());
            }

            //this->debugwriter << "solve system" << std::endl;
            //this->debugwriter << "apply operator" << std::endl;
            bool success = this->m_lvl[l].linearSolver->apply(*sp_u_tstop_approx.get(), *sp_rhs.get());

            if (!success) {
                if (this->m_lvl[l].forceConv) {
                    //this->debugwriter << "Failure Convergence not reached" << std::endl;
                    exit(127);
                }
            }
            //this->debugwriter << "finished linear part" << std::endl;

        } else {
            //this->debugwriter << "non linear case" << std::endl;
            series->push(*sp_u_approx_tstart, t_start);
            //this->debugwriter << "\t after push" << std::endl;

            //prepare_step(u, step, time, time, current_dt);
            SmartPtr<TAssembledOperator> *NLA;
            if(this->m_lvl[l].storeOperator){
                //this->debugwriter << "\t use operator for level" << std::endl;
                NLA = &this->m_lvl[l].NLA;
            } else {
                //this->debugwriter << "\t use shared operator" << std::endl;
                NLA = &this->shared_NLA;
                NLA->get()->set_discretization(this->m_lvl[l].timeDisc);
                NLA->get()->set_level(gridlevel);
            }
            //this->debugwriter << "\t operator assigned" << std::endl;
            size_t stages = this->m_lvl[l].timeDisc->num_stages();
            for (size_t i_stage = 1; i_stage  <= stages; ++i_stage) {
                //this->debugwriter << "Stage - " << i_stage << std::endl;
                //preprocess(u, step, time, currdt);
                this->m_lvl[l].timeDisc->set_stage(i_stage);
                //this->debugwriter << "\t stage set" << std::endl;
                this->m_lvl[l].timeDisc->prepare_step(series, current_dt);
                //this->debugwriter << "\t step prepared" << std::endl;
                // todo newtonSolver:set_line_search(defaultLineSearch)
                this->m_lvl[l].solver->init(*NLA);
                //this->debugwriter << "\t operator init" << std::endl;
                this->m_lvl[l].solver->prepare(*sp_u_tstop_approx.get());
                //this->debugwriter << "\t solver prepared" << std::endl;
                //this->debugwriter << "=== config string" << std::endl;
                //this->debugwriter << "\t config splinsolver" << std::endl;
                //this->debugwriter << this->m_lvl[l].solver->m_spLinearSolver.valid() << std::endl;
                //this->debugwriter << this->m_lvl[l].solver->m_spLinearSolver->config_string() << std::endl;
                //this->debugwriter << "\t config spconvchck" << std::endl;
                //this->debugwriter << this->m_lvl[l].solver->m_spConvCheck.valid() << std::endl;
                //this->debugwriter << this->m_lvl[l].solver->m_spConvCheck->config_string() << std::endl;
                //this->debugwriter << "\t config linesearch" << std::endl;
                //this->debugwriter << this->m_lvl[l].solver->m_spLineSearch.valid() << std::endl;
                //this->debugwriter << "\t config freq" << std::endl;
                //this->debugwriter << this->m_lvl[l].solver->m_reassembe_J_freq << std::endl;
                //this->debugwriter << "\t complete" << std::endl << std::flush;

                //this->debugwriter << this->m_lvl[l].solver->config_string() << std::endl;
                //this->debugwriter << "<<< config string" << std::endl;
                success = this->m_lvl[l].solver->apply(*sp_u_tstop_approx.get());
                //this->debugwriter << "\t solver applied" << std::endl;
                //this->debugwriter << "\t success? : \t " << success << std::endl;
                // postprocess(u,step,timeDisc.futureTime(),currdt);
            }
        }

#if TRACE_GRIDFUNCTION == 1
        MATLAB(sp_u_tstop_approx.get(), u->index, t_start);
#endif
        *sp_u_approx_tstart = sp_u_tstop_approx;
        series->clear();
        //this->debugwriter << "finished step" << std::endl;
        StopRedoranLevel(LevelObserver::T_STEP, l);

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
        StartRedoranLevel(LevelObserver::T_RESIDUAL, l);

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

        SPAssembledLinearOperator *spA;
        bool reassemble = false;
        if (this->m_lvl[l].storeOperator) {
            spA = &this->m_lvl[l].A;
            if (this->m_lvl[l].assembled_dt == -1.0) {
                this->debugwriter << "create operator for level Residual " << l << std::endl;
                this->m_lvl[l].assembled_dt = current_dt;
                reassemble = true;
            }
        } else {
            spA = &this->shared_A;
            if (fabs(this->shared_assembled_dt - current_dt) > 1e-14) {
                this->debugwriter << "create shared operator Residual" << std::endl;
                (*spA)->set_discretization(this->m_lvl[l].timeDisc);
                (*spA)->set_level(gridlevel);
                this->shared_assembled_dt = current_dt;
                reassemble = true;
            }
        }
        if (reassemble) {
            std::cout << "reassemble residual " << std::endl;
            StartRedoranLevel(LevelObserver::T_ASSEMBLE_OP, l);

            this->m_lvl[l].timeDisc->assemble_linear(*spA->get(), *sp_rhs.get(), gridlevel);
            this->m_lvl[l].linearSolver->init(*spA, *u_approx_tstart->get());
            StopRedoranLevel(LevelObserver::T_ASSEMBLE_OP, l);
        } else {
            StartRedoranLevel(LevelObserver::T_ASSEMBLE_RHS, l);
            this->m_lvl[l].timeDisc->assemble_rhs(*sp_rhs.get(), gridlevel);
            StopRedoranLevel(LevelObserver::T_ASSEMBLE_RHS, l);
        }


        this->m_lvl[l].linearSolver->linear_operator()->apply_sub(
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
        StopRedoranLevel(LevelObserver::T_RESIDUAL, l);
        return 0;

    };


    braid_Int SpatialNorm(braid_Vector u, braid_Real *norm_ptr) override {
        auto *uref = (SPGridFunction *) u->value;
        // todo clone for non residual
        if(m_residual){ // xbraid residual
            SPGridFunction copy = (*uref)->clone();
            *norm_ptr = copy->norm();
        } else {
            *norm_ptr = (*uref)->norm();
        }

#if TRACE_INDEX == 1
        if (this->m_verbose) {
            this->debugwriter << "norm( u_" << u->index << ") % " << *norm_ptr << std::endl;
        }
#endif
        return 0;
    };


    braid_Int Coarsen(braid_Vector fu, braid_Vector *cu, BraidCoarsenRefStatus &status)
    override {
    	/*Clone(fu, cu);

        auto *sp_fu = (SPGridFunction *) fu->value;
        auto *sp_cu = (SPGridFunction *) (*cu)->value;

        double t_upper;
        double t_lower;
        status.GetT(&t_lower);
        status.GetCTstop(&t_upper);
        //status.GetFTstop(&t_upper);
*/
        // m_spIntegratorC->apply(*sp_fu, t_upper, sp_cu->cast_const(), t_lower); // todo check fu, cu order
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

        // m_spIntegratorF->apply(*sp_cu, t_upper, sp_fu->cast_const(), t_lower); // todo check fu, cu order
*/
        return 0;
    }
};


#endif //UG_PLUGIN_XBRAIDFORUG4_ITSGFBRAIDAPP_H
