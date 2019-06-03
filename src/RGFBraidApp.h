//
// Created by parnet on 24.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_RGFBRAIDAPP_H
#define UG_PLUGIN_XBRAIDFORUG4_RGFBRAIDAPP_H

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
class RGFBraidApp : public GFBraidApp<TDomain, TAlgebra> {

public:

    /****************************************************************************
    * Typedefs
    ***************************************************************************/
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

    /****************************************************************************
    * Members
    ***************************************************************************/
    SPTimeStep timeDisc;
    SPSolver linSolver;
    SPTimeIntegrator m_spIntegratorC;
    SPTimeIntegrator m_spIntegratorF;
    SPOutput out;

    SPGridFunction ux; // for t > tstart
    SmartPtr<XCommunicator> comm;

    const char *filename{};

    /****************************************************************************
    * Functions to outsource // todo
    ***************************************************************************/
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

public:

    RGFBraidApp(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) :
            GFBraidApp<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {}

    ~RGFBraidApp() override = default;


    void setStartVector(SPGridFunction p_u0) {
        this->u0 = p_u0;
    }

    void setRemainingVector(SPGridFunction p_ux) {
        this->ux = p_ux;
    }

    void setTimeDisc(SPTimeStep p_timeDisc, size_t level = 0) {
        this->timeDisc = p_timeDisc;
    }

    void setLinSolver(SPSolver p_linSolver, size_t level = 0) {
        this->linSolver = p_linSolver;
    }


    double assembled_dt = 0;
    SmartPtr<ug::AssembledLinearOperator<TAlgebra>> A;
    //SPGridFunction defect;

    void init() {
        const ug::GridLevel gridlevel = this->u0->grid_level();
        A = SmartPtr<ug::AssembledLinearOperator<TAlgebra>>(
                new ug::AssembledLinearOperator<TAlgebra>(timeDisc, gridlevel));
    }

    //todo change for user function
    braid_Int Init(braid_Real t, braid_Vector *u_ptr) override {
        std::cout << "Creating new Vector @t=" << t << "\t" << indexpool << std::endl;
        auto *u = (BraidVector *) malloc(sizeof(BraidVector));
        SPGridFunction *vec;
        if (t == this->tstart) {
            vec = new SPGridFunction(new TGridFunction(*this->u0));
        } else {
            vec = new SPGridFunction(new TGridFunction(*this->ux));
        }
        u->value = vec;
        u->time = t;
        u->index = indexpool;
        indexpool++;
        *u_ptr = u;
        return 0;
    };

    int progress = 0;
    bool writeparam = true;
    bool dispstorageinfo = false;

    /** @brief Apply the time stepping routine to the input vector @a u_
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

    /**
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


    /**
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
        double t_start, t_stop;
        pstatus.GetTstartTstop(&t_start, &t_stop);

        double current_dt = t_stop - t_start;
        int l; // level
        pstatus.GetLevel(&l);

        if (this->verbose) {
            std::cout << progress << "\t" << l << "\tStep with dt=" << current_dt << "\t t = [" << t_start << ", "
                      << t_stop << "] ";

            if (fstop != nullptr) {
                std::cout << "\tR - set";
            }
            std::cout << std::endl;
            std::cout << "u(" << u->time << "): " << u->index << "\tustop(" << ustop->time << "): " << ustop->index;
            if (fstop != nullptr) {
                std::cout << "\tfstop(" << fstop->time << "): " << fstop->index;
            }
            std::cout << std::endl;
        }

        auto *sp_u_tstart = (SPGridFunction *) u->value;
        auto *constsp_u_tstop_approx = (SPGridFunction *) ustop->value;
        SPGridFunction sp_u_tstop_approx = constsp_u_tstop_approx->get()->clone();

        SPGridFunction sp_rhs = this->u0->clone_without_values(); // for rhs
        const ug::GridLevel gridlevel = sp_u_tstart->get()->grid_level();

        if (true) {//this->writeparam) {
            write(sp_u_tstart->get(), progress, t_start, "Sustart");
            write(constsp_u_tstop_approx->get(), progress, t_start, "Suapprox");

            if (fstop != nullptr) {
                auto *sp_rhs_stop = (SPGridFunction *) fstop->value;
                write(sp_rhs_stop->get(), progress, t_start, "Sfstop");
            }
        }

        SPTimeSeries series = SPTimeSeries(new TTimeSeries());
        series->push(sp_u_tstart->get()->clone(), t_start);
        timeDisc->prepare_step(series, current_dt);

        if (current_dt != assembled_dt) { // todo is close?
            timeDisc->assemble_linear(*A, *sp_rhs.get(), gridlevel);
            linSolver->init(A, *sp_u_tstop_approx.get());
            assembled_dt = current_dt;
        } else {
            timeDisc->assemble_rhs(*sp_rhs.get(), gridlevel); // not neccassary
        }

        bool success;
        if (fstop != nullptr) {
            auto *sp_rhs_stop = (SPGridFunction *) fstop->value;
            VecAdd(1, *sp_rhs.get(), 1, *sp_rhs_stop->get());

        }
        success = linSolver->apply(*sp_u_tstop_approx.get(), *sp_rhs.get());


        if (!success) {
            std::cout << "Failure" << std::endl;
            exit(127);
        }

        *sp_u_tstart = sp_u_tstop_approx;
        u->time = t_stop;


        write(sp_u_tstart->get(), progress, t_stop, "Sresult");
        exit(0);
        progress++; // todo delete
        return 0;
    };



    /** @brief Compute the residual at time @a tstop, given the approximate
    solutions at @a tstart and @a tstop. The values of @a tstart and @a tstop
    can be obtained from @a pstatus.

    @param[in]     u_ Input: approximate solution at time @a tstop.
    @param[in,out] r_ Input: approximate solution at time @a tstart.
                      Output: residual at time @a tstop.

    @see braid_PtFcnResidual.

    virtual braid_Int Residual(braid_Vector     u_,
                               braid_Vector     r_,
                               BraidStepStatus &pstatus) = 0;*/
    /**
     * This function (optional) computes the residual *r* at time *tstop*.  On
     * input, *r* holds the value of *u* at *tstart*, and *ustop* is the value of
     * *u* at *tstop*.  If used, set with @ref braid_SetResidual.
     *
     * Query the status structure with *braid_StepStatusGetTstart(status, &tstart)*
     * and *braid_StepStatusGetTstop(status, &tstop)* to get *tstart* and *tstop*.
     *
    typedef braid_Int
    (*braid_PtFcnResidual)(braid_App        app,    / **< user-defined _braid_App structure /
                           braid_Vector     ustop,  / **< input, u vector at *tstop* /
                           braid_Vector     r     , / **< output, residual at *tstop* (at input, equals *u* at *tstart*) /
                           braid_StepStatus status  / **< query this struct for info about u (e.g., tstart and tstop) /
    );*/


    /**
     * Calculates the residual value for given approximate solution at t=tstart and t=tstop by calculating
     *
     * r = - r + A * u  - forcing - boundary
     * r = - ustart - forcing - boundary + A * ustop
     *
     * @param u [in] approximate solution at t=tstop
     * @param r [in / out] input is approximate solution at t=tstart and output residual at t=tstop
     * @param pstatus
     * @return
     */
    braid_Int Residual(braid_Vector u, braid_Vector r,
                       BraidStepStatus &pstatus) override {
        double t_start, t_stop;
        pstatus.GetTstartTstop(&t_start, &t_stop);
        double current_dt = t_stop - t_start;

        int l; // level;
        pstatus.GetLevel(&l);

        if (this->verbose) {
            std::cout << progress << "\t" << l << "\tResidual with dt=" << current_dt << "\t t = [" << t_start << ", "
                      << t_stop << "]" << std::endl;
            std::cout << "u(" << u->time << "): " << u->index << "\tr(" << r->time << "): " << r->index << std::endl;

        }

        auto *constsp_u_start = (SPGridFunction *) u->value; // input
        auto *sp_u_stop = (SPGridFunction *) r->value; // input
        //auto *sp_r_stop = (SPGridFunction *) r->value; // output == u_stop!

        if (this->writeparam) {
            write(constsp_u_start->get(), progress, t_start, "Rustart");
            write(sp_u_stop->get(), progress, t_start, "Rustop");
        }


        const ug::GridLevel gridlevel = constsp_u_start->get()->grid_level();

        SmartPtr<TTimeSeries> series = SmartPtr<TTimeSeries>(new TTimeSeries());

        series->push(*constsp_u_start, t_start); //todo clone sp_current_u_ref ?
        // r = utstart
        // addboundary(r)
        // addforce(tstop, r) [ r+= dt*F(x,y)
        timeDisc->prepare_step(series, current_dt);
        auto sp_rhs = this->u0->clone_without_values();
        if (current_dt != assembled_dt) { // todo check for is close?
            timeDisc->assemble_linear(*A, *sp_rhs.get(), gridlevel);
            linSolver->init(A, *sp_u_stop->get()); // todo linearization point ?
            assembled_dt = current_dt;
        } else {
            timeDisc->assemble_rhs(*sp_rhs.get(), gridlevel);
        }

        /** todo delete docu of ug4
         * This method applies the operator and subracts the result from the input
         * codomain function, i.e. f -= L*u (or d -= J(u)*c in iterative schemes).
         * Note, that the operator must have been initialized once before this
         * method can be used.
         *
         * \param[in]		u		domain function
         * \param[in,out]	f		codomain function
         * \returns			bool	success flag
         * 	virtual void apply_sub(Y& f, const X& u) = 0;
         */
        linSolver->linear_operator()->apply_sub(
                *sp_rhs.get(), // f co domain function [in / out]
                *sp_u_stop->get() // u domain function [in]
        ); // calculates r = r - A * u

        (*sp_rhs) *= -1; // r = -r + A*u
        *sp_u_stop = sp_rhs; // todo avoid copy?

        r->time = t_stop;
        if (this->writeparam) {
            write(sp_rhs.get(), progress, t_start, "Rresult");
            //exit(-1);
        }
        //exit(77);
        std::cout << "norm(r) = \t" << sp_rhs.get()->norm() << std::endl;

        progress++; // todo delete
        return 0;

    };


    braid_Int SpatialNorm(braid_Vector u, braid_Real *norm_ptr) override {
        auto *uref = (SPGridFunction *) u->value;
        *norm_ptr = (*uref)->norm(); // todo replace?
        return 0;
    };

    braid_Int Access(braid_Vector u, BraidAccessStatus &astatus) override {
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
    braid_Int Coarsen(braid_Vector fu, braid_Vector *cu, BraidCoarsenRefStatus &status) override {
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
    braid_Int Refine(braid_Vector cu, braid_Vector *fu, BraidCoarsenRefStatus &status) override {
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


};

#endif //UG_PLUGIN_XBRAIDFORUG4_RGFBRAIDAPP_H
