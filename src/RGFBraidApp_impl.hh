#ifndef _UG_PLUGIN_XBRAIDFORUG4_RGFBRAIDAPP_IMPL_HH_
#define _UG_PLUGIN_XBRAIDFORUG4_RGFBRAIDAPP_IMPL_HH_



    /*
     *
     * @param u [in/out] input approximate solution at t=tstart, output computed solution
     * @param ustop previous approximate solution at t=tstop
     * @param fstop additional source at t=tstop or nullptr if not aviable
     * @param pstatus
     * @return
     */
template<typename TDomain, typename TAlgebra>
braid_Int RGFBraidApp<TDomain,TAlgebra>::Step(braid_Vector u, //
                   braid_Vector ustop,
                   braid_Vector fstop,
                   BraidStepStatus &pstatus)
{
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
        // TGridFunction& ref = TBraidVectorFunctors().as_grid_function(*u);
        auto *sp_u_approx_tstart = (SPGridFunction *) u->value;
        auto *constsp_u_approx_tstop = (SPGridFunction *) ustop->value;
        SPGridFunction sp_u_tstop_approx = constsp_u_approx_tstop->get()->clone();

        SPGridFunction aux = constsp_u_approx_tstop->get()->clone();
       //  TGridFunction &aux_ = *aux;

        SPGridFunction sp_rhs = this->m_u0->clone_without_values(); // for rhs
        TGridFunction &rhs_ = *sp_rhs;

        const ug::GridLevel gridlevel = sp_u_approx_tstart->get()->grid_level();


        series->push(*sp_u_approx_tstart, t_start);
        m_timeDisc->prepare_step(series, current_dt);
        //this->debugwriter << "assemblation" << std::endl;
        if (fabs(current_dt - this->m_assembled_dt) > 1e-14) {
        	// Assemble J,b
            StartRedoranLevel(LevelObserver::T_ASSEMBLE_OP,l);
            if(this->m_verbose){
                this->debugwriter << "% Assemble operator " << current_dt<< std::endl;
            }
            m_timeDisc->assemble_linear(*m_A, rhs_, gridlevel);
            m_linSolver->init(m_A, *aux);
            this->m_assembled_dt = current_dt;
            StopRedoranLevel(LevelObserver::T_ASSEMBLE_OP,l);
        } else {
        	// Assemble b only
            StartRedoranLevel(LevelObserver::T_ASSEMBLE_RHS,l);
            m_timeDisc->assemble_rhs(rhs_, gridlevel);
            StopRedoranLevel(LevelObserver::T_ASSEMBLE_RHS,l);
        }

        // OPTIONAL: Add term to b = b + b'
        if (fstop != nullptr) {

            // auto *sp_fstop = (SPGridFunction *) fstop->value;
        	// usingTBraidVectorFunctors;
            TGridFunction& fstop_ =  typename base_type::TBraidVectorFunctors().as_grid_function(*fstop);
            //VecAdd(a,x,b,y) <=> x = a * x + b * y
            VecAdd(1.0, rhs_, 1.0, fstop_);
        }
#if TRACE_INDEX == 1
        MATLAB(sp_rhs->clone().get(), u->index, t_stop);
#endif

        // Solve Jx=b
        //this->debugwriter << "solve" << std::endl;
        StartRedoranLevel(LevelObserver::T_SOLVE,l);
        bool success = m_linSolver->apply(*sp_u_tstop_approx.get(), rhs_);
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

        // RETURN value
        *sp_u_approx_tstart = sp_u_tstop_approx;
        series->clear();
        //this->debugwriter << "end" << std::endl;
        StopRedoranLevel(LevelObserver::T_STEP,l);
        return 0;
    };

template<typename TDomain, typename TAlgebra>
braid_Int RGFBraidApp<TDomain,TAlgebra>::Residual(braid_Vector u, braid_Vector r,
                       BraidStepStatus &pstatus)
{
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
#endif
