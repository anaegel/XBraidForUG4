#ifndef _UG_PLUGIN_XBRAIDFORUG4_PBRAIDAPP_IMPL_HH_
#define _UG_PLUGIN_XBRAIDFORUG4_PBRAIDAPP_IMPL_HH_



//todo change for user function
template<typename TDomain, typename TAlgebra>
braid_Int PBraidApp<TDomain,TAlgebra>::Init(braid_Real t, braid_Vector *u_ptr)
{
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
template<typename TDomain, typename TAlgebra>
braid_Int PBraidApp<TDomain,TAlgebra>::Step(braid_Vector u, //
		braid_Vector ustop, // estimated solution?
		braid_Vector fstop,
		BraidStepStatus &pstatus)
{
	double tstart;
	double tstop;
	pstatus.GetTstartTstop(&tstart, &tstop);
	double current_dt = tstop - tstart;

	auto *sp_u_start = (SPGridFunction *) u->value;
	auto *sp_u_stop = (SPGridFunction *) u->value;
	auto *sp_u_stop_approx = (SPGridFunction *) ustop->value;

	auto *sp_rhs = new SPGridFunction(new TGridFunction(*this->u0));

	const ug::GridLevel gridlevel = sp_u_stop->get()->grid_level();

	SPTimeSeries series = SPTimeSeries(new TTimeSeries());
	series->push(sp_u_start->get()->clone(), tstart);
	timeDisc->prepare_step(series, current_dt);

	if (current_dt != assembled_dt) { // todo is close?
			// Reassemble operator.
			timeDisc->assemble_linear(*A, *sp_rhs->get(), gridlevel);
			linSolver->init(A, *sp_u_stop_approx->get());
			assembled_dt = current_dt;
	} else {
		// Assemble rhs only.
		timeDisc->assemble_rhs(*sp_rhs->get(), gridlevel);
	}

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
};

template<typename TDomain, typename TAlgebra>
braid_Int PBraidApp<TDomain,TAlgebra>::Residual(braid_Vector u, braid_Vector r, BraidStepStatus &pstatus)
{
	exit(20); // "Not functioning with residual" // todo
	return 0;
}

// todo replace?
template<typename TDomain, typename TAlgebra>
braid_Int PBraidApp<TDomain,TAlgebra>::SpatialNorm(braid_Vector u, braid_Real *norm_ptr)
{
	//a->setSpatialNorm(&ug::VecTwoNormSq<typename TAlgebra::vector_type>);
	//a->setSpatialNorm(&VecNorm2<typename TAlgebra::vector_type>);
	auto *uref = (SPGridFunction *) u->value;
	*norm_ptr = (*uref)->norm();
	return 0;
};

template<typename TDomain, typename TAlgebra>
braid_Int PBraidApp<TDomain,TAlgebra>::Access(braid_Vector u, BraidAccessStatus &astatus)
{
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
template<typename TDomain, typename TAlgebra>
braid_Int PBraidApp<TDomain,TAlgebra>::Coarsen(braid_Vector fu, braid_Vector *cu, BraidCoarsenRefStatus &status)
{
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
template<typename TDomain, typename TAlgebra>
braid_Int  PBraidApp<TDomain,TAlgebra>::Refine(braid_Vector cu, braid_Vector *fu, BraidCoarsenRefStatus &status)
{
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
#endif
