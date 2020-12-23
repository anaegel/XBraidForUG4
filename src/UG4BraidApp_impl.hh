#ifndef UG_PLUGIN_XBRAIDFORUG4_UG4BRAIDAPP_IMPL_HH
#define UG_PLUGIN_XBRAIDFORUG4_UG4BRAIDAPP_IMPL_HH


template<typename TDomain, typename TAlgebra>
braid_Int UG4BraidApp<TDomain,TAlgebra>::Init(braid_Real t, braid_Vector *u_ptr) {

#if TRACE_INDEX == 1
     if (this->m_verbose) {
         this->debugwriter << "u_" << indexpool << " = init(" << t << ")" << std::endl;
     }
#endif
     StartRedoran(Observer::T_INIT);
     auto *u = (BraidVector *) malloc(sizeof(BraidVector));
     SPGridFunction *vec = new SPGridFunction();
     if (t == this->tstart) {
         *vec = this->m_u0->clone();
     } else {
         *vec = this->m_u0->clone_without_values();
         Interpolate(this->m_data, *vec, this->m_cmp, NULL, t);
         m_domainDisc->adjust_solution(*vec->get(), t);
     }

     u->value = vec;


#if TRACE_INDEX == 1
     u->index = indexpool;
     indexpool++;
     MATLAB(vec->get(), u->index, t);
#endif
     *u_ptr = u;

     StopRedoran(Observer::T_INIT);
     return 0;
 };



//! Allocate a new vector in @a *v_ptr, which is a deep copy of @a u.
template<typename TDomain, typename TAlgebra>
braid_Int UG4BraidApp<TDomain,TAlgebra>::Clone(braid_Vector u, braid_Vector *v_ptr)  {
	// Profiling BEGIN.
    StartRedoran(Observer::T_CLONE);
#if TRACE_INDEX == 1
    if (this->m_verbose) {
        this->debugwriter << "u_" << indexpool << " = clone(u_" << u->index << ")" << std::endl;
    }
#endif

    // Create (two-stage).
    SPGridFunction *uref = (SPGridFunction *) u->value;  // Ptr to (existing) SmartPtr.
    SPGridFunction *vref = new SPGridFunction();         // STEP A: Create (invalid) SmartPtr object.
    (*vref) = uref->get()->clone();       				 // Create a new GridFunction object; assign its SmartPtr.

    BraidVector *v = (BraidVector *) malloc(sizeof(BraidVector)); // STEP B: Create a ptr to a (new) BraidVector.
    v->value = vref; // Store ptr

#if TRACE_INDEX == 1
    v->index = indexpool;
    indexpool++;
    MATLAB(vref->get(), v->index, -1.0);
#endif

    // Assign return value.
    *v_ptr = v;

    // Profiling END.
    StopRedoran(Observer::T_CLONE);
    return 0;
};

//! De-allocate the vector @a u.
template<typename TDomain, typename TAlgebra>
braid_Int UG4BraidApp<TDomain,TAlgebra>::Free(braid_Vector u) {
	// Profiling BEGIN.
    StartRedoran(Observer::T_FREE);

#if TRACE_INDEX == 1
    if (this->m_verbose) {
        this->debugwriter << "u_" << u->index << " = null" << std::endl;
    }
#if TRACE_CONST == 1
        if (u->m_const) {
            this->debugwriter << "u_" << u->index << " was const" << std::endl;
            const_free++;
        }
#endif
#endif

    // Delete (two-stage).
    SPGridFunction *u_value = (SPGridFunction *) u->value;
    delete u_value;  // STEP A: Delete SmartPtr object (which deletes real object, if applicable).
    free(u);         // STEP B: Delete BraidVector

    // Profiling END.
    StopRedoran(Observer::T_FREE);
    return 0;
};

//! Perform the operation: @a y = @a alpha * @a x + @a beta * @a y.
template<typename TDomain, typename TAlgebra>
braid_Int UG4BraidApp<TDomain,TAlgebra>::Sum(double alpha, braid_Vector x, double beta, braid_Vector y)  {
	// Profiling BEGIN.
	StartRedoran(Observer::T_SUM);
#if TRACE_INDEX == 1
    if (this->m_verbose) {
        if (alpha == 0) {
            this->debugwriter << "u_" << y->index << " = " << beta << "* u_" << y->index << " % Scale "
                              << std::endl;
        } else if (beta == 0) {
            this->debugwriter << "u_" << y->index << " = " << alpha << "*u_" << x->index << "  % Replace "
                              << std::endl;
        } else {
            this->debugwriter << "u_" << y->index << " = " << alpha << "* u_" << x->index << "  + " << beta
                              << "* u_"
                              << y->index << " % Sum " << std::endl;
        }
    }
#endif
#if TRACE_CONST == 1
    y->m_const = false;
#endif

    const TGridFunction& xvec = TBraidVectorFunctors().as_grid_function(*x);
    TGridFunction& yvec = TBraidVectorFunctors().as_grid_function(*y);
    VecAdd(beta, yvec, alpha, xvec);

    // Profiling END.
    StopRedoran(Observer::T_SUM);
#if TRACE_INDEX ==1
    MATLAB(yref->get(), y->index, -1.0);
#endif
    return 0;
};


//! Some output.
template<typename TDomain, typename TAlgebra>
  braid_Int UG4BraidApp<TDomain,TAlgebra>::Access(braid_Vector u, BraidAccessStatus &astatus)  {
      StartRedoran(Observer::T_ACCESS);
#if TRACE_INDEX == 1
      if (this->m_verbose) {
          this->debugwriter << "% \t Access \t" << u->index << std::endl;
      }
#endif

      int v = 0;
      int index;
      astatus.GetTIndex(&index);
      double timestamp;
      astatus.GetT(&timestamp);

      auto *ref = (SPGridFunction *) u->value;


      int iter;
      int lvl;
      int done;

      astatus.GetIter(&iter);
      astatus.GetLevel(&lvl);
      astatus.GetDone(&done);
#if TRACE_ACCESS == 1

      if(done == 1){
          v = this->m_out->write(ref->get(), index, timestamp);
      } else {
          v = this->m_out->write(ref->get(), index, timestamp, iter, lvl);
      }
#else
      v = this->m_out->write(ref->get(), index, timestamp);
#endif
      StopRedoran(Observer::T_ACCESS);
      return v;

  };

#endif
