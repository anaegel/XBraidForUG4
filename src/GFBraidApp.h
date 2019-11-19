//
// Created by parnet on 24.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_GFBRAIDAPP_H
#define UG_PLUGIN_XBRAIDFORUG4_GFBRAIDAPP_H



// Std lib.
#include <iomanip>
#include <sstream>


// This lib.
#include "trace_tools_config.h"
#include "BraidVectorStruct.h"
#include "SpaceTimeCommunicator.h"
#include "Scriptor.h"
#include "Talasma.h"
#include "Telvanni.h"
// #include "MemoryObserver.h"

#if TRACE_GRIDFUNCTION == 1
    #define MATLAB(u,i,t) this->matlab->write(u,i,t)
#else
    #define MATLAB(u, i, t)
#endif


template<typename TDomain, typename TAlgebra>
class GFBraidApp : public BraidApp {
public:
    // -----------------------------------------------------------------------------------------------------------------
    // type definition to shorten identifier
    // -----------------------------------------------------------------------------------------------------------------
	typedef GFBraidApp<TDomain, TAlgebra> this_type;

	typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef typename TAlgebra::vector_type::value_type vector_value_type;
    typedef typename ug::XBraidForUG4::SpaceTimeCommunicator TSpaceTimeCommunicator;

    typedef SmartPtr<TGridFunction> SPGridFunction;
    typedef SmartPtr<ug::UserData<double, TGridFunction::dim> > SPData;
    typedef SmartPtr<ug::IDomainDiscretization<TAlgebra> > SPDomainDisc;
    typedef SmartPtr<Scriptor<TDomain, TAlgebra> > SPScriptor;

    typedef typename ug::XBraidForUG4::BraidVectorFunctors<TGridFunction> TBraidVectorFunctors;

    SmartPtr <TSpaceTimeCommunicator> m_comm;
protected:
    // -----------------------------------------------------------------------------------------------------------------
    // members for vector creation / initialization
    // -----------------------------------------------------------------------------------------------------------------

    const char *name;
    SPGridFunction m_u0; // for t = tstart
    SPData m_data;
    SPDomainDisc m_domainDisc; // for adjust gridfunction in generator
    const char *m_cmp;  // function component for interpolation

    bool m_timing = true;
    bool m_verbose = true;
    bool m_writeparam = true;

    SPScriptor m_out;


    std::ofstream debugwriter;
    TraceTools::Talasma timer;

#if TRACE_TIMINGS == 1
    TraceTools::Redoran redoran;
#endif

#if TRACE_GRIDFUNCTION == 1
    SmartPtr<MATLABScriptor<TDomain, TAlgebra>> matlab;
#endif
    int m_levels = 15;

    /**
     * Note that this default constructor does not create a consistent object. The parameter t_comm (of type MPI_Comm)
     * for the temporal communication has to be set.
     */
public:
    GFBraidApp() : BraidApp(MPI_COMM_NULL, 0.0, 10.0, 10) {}

    GFBraidApp(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
    : BraidApp(mpi_temporal, tstart, tstop, steps)
    {}

    ~GFBraidApp() {} ;


    // -----------------------------------------------------------------------------------------------------------------
    // Values for time setting
    // -----------------------------------------------------------------------------------------------------------------
    void setTimeValues(double startTime, double endTime, int n) {
        this->tstart = startTime;
        this->tstop = endTime;
        this->ntime = n;
    }

    void setStartTime(double startTime) {
        this->tstart = startTime;
    }

    void setEndTime(double endTime) {
        this->tstop = endTime;
    }

    void setNumberOfTimesteps(int n) {
        this->ntime = n;
    }

    // -----------------------------------------------------------------------------------------------------------------
    // Set Grid functions
    // -----------------------------------------------------------------------------------------------------------------
    void setStartVector(SPGridFunction p_u0) {
        this->m_u0 = p_u0;
    }


    void setGeneratorComponent(const char *cmp) {
        this->m_cmp = cmp;
    }


    void setVectorGenerator(SPData p_data) {
        this->m_data = p_data;
    }

    void setDomainDisc(SPDomainDisc p_domainDisc) {
        m_domainDisc = p_domainDisc;
    }

#ifdef UG_FOR_LUA

    void setVectorGenerator(const char *fctName) {
        setVectorGenerator(ug::LuaUserDataFactory<double, TDomain::dim>::create(fctName));
    }

    void setVectorGenerator(ug::LuaFunctionHandle fct) {
        setVectorGenerator(make_sp(new ug::LuaUserData<double, TDomain::dim>(fct)));
    }

#endif


    void setVerbose(bool p_verbose) {
        this->m_verbose = p_verbose;
    }

    void setMaxLevels(size_t levelcount) {
#if TRACE_TIMINGS == 1
        std::cerr << "Warning: Tracing times is currently limited to 15 level" << std::endl;
#endif
        this->m_levels = levelcount;
    }

    virtual void init() = 0;

    void release() {
#if TRACE_TIMINGS
    {
    	using namespace TraceTools;
        for (Observer i = T_INIT; i != T_RECV; i++) {
            Telvanni &tel = redoran.get(i);
            this->debugwriter << std::setw(20) << ObserverNames[i]
                              << std::setw(5) << " ;"
                              << std::setw(12) << tel.getTime() << ";"
                              << std::setw(12) << tel.getUsage() << ";"
                              << std::setw(12) << (tel.getTime() / tel.getUsage()) << ";"
                              << std::endl;
        }

        for (int i = 0; i != cLevelObserver; i++) {
            for (int l = 0; l < this->m_levels; l++) {
                Telvanni &tel = redoran.get(static_cast<LevelObserver >(i), l);
                this->debugwriter << std::setw(20) << LevelObserverNames[i]
                                  << std::setw(5) << l << ";"
                                  << std::setw(12) << tel.getTime() << ";"
                                  << std::setw(12) << tel.getUsage() << ";"
                                  << std::setw(12) << (tel.getTime() / tel.getUsage()) << ";"
                                  << std::endl;
            }
        }
    }
#endif
    }



    // -----------------------------------------------------------------------------------------------------------------
    // Access
    // -----------------------------------------------------------------------------------------------------------------

    void setScriptor(SPScriptor p_out) {
        this->m_out = p_out;
    }



    // -----------------------------------------------------------------------------------------------------------------
    // Braid Methods
    // -----------------------------------------------------------------------------------------------------------------

    //! Allocate a new vector in @a *u_ptr and initialize it with an initial guess appropriate for time @a t. */
    braid_Int Init(braid_Real t, braid_Vector *u_ptr) override {

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
    braid_Int Clone(braid_Vector u, braid_Vector *v_ptr) override {
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
    braid_Int Free(braid_Vector u) override {
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
    braid_Int Sum(double alpha, braid_Vector x, double beta, braid_Vector y) override {
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
       // SPGridFunction *xref = (SPGridFunction *) x->value;
       // SPGridFunction *yref = (SPGridFunction *) y->value;
       // VecAdd(beta, *(yref->get()), alpha, *(xref->get()));

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
    braid_Int Access(braid_Vector u, BraidAccessStatus &astatus) override {
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

protected:
    inline void pack(void *buffer, TGridFunction *u_ref, int *bufferSize) {

        char *chBuffer = (char *) buffer;
        *bufferSize = 0;
        size_t szVector = u_ref->size();

        memcpy(buffer, &szVector, sizeof(size_t));
        *bufferSize += sizeof(size_t);
        for (size_t i = 0; i < szVector; i++) {
            memcpy(chBuffer + *bufferSize, &(*u_ref)[i], sizeof(vector_value_type));
            *bufferSize += sizeof(vector_value_type);
        }

        int temprank = this->m_comm->getTemporalRank();

        memcpy(chBuffer + *bufferSize, &temprank, sizeof(int));
        *bufferSize += sizeof(int);

    }

    inline void unpack(void *buffer, TGridFunction *u_ref, int *bufferSize) {

        char *chBuffer = (char *) buffer;
        size_t szVector = 0;
        memcpy(&szVector, chBuffer, sizeof(size_t));
        int pos = sizeof(size_t);
        for (size_t i = 0; i < szVector; i++) {
            double val = 0;
            memcpy(&val, chBuffer + pos, sizeof(vector_value_type));
            pos += sizeof(vector_value_type);
            (*u_ref)[i] = val;
        }

        int temprank;
        memcpy(&temprank, chBuffer + pos, sizeof(int));
        pos += sizeof(int);

#if TRACE_INDEX == 1
        int index;
        memcpy(&index, chBuffer + pos, sizeof(int));
        pos += sizeof(int);
        if (this->m_verbose) {
            debugwriter << "rec( v_" << temprank << "_" << index << ")" << std::endl;
        }
#endif

    }

public:
    //! Define buffer size.
    braid_Int BufSize(braid_Int *size_ptr, BraidBufferStatus &bstatus) override {
        *size_ptr = (sizeof(vector_value_type) * (*this->m_u0).size()
                     + sizeof(size_t))
                    + 2 * sizeof(int);
        return 0;
    };

    //! Packing buffer.
    braid_Int BufPack(braid_Vector u, void *buffer,
                      BraidBufferStatus &bstatus) override {
    	// Profiling BEGIN
    	StartRedoran(Observer::T_SEND);

#if TRACE_INDEX == 1
        if (this->m_verbose) {
            debugwriter << "send(u_" << u->index << ")" << std::endl << std::flush;
        }
#endif

        int bufferSize;
        auto *u_ref = (SPGridFunction *) u->value;
        pack(buffer, u_ref->get(), &bufferSize);

#if TRACE_INDEX == 1
        char *chBuffer = (char *) buffer;
        memcpy(chBuffer + bufferSize, &u->index, sizeof(int));
        bufferSize += sizeof(int);
#endif
        bstatus.SetSize(bufferSize);

        // Profiling END
        StopRedoran(Observer::T_SEND);
#if TRACE_RECVTIME == 1
        double diff, total;
        this->timer.now(total, diff);
        this->debugwriter << std::setw(10) << "@time:"
                          << std::setw(12) << total << " ; "
                          << std::setw(12) << diff << " Vector Send" << std::endl;
#endif
        return 0;
    };

    //! Unpacking buffer.
    braid_Int BufUnpack(void *buffer, braid_Vector *u_ptr, BraidBufferStatus &bstatus) override {
#if TRACE_RECVTIME == 1
        double diff, total;
        this->timer.now(total, diff);
        this->debugwriter << std::setw(10) << "@time:"
                          << std::setw(12) << total << " ; "
                          << std::setw(12) << diff << " Vector Received" << std::endl;
#endif
        // Profiling BEGIN.
        StartRedoran(Observer::T_RECV);
#if TRACE_INDEX == 1
        if (this->m_verbose) {
            this->debugwriter << "u_" << indexpool << " = ";
        }
#endif

        int bufferSize;
        BufSize(&bufferSize, bstatus);
        auto *u = (BraidVector *) malloc(sizeof(BraidVector));
        auto *sp_u = new SPGridFunction(new TGridFunction(*this->m_u0)); // todo
        unpack(buffer, sp_u->get(), &bufferSize);
        u->value = sp_u;

#if TRACE_INDEX == 1
        u->index = indexpool;
        indexpool++;
        MATLAB(sp_u->get(),u->index,-1.0);
#endif
        *u_ptr = u;

        // Profiling END.
        StopRedoran(Observer::T_RECV);
        return 0;
    };


};

#endif //UG_PLUGIN_XBRAIDFORUG4_GFBRAIDAPP_H
