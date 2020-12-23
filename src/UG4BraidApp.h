//
// Created by parnet on 24.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_GFBRAIDAPP_H
#define UG_PLUGIN_XBRAIDFORUG4_GFBRAIDAPP_H



// Std lib.
#include <iomanip>
#include <string>
#include <sstream>


// This lib.
#include "BraidVectorStruct.h"
#include "SpaceTimeCommunicator.h"
#include "Scriptor.h"

// TraceTools
#include "tools/trace_tools_config.h"
#include "tools/Talasma.h"
#include "tools/Telvanni.h"

#if TRACE_GRIDFUNCTION == 1
    #define MATLAB(u,i,t) this->matlab->write(u,i,t)
#else
    #define MATLAB(u, i, t)
#endif

/**This is a 'BraidApp' that deals with UG4-type grid functions. */
template<typename TDomain, typename TAlgebra>
class UG4BraidApp : public BraidApp {
public:
    // -----------------------------------------------------------------------------------------------------------------
    // type definition to shorten identifier
    // -----------------------------------------------------------------------------------------------------------------
	typedef UG4BraidApp<TDomain, TAlgebra> this_type;

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

    std::string name;
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


public:
    /** Note that this default constructor does not create a consistent object!
     *  The parameter t_comm (of type MPI_Comm) for the temporal communication must be set.
     */
    UG4BraidApp()
	: BraidApp(MPI_COMM_NULL, 0.0, 10.0, 10)
	{}

    UG4BraidApp(MPI_Comm mpi_temporal, double tstart, double tstop, int steps)
    : BraidApp(mpi_temporal, tstart, tstop, steps)
    {}

    ~UG4BraidApp() {} ;


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
    braid_Int Init(braid_Real t, braid_Vector *u_ptr) override;

    //! Allocate a new vector in @a *v_ptr, which is a deep copy of @a u.
    braid_Int Clone(braid_Vector u, braid_Vector *v_ptr) override;

    //! De-allocate the vector @a u.
    braid_Int Free(braid_Vector u) override;

    //! Perform the operation: @a y = @a alpha * @a x + @a beta * @a y.
    braid_Int Sum(double alpha, braid_Vector x, double beta, braid_Vector y) override;

    //! Some output.
    braid_Int Access(braid_Vector u, BraidAccessStatus &astatus) override;

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

#include "UG4BraidApp_impl.hh"

#endif //UG_PLUGIN_XBRAIDFORUG4_GFBRAIDAPP_H


