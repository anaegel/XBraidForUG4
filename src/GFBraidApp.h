//
// Created by parnet on 24.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_GFBRAIDAPP_H
#define UG_PLUGIN_XBRAIDFORUG4_GFBRAIDAPP_H


#include "BraidVectorStruct.h"
#include "XCommunicator.h"

template<typename TDomain, typename TAlgebra>
class GFBraidApp : public BraidApp {
public:
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr<TGridFunction> SPGridFunction;

    SPGridFunction u0; // for t = tstart

    bool verbose = true;

    void setVerbose(bool p_verbose){
        this->verbose = p_verbose;
    }

    typedef double vector_value_type;

    GFBraidApp(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) : BraidApp(mpi_temporal, tstart, tstop,
                                                                                         steps) {}


    braid_Int Clone(braid_Vector u, braid_Vector *v_ptr) override {

        auto *v = (BraidVector *) malloc(sizeof(BraidVector));

        std::cout << "Cloning Vector @t="<<u->time<<"\t" << u->index << "\t to \t" << indexpool <<  std::endl;

        auto *uref = (SPGridFunction *) u->value;
        v->value = new SPGridFunction(new TGridFunction(*(*uref)));
        v->index = indexpool;
        v->time = u->time;

        *v_ptr = v;
        indexpool++;
        return 0;
    };

    braid_Int Free(braid_Vector u) override {
        std::cout << "Free Vector @t=" << u->time << "\t" << u->index <<std::endl;
        free(u->value);
        free(u);
        return 0;
    };

    braid_Int Sum(double alpha, braid_Vector x, double beta, braid_Vector y) override {
        auto *xref = (SPGridFunction *) x->value;
        auto *yref = (SPGridFunction *) y->value;
        VecAdd(beta, *yref->get(), alpha, *xref->get());
        return 0;
    };

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
    }

    braid_Int BufSize(braid_Int *size_ptr, BraidBufferStatus &bstatus) override {
        *size_ptr = (sizeof(vector_value_type) * (*u0).size() + sizeof(size_t));
        return 0;
    };


    braid_Int BufPack(braid_Vector u, void *buffer,
                      BraidBufferStatus &bstatus) override {
        auto *u_ref = (SPGridFunction *) u->value;
        int bufferSize;
        pack(buffer, u_ref->get(), &bufferSize);
        bstatus.SetSize(bufferSize);
        return 0;
    };


    braid_Int BufUnpack(void *buffer, braid_Vector *u_ptr, BraidBufferStatus &bstatus) override {
        int bufferSize;
        BufSize(&bufferSize, bstatus);
        auto *u = (BraidVector *) malloc(sizeof(BraidVector));
        auto *sp_u = new SPGridFunction(new TGridFunction(*this->u0));
        unpack(buffer, sp_u->get(), &bufferSize);
        u->value = sp_u;
        u->index = indexpool; indexpool++;
        *u_ptr = u;
        return 0;
    };
};

#endif //UG_PLUGIN_XBRAIDFORUG4_GFBRAIDAPP_H