//
// Created by parnet on 12.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_PBRAIDAPP_H
#define UG_PLUGIN_XBRAIDFORUG4_PBRAIDAPP_H

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
class PBraidApp : public GFBraidApp<TDomain, TAlgebra> {

public:
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

    SPTimeStep timeDisc;
    SPSolver linSolver;
    SPTimeIntegrator m_spIntegratorC;
    SPTimeIntegrator m_spIntegratorF;
    SPOutput out;


    SPGridFunction ux; // for t > tstart
    SmartPtr<SpaceTimeCommunicator> comm;

    const char *filename{};

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



    // bool verbose = true; todo

public:

    PBraidApp(MPI_Comm mpi_temporal, double tstart, double tstop, int steps) :
            GFBraidApp<TDomain, TAlgebra>(mpi_temporal, tstart, tstop, steps) {}

    ~PBraidApp() override = default;


    void setStartVector(SPGridFunction p_u0) {
        this->u0 = p_u0;
    }

    void setRemainingVector(SPGridFunction p_ux) {
        this->ux = p_ux;
    }

    void setTimeDisc(SPTimeStep p_timeDisc) {
        this->timeDisc = p_timeDisc;
    }

    void setLinSolver(SPSolver p_linSolver) {
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
    virtual braid_Int Step(braid_Vector u, //
                           braid_Vector ustop, // estimated solution?
                           braid_Vector fstop,
                           BraidStepStatus &pstatus) {
        double tstart;
        double tstop;
        pstatus.GetTstartTstop(&tstart, &tstop);
        double current_dt = tstop - tstart;

        auto *sp_u_start = (SPGridFunction *) u->value;
        auto *sp_u_stop = (SPGridFunction *) u->value;
        auto *sp_u_stop_approx = (SPGridFunction *) ustop->value;

        auto *sp_rhs = new SPGridFunction(new TGridFunction(*this->u0));
        /*std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << pointerToString((void*)sp_rhs) << std::endl;
        std::cout << pointerToString((void*)sp_u_start) << std::endl;
        std::cout << pointerToString((void*)sp_u_stop) << std::endl;
        std::cout << pointerToString((void*)sp_u_stop_approx) << std::endl;*/



        //auto *rhs = (BraidVector *) malloc(sizeof(BraidVector));
        //auto *sp_rhs_val_ref = new SPGridFunction(new TGridFunction(*u0));
        //rhs->value = sp_rhs_val_ref;
        //fstop = rhs;

        const ug::GridLevel gridlevel = sp_u_stop->get()->grid_level();


        SPTimeSeries series = SPTimeSeries(new TTimeSeries());
        series->push(sp_u_start->get()->clone(), tstart);
        timeDisc->prepare_step(series, current_dt);

        if (current_dt != assembled_dt) { // todo is close?

            timeDisc->assemble_linear(*A, *sp_rhs->get(), gridlevel);
            linSolver->init(A, *sp_u_stop_approx->get());
            assembled_dt = current_dt;
        } else {
            timeDisc->assemble_rhs(*sp_rhs->get(), gridlevel);
        }

        /*if (fstop != nullptr) {
            auto *sp_rhs_stop = (SPGridFunction *) fstop->value;

            write(sp_u_stop_approx->get(),7055,0.01,"u_stop_approx");
            write(sp_u_start->get(),7055,0.01,"u_start");
            write(sp_rhs->get(),75055,0.01, "rhs");
            write(sp_rhs_stop->get(),75055,0.01,"rhs_stop");
            write(sp_u_stop->get(),75055,0.01,"u_stop");
            exit(80);
            VecAdd(1,*sp_u_stop_approx->get(),1,*sp_rhs_stop->get());
        }*/

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

    int step_ignoring_residual(braid_Vector u, //
                               braid_Vector ustop, // estimated solution?
                               braid_Vector fstop,
                               BraidStepStatus &pstatus) {
        double tstart;
        double tstop;
        pstatus.GetTstartTstop(&tstart, &tstop);
        double current_dt = tstop - tstart;

        auto *sp_u_start = (SPGridFunction *) u->value;
        auto *sp_u_stop = (SPGridFunction *) u->value;
        auto *sp_u_stop_approx = (SPGridFunction *) ustop->value;

        auto *sp_rhs = new SPGridFunction(new TGridFunction(*this->u0));
        /*std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << pointerToString((void*)sp_rhs) << std::endl;
        std::cout << pointerToString((void*)sp_u_start) << std::endl;
        std::cout << pointerToString((void*)sp_u_stop) << std::endl;
        std::cout << pointerToString((void*)sp_u_stop_approx) << std::endl;*/



        //auto *rhs = (BraidVector *) malloc(sizeof(BraidVector));
        //auto *sp_rhs_val_ref = new SPGridFunction(new TGridFunction(*u0));
        //rhs->value = sp_rhs_val_ref;
        //fstop = rhs;

        const ug::GridLevel gridlevel = sp_u_stop->get()->grid_level();


        SPTimeSeries series = SPTimeSeries(new TTimeSeries());
        series->push(sp_u_start->get()->clone(), tstart);
        timeDisc->prepare_step(series, current_dt);

        if (current_dt != assembled_dt) { // todo is close?

            timeDisc->assemble_linear(*A, *sp_rhs->get(), gridlevel);
            linSolver->init(A, *sp_u_stop_approx->get());
            assembled_dt = current_dt;
        } else {
            timeDisc->assemble_rhs(*sp_rhs->get(), gridlevel);
        }

        /*if (fstop != nullptr) {
            auto *sp_rhs_stop = (SPGridFunction *) fstop->value;

            write(sp_u_stop_approx->get(),7055,0.01,"u_stop_approx");
            write(sp_u_start->get(),7055,0.01,"u_start");
            write(sp_rhs->get(),75055,0.01, "rhs");
            write(sp_rhs_stop->get(),75055,0.01,"rhs_stop");
            write(sp_u_stop->get(),75055,0.01,"u_stop");
            exit(80);
            VecAdd(1,*sp_u_stop_approx->get(),1,*sp_rhs_stop->get());
        }*/

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
    }


    braid_Int Residual(braid_Vector u, braid_Vector r,
                       BraidStepStatus &pstatus) override {
        exit(20); // "Not functioning with residual" // todo
        return 0;

    };

    // todo replace?
    braid_Int SpatialNorm(braid_Vector u, braid_Real *norm_ptr) override {
        //a->setSpatialNorm(&ug::VecTwoNormSq<typename TAlgebra::vector_type>);
        //a->setSpatialNorm(&VecNorm2<typename TAlgebra::vector_type>);
        auto *uref = (SPGridFunction *) u->value;
        *norm_ptr = (*uref)->norm();
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

/*
      // delete function
    std::string pointerToString(void *c) {
        std::stringstream ss;
        ss << c;
        return ss.str();
    }
 */

// buffer pack

// <<<<<---------------------------------------------------------ug::Serialize(tmpbuffer,(vector_type&) *u_ref);

//std::cout << "original size\t" << s << std::endl;
//size_t rs;
//tmpbuffer->read((char *)&rs,sizeof(size_t));
//std::cout << "controll size\t" << rs << std::endl;
//char * buf = tmpbuffer->buffer();
//size_t * bufsz = (size_t *) buf;
//std::cout << "controll size*\t" << bufsz[0] << std::endl;

//std::cout << " === === === BGN SZV=== === ===" << std::endl;
//char * testbf = tmpbuffer->buffer();
//for(int i = 0; i < sizeof(size_t); i++){
//   std::cout << testbf[i] << std::endl;
//}
//std::cout << " === === === END SZV=== === ===" << std::endl;



// >>>>>--------------------------------------------------------ug::Serialize(tmpbuffer,(vector_type&) *u_ref);

/*
    inline void packWBuffer(void *buffer, TGridFunction *u_ref, int *bufferSize) {
        ug::BinaryBuffer handler = ug::BinaryBuffer();
        size_t szVector = u_ref->size();
        handler.write((char *) &szVector, sizeof(size_t));
        for (size_t i = 0; i < szVector; i++) {
            std::cout << i << "\t" << (*u_ref)[i] << std::endl << std::flush;
            handler.write((char *) &(*u_ref)[i], sizeof(double));
        }


        char *bufferOfBinaryBuffer = handler.buffer();//handler->detach((char *) buffer);
        memcpy(buffer, bufferOfBinaryBuffer, handler.capacity()); // copy is not neccessary
        size_t *chBuffer = (size_t *) buffer;
        *bufferSize = handler.capacity();

}

     inline void unpackWBuffer(void *buffer, TGridFunction *u_ref, int *bufferSize) {
        std::cout << "size(buf): " << *bufferSize << std::endl;
        ug::BinaryBuffer *handler = new ug::BinaryBuffer();
        handler->write((char *) buffer, *bufferSize);
        size_t szVector = 0;
        handler->read((char *) &szVector, sizeof(size_t));
        std::cout << "size(u): " << szVector << std::endl;
        for (size_t i = 0; i < szVector; i++) {
            double val;
            handler->read((char *) &val, sizeof(double));
            //ug::Deserialize(*tmpbuffer, val);
            (*u_ref)[i] = val;
            std::cout << i << "\t" << val << std::endl << std::flush;
        }

        free(handler);
}
 */

//int buffersize;
//BufSize(&buffersize, bstatus);
//std::cout << tmpbuffer->capacity() << std::endl;
//std::cout << buffersize << std::endl;
//size_t * sz_buffer = (size_t *) buffer;
//std::cout << "RBuffer Size*\t" << sz_buffer[0] << std::endl;
//std::cout << " === === === BGN TST=== === ===" << std::endl;
//char * scndTBuffer =(char*) buffer;
//for(int i = 0; i < sizeof(size_t); i++){
//  std::cout << scndTBuffer[i] << std::endl;
//}
//std::cout << " === === === END TST=== === ===" << std::endl;
//exit(207);
//free(buffer);
//exit(20);
/*
write(u_ref,cntS,0,"send"); cntS++;
auto *vec_ref = (vector_type *) u_ref;
memcpy(buffer, &(*vec_ref)[0], vec_ref->size());
*/

//int * intbuffer = (int*)buffer;
//intbuffer[0] = this->comm->getGlobalRank();
//intbuffer[1] = this->comm->getSpatialRank();
//intbuffer[2] = this->comm->getTemporalRank();

/*this->comm->print("Start Buffer pack");
auto *u_refa = (TGridFunction *) u->value;
this->comm->print(std::to_string(u_refa->dim));
this->comm->print(std::to_string(u_refa->size()));
this->comm->print(std::to_string(u_refa->num_indices()));

if(this->comm->getSpatialSize() > 1) {
    if (this->comm->getSpatialRank() == 0) {
        this->comm->print(std::to_string(u_refa->num_indices(1)));
    } else {
        this->comm->print(std::to_string(u_refa->num_indices(0)));
    }
}
const vector_type * vct =  ug::getVector(u_refa);
this->comm->print(std::to_string(vct->size()));

exit(-1);
this->comm->print("Start Buffer pack");
auto *u_ref = (TGridFunction *) u->value;
this->comm->print("BufPack   1");
std::cout << (*u_ref).norm() << std::endl;
this->comm->print("BufPack   2");
auto v = static_cast <vector_type &>(*u_ref);
this->comm->print("BufPack   3");
auto *dbuffer = (vector_type *) buffer;
this->comm->print("BufPack   4");
new((void *) &dbuffer[0]) vector_type(v);
this->comm->print("BufPack   5");
bstatus.SetSize(sizeof(v));
auto *dbuffer = (TGridFunction *) buffer;
auto *u_ref = (TGridFunction *) u->value;
new((void *) &dbuffer[0]) TGridFunction(*u_ref);
bstatus.SetSize(sizeof(TGridFunction));
this->comm->print("End Buffer Pack");*/
/*if(comm->getGlobalRank() == 0){
    std::cout << "Enter Process" << std::endl;
    comm->processPrint();
    std::cout << "Leaving Process" << std::endl;
}*/

// UNPACK


//DetachableBuffer *tmpbuffer = new DetachableBuffer(buffersize);

//buffer = tmpbuffer->attachAndHold((char *) buffer,buffersize);

//<<<<<<<<<<<-----------------------------------------------------------------ug::Deserialize(tmpbuffer, u_ref);


//ug::Deserialize(*tmpbuffer, s);
//u_ref->operator[](7) = s;

//(*u_ref).resize(s); // should have the needed size


//>>>>>>>>>>>-----------------------------------------------------------------ug::Deserialize(tmpbuffer, u_ref);
//exit(127);


//exit(-99);

//write(u_ref,cntR,0,"recv"); cntR++;
/*
auto *vec_ref = (vector_type *) u_ref;
memcpy(&(*vec_ref)[0], buffer, vec_ref->size());
u->value = u_ref;
*u_ptr = u;
write(u_ref,cntR,0,"recv"); cntR++;
*/
/*std::cout <<  global_rank << "\tStart Buffer unpack" << std::endl;
this->comm->print("Start Buffer Unpack");
BraidVector *u;
this->comm->print("BufUnpack   1");
u = (BraidVector *) malloc(sizeof(BraidVector));
this->comm->print("BufUnpack   2");
auto *dbuffer = (vector_type *) buffer;
this->comm->print("BufUnpack   3");
auto *u_ref = new TGridFunction(*u0);
this->comm->print("BufUnpack   4");
auto uu_ref = new TGridFunction(*u_ref);
this->comm->print("BufUnpack   4.1");

this->comm->print("BufUnpack   4.1");

u_ref->assign(*dbuffer);
this->comm->print("BufUnpack   5");

u->value = u_ref;
this->comm->print("BufUnpack   6");
std::cout << (*u_ref).norm() << std::endl;
this->comm->print("BufUnpack   7");
*u_ptr = u;

this->comm->print("BufUnpack   8");
auto *dbuffer = (TGridFunction *) buffer;
BraidVector *u;
u = (BraidVector *) malloc(sizeof(BraidVector));
auto *u_ref = new TGridFunction(*dbuffer);
u->value = u_ref;
*u_ptr = u;
//std::cout <<  global_rank << "\tEnd Buffer unpack";
this->comm->print("End Buffer unpack");*/



#endif //UG_PLUGIN_XBRAIDFORUG4_PBRAIDAPP_H
