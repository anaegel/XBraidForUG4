//
// Created by parnet on 12.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_PBRAID_H
#define UG_PLUGIN_XBRAIDFORUG4_PBRAID_H

#include "common/assert.h"
#include "pcl/pcl_comm_world.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/time_disc/solution_time_series.h"

#include "PBraidApp.h"
#include "RGFBraidApp.h"
#include "XCommunicator.h"

#include "../libs/xbraid/braid/_braid_status.h"

template<typename TDomain, typename TAlgebra>
class PBraid {

    typedef ug::ThetaTimeStep<TAlgebra> TTimeStep;
    typedef ug::ILinearOperatorInverse<typename TAlgebra::vector_type> TSolver;
    typedef typename TAlgebra::vector_type vector_type;
private:
    RGFBraidApp<TDomain, TAlgebra> *a;

    SmartPtr<XCommunicator> comm;

    SmartPtr<ug::VTKOutput<TDomain::dim>> out;
    const char *filename;
public:
    typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;


    PBraid() {
        a = new RGFBraidApp<TDomain, TAlgebra>(PCL_COMM_WORLD, 0, 1, 10);
    }


    void setStartTime(double startTime) {
        a->tstart = startTime;
    }

    void setEndTime(double endTime) {
        a->tstop = endTime;
    }

    void setNumberOfTimesteps(int n) {
        a->ntime = n;
    }

    void setStartVector(SmartPtr<TGridFunction> u0) {
        a->setStartVector(u0);
    }

    void setRemainingVector(SmartPtr<TGridFunction> ux) {
        a->setRemainingVector(ux);
    }

    void setTimeDisc(SmartPtr<TTimeStep> timeDiscretization) {
        a->setTimeDisc(timeDiscretization);
    }

    void setLinearSolver(SmartPtr<TSolver> linSolver) {
        a->setLinSolver(linSolver);
    }

    void setCommunicator(SmartPtr<XCommunicator> comm) {
        this->comm = comm;
        a->comm = this->comm;
        a->comm_t = this->comm->TEMPORAL;
    }

    bool testResidual(RGFBraidApp<TDomain, TAlgebra> *a, double t, double dt) {
        braid_Vector u, unext, ustop, fstop;
        braid_Real result1;
        braid_Int myid_x, result_int;
        braid_Status status = _braid_CTAlloc(_braid_Status, 1);
        braid_AccessStatus astatus = (braid_AccessStatus) status;
        braid_StepStatus sstatus = (braid_StepStatus) status;

        BraidAccessStatus bas = BraidAccessStatus(astatus);
        BraidStepStatus bss = BraidStepStatus(sstatus);

        _braid_StepStatusInit(t, t + dt, 0, 1e-16, 0, 0, 0, 2, sstatus);
        _braid_AccessStatusInit(t, 0, 0.0, 0, 0, 0, 2, 0, 1, -1, astatus);

        /* Print intro */
        std::cout << "\nStarting braid_TestResidual\n" << std::endl;

        /* Test 1 */
        std::cout << "   braid_TestResidual:   Starting Test 1" << std::endl;
        std::cout << "   braid_TestResidual:   u = init(" << t << ")" << std::endl;
        a->Init(t, &u);
        std::cout << "   braid_TestResidual:   spatialnorm(u)" << std::endl;

        a->SpatialNorm(u, &result1);
        std::cout << "   braid_TestResidual:   " << result1<< std::endl;
        if (fabs(result1) == 0.0) {
            std::cout << "   braid_TestResidual:   Warning:  spatialnorm(u) = 0.0" << std::endl;
        } else if (braid_isnan(result1)) {
            std::cout << "   braid_TestResidual:   Warning:  spatialnorm(u) = nan" << std::endl;
        }

        std::cout << "   braid_TestResidual:   unext = clone(u)" << std::endl;
        a->Clone(u, &unext);

        std::cout << "   braid_TestResidual:   ustop = clone(u)"<<std::endl;
        a->Clone(u, &ustop);

        std::cout << "   braid_TestResidual:   fstop = clone(u)"<<std::endl;
        a->Clone(u, &fstop);

        std::cout << "   braid_TestResidual:   fstop = ustop + fstop"<<std::endl;
        a->Sum(1.0, ustop, 1.0, fstop);

        std::cout << "   braid_TestResidual:   unext = step(ustop, fstop, unext)"<<std::endl;
        //a->Step(unext, ustop, nullptr, bss);

        std::cout << "   braid_TestResidual:   ustop = clone(u)"<<std::endl;
        a->Clone(u, &ustop);

        std::cout << "   braid_TestResidual:   fstop = clone(u)"<<std::endl;
        a->Clone(u, &fstop);

        std::cout << "   braid_TestResidual:   fstop = ustop + fstop"<<std::endl;
        a->Sum(1.0, ustop, 1.0, fstop);

        std::cout << "   braid_TestResidual:   r = residual(unext, u)"<<std::endl;
        a->Residual(unext, u, bss);

        std::cout << "   braid_TestResidual:   r = fstop - r "<<std::endl;
        a->Sum(1.0, fstop, -1.0, u);

        std::cout << "   braid_TestResidual:   spatialnorm(r)"<<std::endl;
        a->SpatialNorm(u, &result1);


        std::cout << "   braid_TestResidual:   access(r) "<<std::endl;
        a->Access(u, bas);


        /* We expect the result to be close to zero */
        result_int = (braid_Int) round(-log10(result1));
        std::cout << "   braid_TestResidual:   actual output:    spatialnorm(r) approx." << result_int << std::endl;
        std::cout << "   braid_TestResidual:   actual output:    spatialnorm(r) approx." << result1 << std::endl;
        std::cout << "   braid_TestResidual:   expected output:  spatialnorm(r) approx. 0.0 " << std::endl << std::endl;

        /* Free variables */
        std::cout << "   braid_TestResidual:   free(u)" << std::endl;
        a->Free(u);

        std::cout << "   braid_TestResidual:   free(unext)" << std::endl;
        a->Free(unext);

        std::cout << "   braid_TestResidual:   free(ustop)" << std::endl;
        a->Free(ustop);

        std::cout << "   braid_TestResidual:   free(fstop)" << std::endl;
        a->Free(fstop);

        _braid_StatusDestroy(status);

        std::cout << "Finished braid_TestResidual" << std::endl;

    }

    int test() {
        BraidUtil bu = BraidUtil();
        FILE *file;
        file = fopen("myfile", "w");

        MPI_Comm comm_x = this->comm->SPATIAL;
        a->init();

        /*std::cout << "Start Test Init, Access" << std::endl;
        bu.TestInitAccess(a, comm_x, file, 1.0);
        std::cout << "Start Test Clone" << std::endl;
        bu.TestClone(a, comm_x,file,1.0);
        std::cout << "Start Test Sum" << std::endl;
        bu.TestSum(a, comm_x,file,1.0);
        std::cout << "Start Test Spatial Norm" << std::endl;
        bu.TestSpatialNorm(a, comm_x,file,1.0);
        std::cout << "Start Test Buffer" << std::endl;
        bu.TestBuf(a, comm_x, file, 0.0);*/
        std::cout << "Start Residual" << std::endl;
        //bu.TestResidual(a, comm_x, file, 1.0, 0.1);
        this->testResidual(a, 0.0, 0.1);
        std::cout << "[[[[Start Test Coarsen and Refine]]]]" << std::endl;
        // todo bu.TestCoarsenRefine()
        fclose(file);
        std::cout << "End Test" << std::endl;

        return 1;
    }

    void createAccess(SmartPtr<ug::VTKOutput<TDomain::dim>> out) {
        this->out = out;
        a->out = this->out;
    }

    void setFilename(const char *filename) {
        a->filename = filename;
        this->filename = filename;
    }

    void setVerbose(bool flag) {
        //this->a->verbose = flag;
    }

    int run() {
        a->init();

        //int rank;
        //MPI_Comm_rank(comm->GLOBAL, &rank);
        BraidCore *bc = new BraidCore(comm->GLOBAL, a);

        /*/////////////////////////////////////////////////////
        bc->Drive(); // run simulation
        bc->SetAbsTol(); // absolute stopping tolerance
        bc->SetRelTol(); // relative stopping tolerance

        bc->SetAggCFactor()
        bc->SetCFactor(); // coarsening factor (*2) with -1 address every level

        bc->SetMinCoarse()
        bc->SetMaxIter() // maximum number of iterations
        bc->SetMaxLevels() // max num of levels
        bc->SetMaxRefinements() // number of refinements


        bc->GetNumIter()


        bc->SetTemporalNorm(); // Norm 1-norm 2-norm, oo-norm
        bc->SetResidual(); // function pointer residual routine
        bc->SetRefine() // time refinement on = 1, off = 0


        bc->SetStorage(); // store points on level (0 for every level, -1 only C points)
        bc->SetSeqSoln()
        bc->SetSkip() // skip all work on the first down cycle (default skip = 1)
        bc->SetSpatialCoarsenAndRefine(); // coarse, refine

        bc->SetNRelax()

        bc->SetFMG(); // FMG - Cycle
        bc->SetNFMG(); // FMG - Cycle n times
        bc->SetNFMGVcyc(); // FMGV Cycle

        bc->SetAccessLevel()
        bc->SetPrintFile(); // 0 no output
        bc->SetPrintLevel(); // 0 no output

        bc->GetRNorms();
        bc->GetNLevels();
        */////////////////////////////////////////////////////
        bc->SetPrintLevel(2);
        bc->SetMaxLevels(7);
        bc->SetMaxIter(8);
        //bc->SetSeqSoln(1);

        //bc->SetRelTol();

        //bc->SetAbsTol(0);
        bc->SetAbsTol(1.0e-04);
        //bc->SetRelTol(1.0e-04);

        bc->SetCFactor(-1, 2);

        bc->SetResidual();

        int num = 4;
        double *nrm_arr = (double *) malloc(num * sizeof(double));
        bc->GetRNorms(&num, nrm_arr);

        bc->Drive();

        for (int i = 0; i < num; i++) {
            std::cout << i << ":\t" << nrm_arr[i] << std::endl;
        }

        free(bc);
        free(a);

        out->write_time_pvd(this->filename, *this->a->ux);
        return 0;
    }
};


#endif //UG_PLUGIN_XBRAIDFORUG4_PBRAID_H
