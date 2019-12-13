//
// Created by parnet on 12.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_PBRAID_H
#define UG_PLUGIN_XBRAIDFORUG4_PBRAID_H

// UG4 lib.
#include "common/assert.h"
#include "pcl/pcl_comm_world.h"
#include "lib_disc/time_disc/theta_time_step.h"
#include "lib_disc/time_disc/solution_time_series.h"
#include "bindings/lua/lua_function_handle.h"
#include "bindings/lua/lua_user_data.h"

// This lib.
#include "PBraidApp.h"
#include "RGFBraidApp.h"
#include "MGFBraidApp.h"
#include "ITSGFBraidApp.h"
#include "SpaceTimeCommunicator.h"
#include "Scriptor.h"

template<typename TDomain, typename TAlgebra>
class PBraid {
public:
    typedef ug::ThetaTimeStep <TAlgebra> TTimeStep;
    typedef ug::ILinearOperatorInverse<typename TAlgebra::vector_type> TSolver;
    typedef typename TAlgebra::vector_type vector_type;

    typedef ug::GridFunction <TDomain, TAlgebra> TGridFunction;
    typedef SmartPtr <TGridFunction> SPGridFunction;

    typedef SmartPtr<ug::UserData<double, TGridFunction::dim> > SPData;
    typedef RGFBraidApp<TDomain, TAlgebra> TBraidApp;

    typedef SmartPtr<UG4BraidApp<TDomain, TAlgebra> > SPBraidApp;

    typedef SmartPtr<ug::DomainDiscretization<TDomain, TAlgebra> > SPDomainDisc;

    typedef typename ug::XBraidForUG4::SpaceTimeCommunicator TSpaceTimeCommunicator;
    typedef SmartPtr<TSpaceTimeCommunicator> SPXCommunicator;
    typedef SmartPtr<Scriptor<TDomain, TAlgebra> > SPScriptor;

private:
    BraidCore *m_bc = nullptr;
    SPBraidApp m_app;
    SmartPtr <TSpaceTimeCommunicator> m_comm;
    const char *m_filename;

public:


    // -----------------------------------------------------------------------------------------------------------------
    // constructor and destructor
    // -----------------------------------------------------------------------------------------------------------------
    PBraid(SPXCommunicator p_comm, SPBraidApp p_app) {
        this->m_comm = p_comm;
        this->setApp(p_app);
    }

    ~PBraid() {
        if (m_bc != nullptr) {
            delete m_bc;
        }
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------------------------------------------------
    inline void setApp(SPBraidApp p_app) {
        if (m_bc != nullptr) {
            delete m_bc;
        }
        this->m_app = p_app;
        this->m_app->m_comm = this->m_comm;
        this->m_app->comm_t = this->m_comm->TEMPORAL;
        m_bc = new BraidCore(this->m_comm->GLOBAL, this->m_app.get());
    }

    SPBraidApp getApp() {
        return this->m_app;
    }


    // -----------------------------------------------------------------------------------------------------------------
    //
    // -----------------------------------------------------------------------------------------------------------------

    void createAccess(SmartPtr<ug::VTKOutput<TDomain::dim> > out) {
        // this->m_out = out;-
        // mgf->m_out = this->m_out;
    }

    void setFilename(const char *filename) {
        // mgf->m_filename = filename;
        // this->m_filename = filename;
    }


    void setMaxIterations(int maxIter) {
        this->m_bc->SetMaxIter(maxIter);
    }

    void setResidual() {
        this->m_bc->SetResidual();
    }

    void setNRelax(int level, int number) {
        this->m_bc->SetNRelax(level, number);
    }

    void setCFactor(int level, int factor) {
        this->m_bc->SetCFactor(level, factor);
    }

    void setMaxLevels(int maxLevel) {
        this->m_bc->SetMaxLevels(maxLevel); // todo check app max level?
    }

    void setSkip(bool skip) {
        this->m_bc->SetSkip(skip);

    }


    void setMinCoarse(int minCoarse) {
        this->m_bc->SetMinCoarse(minCoarse);
    }

    void setAbsoluteTol(double tol) {
        this->m_bc->SetAbsTol(tol);
    }

    void setReduction(double tol) {
        this->m_bc->SetRelTol(tol);
    }

    void setTemporalNorm(int nrm) {
        this->m_bc->SetTemporalNorm(nrm);
    }

    void setSequential() {
        this->m_bc->SetSeqSoln(1);
    }

    void setStoreValues(int level) {
        this->m_bc->SetStorage(level);
    }

    void setSpatialCoarsenAndRefine() {
        this->m_bc->SetSpatialCoarsenAndRefine();
    }

    void setRefine() {
        this->m_bc->SetRefine(1);
    }

    void setMaxRefinements(int number) {
        this->m_bc->SetMaxRefinements(number);
    }


    void setAccessLevel(int level) {
        this->m_bc->SetAccessLevel(level);
    }


    void setPrintLevel(int level) {
        this->m_bc->SetPrintLevel(level);
    }

    void setPrintFile(const char *file) {
        this->m_bc->SetPrintFile(file);
    }

    void setCycleFMG() {
        this->m_bc->SetFMG();
    }

    void setCycleNFMG(int vu) {
        this->m_bc->SetNFMG(vu);
    }

    void setCycleNFMGV(int mu) {
        this->m_bc->SetNFMGVcyc(mu);
    }

    void setCycleType(const char *ctyp) {
        if (strcmp(ctyp, "V FCF") == 0) {
            this->m_bc->SetNRelax(-1, 1);
            this->m_bc->SetNRelax(0, 1);

        } else if (strcmp(ctyp, "V F-FCF") == 0) {
            this->m_bc->SetNRelax(-1, 1);
            this->m_bc->SetNRelax(0, 0);

        } else if (strcmp(ctyp, "F F") == 0) {
            this->m_bc->SetNRelax(-1, 0);
            this->m_bc->SetFMG();

        } else if (strcmp(ctyp, "F F-FCF") == 0) {
            this->m_bc->SetNRelax(-1, 1);
            this->m_bc->SetNRelax(0, 0);
            this->m_bc->SetFMG();
        }
    }


    int run() {
        this->m_app->init();
        this->m_bc->SetAccessLevel(2);

        this->m_bc->SetPrintLevel(3);
        this->m_bc->Drive();

        // GetRNorms
        // GetNLevels
        // GetNumIter
        return 0;
    }

#ifdef UG_FOR_LUA

    int test(SPGridFunction u0,
            const char *generator,
            const char *cmp,
            SPScriptor output){
        this->m_app->setStartVector(u0);
        this->m_app->setVectorGenerator(generator);
        this->m_app->setGeneratorComponent(cmp);
        this->m_app->setScriptor(output);
        this->m_app->init();

         FILE *file;
        file = fopen("myfile", "w");

         BraidUtil bu = BraidUtil();
         bu.TestResidual(this->m_app.get(),
                          this->m_comm->SPATIAL,
                          file,
                          0,
                          0.1);

        this->m_app->release();

    }

    int run(SPGridFunction u0,
            const char *generator,
            const char *cmp,
            SPScriptor output) {
        std::cout << "PBraid run" << std::endl;
        this->m_app->setStartVector(u0);
        this->m_app->setVectorGenerator(generator);
        this->m_app->setGeneratorComponent(cmp);
        this->m_app->setScriptor(output);

        UG_LOG("Init...")
        this->m_app->init();

        UG_LOG("Drive...")
        this->m_bc->Drive();

        UG_LOG("Release...")
        this->m_app->release();
        return 0;
    }

#endif
};


#endif //UG_PLUGIN_XBRAIDFORUG4_PBRAID_H
