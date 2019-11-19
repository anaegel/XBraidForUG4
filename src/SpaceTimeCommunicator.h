//
// Created by parnet on 16.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_SPACETIMECOMMUNICATOR_H
#define UG_PLUGIN_XBRAIDFORUG4_SPACETIMECOMMUNICATOR_H

#include <unistd.h>

#include "common/assert.h"
#include "pcl/pcl_comm_world.h"

#include "../libs/xbraid/braid/braid.hpp"  // for BraidUtil

namespace ug {
namespace XBraidForUG4 {
/**
 * \brief class for splitting the global communicator into a spatial and a temporal communicator
 *
 * This class splits the MPI_Comm of ug4 (named PCL_COMM_WORLD) into a temporal and a spatial communicator.
 * The ug4 MPI_Comm will then be replaced by the SPATIAL Communicator.
 */
class SpaceTimeCommunicator {

public:
    MPI_Comm GLOBAL;
    MPI_Comm TEMPORAL;
    MPI_Comm SPATIAL;
protected:
    int m_globalsize = 1;
    int m_temporalsize = 1;
    int m_spatialsize = 1;

    bool m_verbose;

public:
    SpaceTimeCommunicator()
    : GLOBAL(PCL_COMM_WORLD), TEMPORAL(PCL_COMM_WORLD),  SPATIAL(PCL_COMM_WORLD), m_verbose(true)
    {}

    ~SpaceTimeCommunicator() = default;

    void split(int numSpatialProcesses) { // nproc = x_procs * t_procs
        int world_size;
        MPI_Comm_size(PCL_COMM_WORLD, &world_size);
        GLOBAL = PCL_COMM_WORLD;

        UG_ASSERT(world_size % numSpatialProcesses == 0, "process_x * process_t != total_process");
        m_globalsize = world_size;
        m_spatialsize = numSpatialProcesses;
        m_temporalsize = world_size / numSpatialProcesses;

        BraidUtil bu = BraidUtil();
        if (m_verbose) {
            std::cout << "World size before splitting is:\t" << world_size << std::endl;
        }


        bu.SplitCommworld(&GLOBAL, numSpatialProcesses, &SPATIAL, &TEMPORAL);

        if (m_verbose) {
            MPI_Comm_size(GLOBAL, &world_size);
            std::cout << "World size after splitting is:\t" << world_size << std::endl;
            MPI_Comm_size(TEMPORAL, &world_size);
            std::cout << "... with temporal world size:\t" << world_size << std::endl;
            MPI_Comm_size(SPATIAL, &world_size);
            std::cout << "... and spatial world size:\t" << world_size << std::endl << std::endl;
        }

        PCL_COMM_WORLD = SPATIAL; // replaces ugs world communicator with the communicator for spatial
    }

    void unsplit(){
        PCL_COMM_WORLD = GLOBAL;
        SPATIAL = PCL_COMM_WORLD;
        TEMPORAL = PCL_COMM_WORLD;
    }

    int getGlobalSize() {
        return m_globalsize;
    }

    int getTemporalSize() {
        return m_temporalsize;
    }

    int getSpatialSize() {
        return m_spatialsize;
    }

    int getTemporalRank() {
        int rank = 0;
        MPI_Comm_rank(TEMPORAL, &rank);
        return rank;
    }

    int getSpatialRank() {
        int rank = 0;
        MPI_Comm_rank(SPATIAL, &rank);
        return rank;
    }

    int getGlobalRank() {
        int rank = 0;
        MPI_Comm_rank(GLOBAL, &rank);
        return rank;
    }
};

} // namespace XBraidForUG4
} // namespace ug

#endif //UG_PLUGIN_XBRAIDFORUG4_SPACETIMECOMMUNICATOR_H
