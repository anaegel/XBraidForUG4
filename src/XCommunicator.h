//
// Created by parnet on 16.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_XCOMMUNICATOR_H
#define UG_PLUGIN_XBRAIDFORUG4_XCOMMUNICATOR_H

#include <unistd.h>

#include "common/assert.h"
#include "pcl/pcl_comm_world.h"

#include "../libs/xbraid/braid/braid.hpp"

/**
 * \brief class for splitting the global communicator into a spatial and a temporal communicator
 *
 * This class splits the MPI_Comm of ug4 (named PCL_COMM_WORLD) into a temporal and a spatial communicator.
 * The ug4 MPI_Comm will then be replaced by the SPATIAL Communicator.
 */
class XCommunicator {

public:
    MPI_Comm GLOBAL = PCL_COMM_WORLD;
    MPI_Comm TEMPORAL = PCL_COMM_WORLD;
    MPI_Comm SPATIAL = PCL_COMM_WORLD;

    int globalsize = 1;
    int temporalsize = 1;
    int spatialsize = 1;

    bool verbose = true;

    XCommunicator() = default;

    ~XCommunicator() = default;

    void split(int numSpatialProcesses) { // nproc = x_procs * t_procs
        int world_size;
        MPI_Comm_size(PCL_COMM_WORLD, &world_size);
        GLOBAL = PCL_COMM_WORLD;

        UG_ASSERT(world_size % numSpatialProcesses == 0, "process_x * process_t != total_process");
        globalsize = world_size;
        spatialsize = numSpatialProcesses;
        temporalsize = world_size / numSpatialProcesses;

        BraidUtil bu = BraidUtil();
        if (verbose) {
            std::cout << "World size before splitting is:\t" << world_size << std::endl;
        }


        bu.SplitCommworld(&GLOBAL, numSpatialProcesses, &SPATIAL, &TEMPORAL);

        if (verbose) {
            MPI_Comm_size(GLOBAL, &world_size);
            std::cout << "World size after splitting is:\t" << world_size << std::endl;
            MPI_Comm_size(TEMPORAL, &world_size);
            std::cout << "... with temporal world size:\t" << world_size << std::endl;
            MPI_Comm_size(SPATIAL, &world_size);
            std::cout << "... and spatial world size:\t" << world_size << std::endl << std::endl;
        }

        PCL_COMM_WORLD = SPATIAL; // replaces ugs world communicator with the communicator for spatial
    }


    int getGlobalSize() {
        return globalsize;
    }

    int getTemporalSize() {
        return temporalsize;
    }

    int getSpatialSize() {
        return spatialsize;
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

    // todo delete
    int outputRank = 0;

    void setOutput(int globalRank) {
        outputRank = globalRank;
    }

    bool communication_on = true;

    void print(std::string v) {
        if (communication_on) {
            bool s = true;
            if (getGlobalRank() != outputRank) {
                const char *buffer = v.c_str();
                int size = v.size();
                MPI_Send(buffer, size, MPI_CHAR, outputRank, 0, this->GLOBAL);
                usleep(100000);
            } else if (s) {
                std::cout << outputRank << "\t" << v << std::endl;
            }
        } else {
            std::cout << outputRank << "\t" << v << std::endl;
        }
    }

    void processPrint() {
        if (communication_on) {
            char buffer[2048];
            memset(buffer, 0, sizeof(buffer));
            char exittoken[200] = "exit";

            int rem = globalsize - 1;
            std::cout << "globalsize: " << globalsize << "\t " << rem << std::endl;
            MPI_Status status;
/*        for (int i = 0; i < this->globalsize; i++) {
            std::cout << "Process " << i << std::endl;
            if (i == outputRank) {
                continue;
            }*/

            MPI_Recv(&buffer, 2048, MPI_CHAR, MPI_ANY_SOURCE, 0, this->GLOBAL, &status);
            while (rem != 0) {
                if (strcmp(buffer, exittoken) == 0) {
                    rem--;
                }
                std::cout << status.MPI_SOURCE << "\t" << buffer << std::endl;
                memset(buffer, 0, sizeof(buffer));
                MPI_Recv(&buffer, 2048, MPI_CHAR, MPI_ANY_SOURCE, 0, this->GLOBAL, &status);
            }
        }

    }

};


#endif //UG_PLUGIN_XBRAIDFORUG4_XCOMMUNICATOR_H
