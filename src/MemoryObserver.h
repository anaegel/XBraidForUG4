//
// Created by mparnet on 23.06.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_MEMORYOBSERVER_H
#define UG_PLUGIN_XBRAIDFORUG4_MEMORYOBSERVER_H


#include "sys/types.h"
#ifdef __APPLE__ && __MACH__
#include <sys/sysctl.h>
#else
#include <sys/sysinfo.h>
#endif

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "trace_tools.h"

// basing on https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
namespace TraceTools {

class MemoryObserver {
public:
    static long long getVirtualMemoryTotal() {
        struct sysinfo memInfo;
        sysinfo (&memInfo);
        long long totalVirtualMem = memInfo.totalram;
        totalVirtualMem += memInfo.totalswap;
        totalVirtualMem *= memInfo.mem_unit;
        return totalVirtualMem;
    }

    static long long getVirtualMemoryUsed() {
        struct sysinfo memInfo;
        sysinfo (&memInfo);
        long long virtualMemUsed = memInfo.totalram - memInfo.freeram;
        virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
        virtualMemUsed *= memInfo.mem_unit;
        return virtualMemUsed;
    }


    static int parseLine(char *line) {
        int i = strlen(line);
        const char *p = line;
        while (*p < '0' || *p > '9') p++;
        line[i - 3] = '\0';
        i = atoi(p);
        return i;
    }

    static int getVirtualMemoryConsumed() {
        FILE *file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != nullptr) {
            if (strncmp(line, "VmSize:", 7) == 0) {
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    }

    static long long getPhysicalMemoryTotal() {
        struct sysinfo memInfo;
        long long totalPhysMem = memInfo.totalram;
        totalPhysMem *= memInfo.mem_unit;
        return totalPhysMem;

    }

    static long long getPhysicalMemoryUsed() {
        struct sysinfo memInfo;
        long long physMemUsed = memInfo.totalram - memInfo.freeram;
        physMemUsed *= memInfo.mem_unit;
        return physMemUsed;
    }

    static int getPhysicalMemoryConsumed() {
        FILE *file = fopen("/proc/self/status", "r");
        int result = -1;
        char line[128];

        while (fgets(line, 128, file) != nullptr) {
            if (strncmp(line, "VmRSS:", 6) == 0) {
                result = parseLine(line);
                break;
            }
        }
        fclose(file);
        return result;
    }
};

}

#endif //UG_PLUGIN_XBRAIDFORUG4_MEMORYOBSERVER_H
