//
// Created by parnet on 08.06.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_TALASMA_H
#define UG_PLUGIN_XBRAIDFORUG4_TALASMA_H

#include <chrono>
#include <ctime>
#include <ratio>

namespace TraceTools {

class Talasma {
public:
    std::chrono::high_resolution_clock::time_point t0;
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point t2;

    Talasma(){
        t0 = std::chrono::high_resolution_clock::now();
        t1 = std::chrono::high_resolution_clock::now();
        t2 = std::chrono::high_resolution_clock::now();
    }
    void start(){
        t0 = std::chrono::high_resolution_clock::now();

    }
    void stop(){
        t1 = std::chrono::high_resolution_clock::now();

    }

    void now(double &total,double &diff){
        t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> t0diff = std::chrono::duration_cast<std::chrono::duration<double> >(t2 - t0);
        total = t0diff.count();
        diff = 0;
    }
    double get(){
        std::chrono::duration<double> difference = std::chrono::duration_cast<std::chrono::duration<double> >(t1- t0);
        return difference.count();
    }
};

}

#endif //UG_PLUGIN_XBRAIDFORUG4_TALASMA_H
