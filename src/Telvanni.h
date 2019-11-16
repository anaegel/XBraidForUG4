//
// Created by parnet on 08.06.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_TELVANNI_H
#define UG_PLUGIN_XBRAIDFORUG4_TELVANNI_H

#include <chrono>
#include <ctime>
#include <ratio>

#include "trace_tools_config.h"

namespace TraceTools {

class Telvanni {
public:
    std::chrono::high_resolution_clock::time_point t0;
    std::chrono::high_resolution_clock::time_point t1;

    double time = 0;
    int usage = 0;

    Telvanni(){

    }
    void start(){
        t0 = std::chrono::high_resolution_clock::now();
    }
    void stop(){
        t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> difference = std::chrono::duration_cast<std::chrono::duration<double> >(t1- t0);
        time += difference.count();
        usage++;
    }

    double getTime(){
        return time;
    }

    double getUsage(){
        return usage;
    }
};

enum Observer { T_INIT,
    T_CLONE,
    T_FREE,
    T_ACCESS,
    T_SUM,
    T_SEND,
    T_RECV
};
std::string ObserverNames[] = {"init", "clone", "free", "access", "sum", "send", "recv"};

enum LevelObserver {
    T_STEP,
    T_RESIDUAL,
    T_ASSEMBLE_OP,
    T_ASSEMBLE_RHS,
    T_ADAPTIVE_TOL,
    T_SOLVE
};
int cLevelObserver = 6;
std::string LevelObserverNames[] = {"step", "residual", "assemble_op", "assemble_rhs", "adaptive_tol","solve"};

class Redoran {
public:

    std::vector<Telvanni> timer = std::vector<Telvanni>();
    std::vector<std::vector<Telvanni> > leveltimer = std::vector<std::vector<Telvanni> >();

    Redoran(){

    }

    Redoran(int maxlevel){
        for(int i = Observer::T_INIT; i != Observer::T_RECV;i++){
          //  timer.emplace_back(Telvanni());
            timer.push_back(Telvanni());
        }


        leveltimer = std::vector<std::vector<Telvanni> >(cLevelObserver,
                std::vector<Telvanni>(maxlevel, Telvanni()));
    }

    Telvanni& get(Observer o){
        //std::cout <<"Observer: " <<  o << std::endl;
        return this->timer[o];
    }

    Telvanni& get(LevelObserver o, int level){
        //std::cout <<"LObserver: " <<  o << "&" << level << "be"<< std::endl;
        Telvanni& v= this->leveltimer[o][level];
        //std::cout <<"LObserver: " <<  o << "&" << level << "af"<< std::endl;
        return v;
    };

};

}
#endif //UG_PLUGIN_XBRAIDFORUG4_TELVANNI_H
