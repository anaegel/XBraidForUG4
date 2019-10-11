//
// Created by parnet on 12.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H

#include "../libs/xbraid/braid/braid.hpp"

#define TRACE_INDEX 1
#define TRACE_CONST 0
#define TRACE_DEFECT 1
#define TRACE_ACCESS 1

/**
 * To display the grid functions as matlab vectors (for debugging)
 */
#define TRACE_GRIDFUNCTION 0

/**
 * To display the time points to which the BufferUnpack and BufferPack functions are requested to
 */
#define TRACE_RECVTIME 1
#define TRACE_TIMINGS 1

#if TRACE_TIMINGS == 1
    #define StartRedoran(opt) this->redoran.get(opt).start()
    #define StopRedoran(opt) this->redoran.get(opt).stop()
    #define StartRedoranLevel(opt,l) this->redoran.get(opt,l).start()
    #define StopRedoranLevel(opt,l) this->redoran.get(opt,l).stop()
#else
    #define StartRedoran(opt)
    #define StopRedoran(opt)
    #define StartRedoranLevel(opt, l)
    #define StopRedoranLevel(opt, l)
#endif


#if TRACE_GRIDFUNCTION == 1
    #define MATLAB(u,i,t) this->matlab->write(u,i,t)
#else
    #define MATLAB(u, i, t)
#endif


/*
 * this structure represents the grid functions that are used within xbraid as vector types. It is possible to create
 * a class which stores additional information in the void pointer to use this structure within different context without
 * changing it.
 */
typedef struct _braid_Vector_struct {
    void *value;


#if TRACE_INDEX == 1
    size_t index = 0;
#endif

#if TRACE_CONST == 1
    bool m_const = true;
#endif
} BraidVector;

#if TRACE_INDEX == 1
size_t indexpool = 0;
#endif

#if TRACE_CONST == 1
size_t const_free = 0;
size_t const_clone = 0;
#endif

#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H
