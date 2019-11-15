//
// Created by parnet on 12.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H

#include "../libs/xbraid/braid/braid.hpp"
#include "trace_tools.h"

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
