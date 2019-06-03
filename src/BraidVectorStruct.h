//
// Created by parnet on 12.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H

#include "../libs/xbraid/braid/braid.hpp"

typedef struct _braid_Vector_struct {
    void *value;
    size_t index = 0;
    double time = 0.0;
} BraidVector;


size_t indexpool = 0;

#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H
