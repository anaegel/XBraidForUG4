//
// Created by parnet on 12.05.19.
//

#ifndef UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H
#define UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H


// XBraid lib.
#include "../libs/xbraid/braid/braid.hpp"

// UG4 lib.
#include "common/util/smart_pointer.h"  		// for SmartPtr<T>

// Plugin lib.
#include "trace_tools_config.h"					// for TRACE_*

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

namespace ug {
namespace XBraidForUG4 {

//! This functors provide the interface for BraidVectors
template <typename TGridFunction>
class BraidVectorFunctors
{
public:
	// typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
	typedef SmartPtr<TGridFunction> SPGridFunction;

	static SPGridFunction as_smart_pointer(BraidVector &x)
	{ return *((SPGridFunction *) x.value); }

	static TGridFunction& as_grid_function(BraidVector &x)
	{ return *as_smart_pointer(x); }

	static const TGridFunction& as_grid_function(const BraidVector &x)
	{ return *as_smart_pointer(x); }

};

} // namespace XBraidForUG4
} // namespace ug

#endif //UG_PLUGIN_XBRAIDFORUG4_BRAIDVECTORSTRUCT_H
