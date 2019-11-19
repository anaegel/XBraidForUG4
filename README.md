# XBraidForUG4Lib

use ughub to install this plugin for ug4

xbraid should be compiled while the plugin compiles.
it could be necessary to set the compiler (for instance as a CC, CXX, MPICC, MPICXX variable) manual because an additional make is called during the cmake and this may not inherit the compiler from toolchains for ug4.

When using this plugin the SpaceTimeCommunicator has to be set so that the routines of UG4 and those of XBraid are working properly.
By doing this the Process Communication Layer will be overwritten, and therefore it is not possible to use -outproc for every proc, but only for those that are in the same spatial processor group.
Other feature like the profiling (or profiling pcl) might not work as expected.

It is possible to define variables in  BraidVectorStruct.h that are respected by the compiler to define the behavior of this plugin.
For example to
 * get execution times for some operations (clone, init, step, assemble operator, assemble rhs, solve, residual).
 * write out every result as a vector which can be used in matlab
 * write out a script in matlab notation
 * get timestamps when receiving and sending is executed
 
## Possible extensions:
Here are some ideas:
* Provide a functor class for converting BraidVectors to GridFunctions (cf. BraidVectorStruct.h for a draft).
* Make `GFBraidApp` a `TimeIntegratorSubject`. Also replace scriptors by `ITimeIntegratorObserver`?
* Can we use `DebugVectorWritingObject`?
* For systems (Navier-Stokes, Biot, ...), we should consider an alternative to `Interpolate` calls which does not depend on components. 

 
