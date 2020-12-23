#ifndef _XBRAID4UG__TRACE_TOOLS_H_
#define _XBRAID4UG__TRACE_TOOLS_H_


#define TRACE_INDEX 0 // 1
#define TRACE_CONST 0
#define TRACE_DEFECT 0 // 1
#define TRACE_ACCESS 0 // 1

/**
 * To display the grid functions as matlab vectors (for debugging)
 */
#define TRACE_GRIDFUNCTION 0

/**
 * To display the time points to which the BufferUnpack and BufferPack functions are requested to
 */
#define TRACE_RECVTIME 0 // 1
#define TRACE_TIMINGS 0 // 1

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




#endif
