//
//  backtrace.hpp
//  RNAMake
//
//  Created by Joseph Yesselman on 6/16/16.
//  Copyright © 2016 Joseph Yesselman. All rights reserved.
//

#include <base/backtrace.h>

inline
std::string
base::demangle( std::string trace ) {
    
    std::string::size_type begin, end;
    
    // find the beginning and the end of the useful part of the trace
    
    // On a Mac, part to be demangled starts with underscore, and then is followed by "+" increment,
    //  separated by spaces from surrounding.
    begin = trace.find(" _") + 1;
    end   = trace.find(" +",begin);
    
    //How it looks for Linux, with parentheses around part to be demangled.
    // /home/rhiju/src/rosetta/main/source/cmake/build_release/libutility.so(_Z15print_backtracev+0x23) [0x7fb75e081a3]
    if ( begin == std::string::npos || end == std::string::npos ) {
        begin = trace.find("(_") + 1;
        end   = trace.find("+",begin);
    }
    
    // if begina and end were found, we'll go ahead and demangle
    if ( begin != std::string::npos && end != std::string::npos ) {
        std::string mangled_trace = trace.substr(begin, end - begin);
        size_t maxName = 1024;
        int demangleStatus;
        
        char* demangledName = (char*) malloc(maxName);
        if ( (demangledName = abi::__cxa_demangle(mangled_trace.c_str(), demangledName, &maxName,
                                                  &demangleStatus)) && demangleStatus == 0 ) {
            trace = trace.substr(0,begin) + demangledName + trace.substr(end ); // the demangled name is now in our trace string
        }
        free(demangledName);
    }
    return trace;
}

////////////////////////////////////////////////////////////////////
//
// See note above if this is causing your build to not compile.
//
// stolen directly from
//
// https://developer.apple.com/library/mac/documentation/Darwin/Reference/ManPages/man3/backtrace.3.html
//
// see also:
//
// http://man7.org/linux/man-pages/man3/backtrace.3.html#NOTES
//
//   -- rhiju, 2014.
////////////////////////////////////////////////////////////////////

void
base::print_backtrace() {
#if defined(_WIN32) || defined(_WIN64)
      HANDLE process = GetCurrentProcess();
      HANDLE thread = GetCurrentThread();
      
      CONTEXT context;
      memset(&context, 0, sizeof(CONTEXT));
      context.ContextFlags = CONTEXT_FULL;
      RtlCaptureContext(&context);
      
      SymInitialize(process, NULL, TRUE);
      
      DWORD image;
      STACKFRAME64 stackframe;
      ZeroMemory(&stackframe, sizeof(STACKFRAME64));
      
    #ifdef _M_IX86
      image = IMAGE_FILE_MACHINE_I386;
      stackframe.AddrPC.Offset = context.Eip;
      stackframe.AddrPC.Mode = AddrModeFlat;
      stackframe.AddrFrame.Offset = context.Ebp;
      stackframe.AddrFrame.Mode = AddrModeFlat;
      stackframe.AddrStack.Offset = context.Esp;
      stackframe.AddrStack.Mode = AddrModeFlat;
    #elif _M_X64
      image = IMAGE_FILE_MACHINE_AMD64;
      stackframe.AddrPC.Offset = context.Rip;
      stackfram.AddrPC.Mode = AddrModeFlat;
      stackframe.AddrFrame.Offset = context.Rsp;
      stackframe.AddrFrame.Mode = AddrModeFlat;
      stackframe.AddrStack.Offset = context.Rsp;
      stackframe.AddrStack.Mode = AddrModeFlat;
    #elif _M_IA64
      image = IMAGE_FILE_MACHINE_IA64;
      stackframe.AddrPC.Offset = context.StIIP;
      stackframe.AddrPC.Mode = AddrModeFlat;
      stackframe.AddrFrame.Offset = context.IntSp;
      stackframe.AddrFrame.Mode = AddrModeFlat;
      stackframe.AddrBStore.Offset = context.RsBSP;
      stackframe.AddrBStore.Mode = AddrModeFlat;
      stackframe.AddrStack.Offset = context.IntSp;
      stackframe.AddrStack.Mode = AddrModeFlat;
    #endif
    
      for (size_t i = 0; i < 25; i++) {
        
        BOOL result = StackWalk64(
          image, process, thread,
          &stackframe, &context, NULL, 
          SymFunctionTableAccess64, SymGetModuleBase64, NULL);
        
        if (!result) { break; }
        
        char buffer[sizeof(SYMBOL_INFO) + MAX_SYM_NAME * sizeof(TCHAR)];
        PSYMBOL_INFO symbol = (PSYMBOL_INFO)buffer;
        symbol->SizeOfStruct = sizeof(SYMBOL_INFO);
        symbol->MaxNameLen = MAX_SYM_NAME;
        
        DWORD64 displacement = 0;
        if (SymFromAddr(process, stackframe.AddrPC.Offset, &displacement, symbol)) {
          printf("[%i] %s\n", i, symbol->Name);
        } else {
          printf("[%i] ???\n", i);
        }
        
      }
      
      SymCleanup(process);
    
#else
    static int tried_throw = -1;
    
    try {
        // try once to re-throw currently active exception
        tried_throw++;
        if (tried_throw == 0) throw;
    }
    catch (const std::exception &e) {
        std::cerr << __FUNCTION__ << " caught unhandled exception. what(): "
        << e.what() << std::endl;
    }
    catch (...) {
        std::cerr << __FUNCTION__ << " caught unknown/unhandled exception."
        << std::endl;
    }
    
    size_t const callstack_size = 128;
    void* callstack[callstack_size];
    const int nMaxFrames = sizeof(callstack) / sizeof(callstack[0]);
    int i = 0;
    int frames = backtrace(callstack, nMaxFrames);
    std::cout << frames << std::endl;
    char** strs = backtrace_symbols(callstack, frames);
    //std::cerr << utility::CSI_Magenta; // set color of cerr to magenta
    for ( i = 3; i < frames; ++i ) {
        std::cerr << demangle( strs[i] ).c_str() << std::endl;
    }
    //std::cerr << utility::CSI_Reset; // reset color of cerr
    free(strs);
#endif
}
