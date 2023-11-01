#ifndef STB_IMAGE_WRAPPER_H
#define STB_IMAGE_WRAPPER_H

// disable pedantic warnings for this external library
#ifdef _MSC_VER
    // Microsoft Visual C++ Compiler
    #pragma warning (push, 0)
#endif

#define STB_IMAGE_IMPLEMENTATION
#include "external/stb_image.hpp"

// restore warning levels
#ifdef _MSC_VER
    // Microsoft Visual C++ Compiler
    #pragma warning (pop)
#endif

#endif
