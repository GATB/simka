#include "SimkaAlgorithm.hpp"

// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html
// (last example)
// also, to reduce compilation time, I'm splitting it into several (8) files that will be compiled in parallel

template class SimkaAlgorithm <KSIZE_1>;
