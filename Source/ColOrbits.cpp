
#include "matrix.h"
#ifndef USE_CUDA		// NOT yet implemented for GPU
template <>int CColOrbitCS<SIZE_TYPE>::m_maxElement = 0;
#endif




