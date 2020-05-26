#include "hitable.h"

hitable::~hitable() {}

const __m128 hitable::EPS4 = _mm_set_ps1(1e-4);
const __m128 hitable::MINUSEPS4 = _mm_set_ps1(-1e-4);
const __m128 hitable::ONE4 = _mm_set_ps1(1.0f);