#ifndef FHE_COMMONS_H
#define FHE_COMMONS_H

#include <NTL/LLL.h>
#include "BigFixP.h"

NTL::RR to_RR(const BigTorusRef& a);

NTL::RR to_RR(const BigFixPRef& a);

void to_torus(BigTorusRef reps, const NTL::RR& a);

void to_fixP(BigFixPRef reps, const NTL::RR& a);


#endif //FHE_COMMONS_H
