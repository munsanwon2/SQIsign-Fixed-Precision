#include <ec_params.h>

_Static_assert(FP_TORSION_POWER == TORSION_EVEN_POWER, "field and curve torsion mismatch");
// p+1 divided by the power of 2
const digit_t p_cofactor_for_2f[1] = {3};

