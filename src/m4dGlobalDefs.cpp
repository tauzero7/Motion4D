
#include "m4dGlobalDefs.h"

namespace m4d {

double radians(double phi) {
    return phi * DEG_TO_RAD;
}

double degree(double phi) {
    return phi * RAD_TO_DEG;
}


int find_nat_tetrad_type( const char* name) {
    unsigned int n = 0;
    while(n < NUM_ENUM_NAT_TETRAD_TYPES) {
        if (strcmp(name,stl_nat_tetrad_types[n])==0) {
            return n;
        }
        n++;
    }
    return -1;
}

} // end of namespace m4d
