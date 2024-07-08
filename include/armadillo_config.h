#ifndef CONFIG_H
#define CONFIG_H

#ifdef USE_RCPPARMADILLO
#include <RcppArmadillo.h>
#else // USE_RCPPARMADILLO
#include <armadillo>
#endif // USE_RCPPARMADILLO

#endif // CONFIG_H
