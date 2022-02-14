#include "helper.h"

void mat_print(const double* mem, int rows, int cols, const char* title)
{
    arma::mat A(mem, rows, cols);
    int mr = (rows > 5 ? 5 : rows) - 1;
    int mc = (cols > 5 ? 5 : cols) - 1;
    A(arma::span(0, mr), arma::span(0, mc)).print(std::cout, title);
}