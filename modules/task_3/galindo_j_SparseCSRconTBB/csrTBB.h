// Copyright 2021 Galindo Buitrago Julio Cesar
#ifndef MODULES_TASK_3_GALINDO_J_SPARSECSRCONTBB_CSRTBB_H_
#define MODULES_TASK_3_GALINDO_J_SPARSECSRCONTBB_CSRTBB_H_

#include <vector>
#include <complex>
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/tick_count.h"

struct DispersaMatrix {
    int n;
    std::vector<std::complex<double>> valor;
    std::vector<int> columnaN, filaN;
};

DispersaMatrix transposicion(DispersaMatrix AT);
DispersaMatrix multiplicacion(DispersaMatrix p, DispersaMatrix q);
DispersaMatrix multiplicacion_tbb(DispersaMatrix p, DispersaMatrix q);
DispersaMatrix ma_generacion(const std::vector<std::complex<double>> &vectorial, int nn);
DispersaMatrix ma_generacion_dia(int nn);

#endif  // MODULES_TASK_3_GALINDO_J_SPARSECSRCONTBB_CSRTBB_H_
