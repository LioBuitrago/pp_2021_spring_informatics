// Copyright 2021 Galindo Buitrago Julio Cesar
#include <float.h>
#include <random>
#include <ctime>
#include <vector>
#include <string>
#include "../../../modules/task_3/galindo_j_SparseCSRconTBB/csrTBB.h"

DispersaMatrix transposicion(DispersaMatrix AT) {
    DispersaMatrix TrMat;
    std::vector<std::vector<std::complex<double>>> storageValue(AT.n);
    std::vector<std::vector<int>> storageCol(AT.n);
    TrMat.filaN.push_back(0);
    for (int i = 0; i < AT.n; i++) {
        for (int j = AT.filaN[i]; j < AT.filaN[i + 1]; j++) {
            storageValue[AT.columnaN[j]].push_back(AT.valor[j]);
            storageCol[AT.columnaN[j]].push_back(i);
        }
    }
    for (int i = 0; i < AT.n; i++)
        for (int j = 0; j < static_cast<int>(storageCol[i].size()); j++) {
            TrMat.columnaN.push_back(storageCol[i][j]);
            TrMat.valor.push_back(storageValue[i][j]);
        }
    int sum = 0;
    for (int i = 0; i < AT.n; i++) {
        sum += static_cast<int>(storageCol[i].size());
        TrMat.filaN.push_back(sum);
    }
    return TrMat;
}

DispersaMatrix multiplicacion(DispersaMatrix p, DispersaMatrix q) {
    DispersaMatrix r;
    if (p.n != q.n)
        throw(std::string)"incorrect size";
    r.n = p.n;
    std::vector<std::complex<double>> matC_val;
    std::vector<int> columnaN_en_R;
    std::vector<int> filaN_en_R;
    q = transposicion(q);
    filaN_en_R.push_back(0);
    for (int filaP = 0; filaP < p.n; filaP++) {
        for (int filaQ = 0; filaQ < p.n; filaQ++) {
            std::complex<double> escalar_s = { 0.0, 0.0 };
            for (int k = p.filaN[filaP]; k < p.filaN[filaP + 1]; k++) {
                for (int n = q.filaN[filaQ]; n < q.filaN[filaQ + 1]; n++) {
                    if (p.columnaN[k] == q.columnaN[n]) {
                        escalar_s += p.valor[k] * q.valor[n];
                        break;
                    }
                }
            }
            if (std::abs(escalar_s) >= DBL_MIN) {
                matC_val.push_back(escalar_s);
                columnaN_en_R.push_back(filaQ);
            }
        }
        filaN_en_R.push_back(matC_val.size());
    }
    for (int j = 0; j < static_cast<int>(columnaN_en_R.size()); j++) {
        r.columnaN.push_back(columnaN_en_R[j]);
        r.valor.push_back(matC_val[j]);
    }
    for (int i = 0; i <= p.n; i++) {
        r.filaN.push_back(filaN_en_R[i]);
    }
    return r;
}

DispersaMatrix multiplicacion_tbb(DispersaMatrix p, DispersaMatrix q) {
    DispersaMatrix r;
    if (p.n != q.n)
        throw(std::string)"incorrect size";
    r.n = p.n;
    tbb::task_scheduler_init init;
    std::vector<std::vector<std::complex<double>>> matC_val(p.n);
    std::vector<std::vector<int>> columnaN_en_R(p.n);
    std::vector<int> filaN_en_R;
    q = transposicion(q);
    filaN_en_R.push_back(0);

    tbb::parallel_for(tbb::blocked_range<int>(0, p.n, 50),
        [&](const tbb::blocked_range<int> & r) {
        int start = r.begin();
        int finish = r.end();
        std::vector<int>tmp(p.n);

        for (int filaP = start; filaP < finish; filaP++) {
            tmp.assign(p.n, -1);

            int indexfirst = p.filaN[filaP];
            int indexSecond = p.filaN[filaP + 1];

            for (int j = indexfirst; j < indexSecond; j++) {
                 int column = p.columnaN[j];
                 tmp[column] = j;
            }

            for (int j = 0; j < p.n; j++) {
                std::complex<double> escalar_s = { 0.0, 0.0 };
                for (int n = q.filaN[j]; n < q.filaN[j + 1]; n++) {
                    if (tmp[q.columnaN[n]] != -1) {
                    escalar_s += p.valor[tmp[q.columnaN[n]]] * q.valor[n];
                    }
                }

                if (std::abs(escalar_s) >= DBL_MIN) {
                    matC_val[filaP].push_back(escalar_s);
                    columnaN_en_R[filaP].push_back(j);
                }
            }
        }
    });
    init.terminate();
    for (int j = 0; j < static_cast<int>(p.n); j++) {
        for (int i = 0; i < static_cast<int>(columnaN_en_R[j].size()); i++) {
            r.columnaN.push_back(columnaN_en_R[j][i]);
            r.valor.push_back(matC_val[j][i]);
        }
        filaN_en_R.push_back(r.valor.size());
    }
    for (int i = 0; i <= p.n; i++) {
        r.filaN.push_back(filaN_en_R[i]);
    }
    return r;
}

DispersaMatrix ma_generacion(const std::vector<std::complex<double>> &vectorial, int nn) {
    DispersaMatrix matrix;
    matrix.n = nn;
    matrix.filaN.push_back(0);
    for (int i = 0; i < nn; i++) {
        for (int j = 0; j < nn; j++) {
            if (std::abs(vectorial[i*nn + j]) < DBL_MIN)
                continue;
            matrix.valor.push_back(vectorial[i*nn + j]);
            matrix.columnaN.push_back(j);
        }
        matrix.filaN.push_back(matrix.columnaN.size());
    }
    return matrix;
}

DispersaMatrix ma_generacion_dia(int nn) {
    DispersaMatrix resultMat;
    resultMat.n = nn;
    std::vector<std::complex<double>> values;
    std::vector<int> columnaN, filaN;
    filaN.push_back(0);

    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 0; i < nn; i++) {
        values.push_back(std::complex<double>(
            static_cast<double>(gen() % 100000),
            static_cast<double>(gen() % 100000)));
        columnaN.push_back(i);
        int a = i + 1;
        filaN.push_back(a);
    }

    resultMat.valor = values;
    resultMat.filaN = filaN;
    resultMat.columnaN = columnaN;

    return resultMat;
}
