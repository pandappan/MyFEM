//
// Created by Administrator on 2026/6/11.
//
#include "Material.h"
#include "../../Core/DenseMatrix.h"

void CMaterial::ComputeStress(const std::vector<double>& strain, std::vector<double>& stress) const {
    unsigned int nComp = GetNumStressComponents();
    DenseMatrix<double> D(nComp, nComp);
    ComputeElasticMatrix(D);
    // 也可直接三重循环，这里用展开方式保证效率
    for (unsigned int i = 0; i < nComp; ++i) {
        stress[i] = 0.0;
        for (unsigned int j = 0; j < nComp; ++j) {
            stress[i] += D(i, j) * strain[j];
        }
    }
}
