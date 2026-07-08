//
// Created by Administrator on 2026/6/11.
//
#include <iostream>
#include <iomanip>
#include "CPlaneStressMaterial.h"
#include "../../Core/DenseMatrix.h"

bool CPlaneStressMaterial::Read(std::ifstream& Input) {
    Input >> nset >> rho >> E >> nu >> thk;
    return true;
}
void CPlaneStressMaterial::Write(std::ostream& Output) const {
    Output << std::setw(6) << rho
    << std::setw(6)<< E
    << std::setw(6) << nu
    << std::setw(6) << thk
    << std::endl;
}
void CPlaneStressMaterial::ComputeElasticMatrix(DenseMatrix<double>& D) const {
    double factor = E / (1.0 - nu * nu);
    D.SetZero();  // D 尺寸应为 3x3
    D(0,0) = factor;  D(0,1) = factor * nu;
    D(1,0) = factor * nu; D(1,1) = factor;
    D(2,2) = factor * (1 - nu) / 2;
}
