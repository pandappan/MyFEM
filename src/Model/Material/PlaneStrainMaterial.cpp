//
// Created by Administrator on 2026/6/11.
//
#include <istream>
#include <iomanip>
#include "PlaneStrainMaterial.h"
#include "../../Core/DenseMatrix.h"

CPlaneStrainMaterial::CPlaneStrainMaterial() :thk(1.0) {
    matType = MaterialTypes::PE;
}

void CPlaneStrainMaterial::Write(std::ostream& output) const {
    output << std::setw(6) << rho
    << std::setw(6) << E
    << std::setw(6) << nu
    << std::setw(6) << thk
    << std::endl;
}

void CPlaneStrainMaterial::ComputeElasticMatrix(DenseMatrix<double>& D) const  {
    D.SetZero();
    double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    D(0,0) = factor * (1.0 - nu);   D(0,1) = factor * nu;
    D(1,0) = factor * nu;            D(1,1) = factor * (1.0 - nu);
    D(2,2) = factor * (1.0 - 2.0 * nu) / 2.0;
}

MaterialCategory CPlaneStrainMaterial::GetCateogory() const {
    return MaterialCategory::MechaincalPlaneStrain;
}
