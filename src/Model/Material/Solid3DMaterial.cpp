//
// Created by Administrator on 2026/6/11.
//
#include <iostream>
#include <iomanip>
#include "Solid3DMaterial.h"
#include "../../Core/DenseMatrix.h"

CSolid3DMaterial::CSolid3DMaterial() {
    matType = MaterialTypes::SOLID;
}

void CSolid3DMaterial::Write(std::ostream& Output) const {
    Output << std::setw(6)  << rho
    << std::setw(6) << E
    << std::setw(6) << nu
    << std::endl;
}

void CSolid3DMaterial::ComputeElasticMatrix(DenseMatrix<double>& D) const  {
    D.SetZero();
    double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    double c1 = factor * (1.0 - nu);
    double c2 = factor * nu;
    double c3 = factor * (1.0 - 2.0 * nu) / 2.0;

    D(0,0)=c1; D(0,1)=c2; D(0,2)=c2;
    D(1,0)=c2; D(1,1)=c1; D(1,2)=c2;
    D(2,0)=c2; D(2,1)=c2; D(2,2)=c1;
    D(3,3)=c3; D(4,4)=c3; D(5,5)=c3;
}
