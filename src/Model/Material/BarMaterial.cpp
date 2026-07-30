//
// Created by Administrator on 2026/6/11.
//
#include <fstream>
#include <iostream>
#include <iomanip>
#include "BarMaterial.h"
#include "../../Core/DenseMatrix.h"


BarMaterial::BarMaterial() {
    matType = MaterialCategory::Mechanical1D;
}

//!	Write material data to Stream
void BarMaterial::Write(std::ostream& output) const {
    output << std::setw(6) << rho
    << std::setw(6) << E
    << std::setw(6) << Area
    << std::endl;
}

//  应力分量个数
unsigned int BarMaterial::GetNumStressComponents() const {return 1;}

//  弹性矩阵
void BarMaterial::ComputeElasticMatrix(DenseMatrix<double>& D) const {
    D.SetZero();
    D(0,0) = E;
}

MaterialCategory BarMaterial::GetCategory() const {
    return MaterialCategory::Mechanical1D;
}
