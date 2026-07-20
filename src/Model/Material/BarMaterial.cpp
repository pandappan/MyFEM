//
// Created by Administrator on 2026/6/11.
//
#include <fstream>
#include <iostream>
#include <iomanip>
#include "BarMaterial.h"
#include "../../Core/DenseMatrix.h"


CBarMaterial::CBarMaterial() {
    matType = MaterialTypes::Bar;
}

//!	Write material data to Stream
void CBarMaterial::Write(std::ostream& output) const {
    output << std::setw(6) << rho
    << std::setw(6) << E
    << std::setw(6) << Area
    << std::endl;
}

//  应力分量个数
unsigned int CBarMaterial::GetNumStressComponents() const {return 1;}

//  弹性矩阵
void CBarMaterial::ComputeElasticMatrix(DenseMatrix<double>& D) const {
    D.SetZero();
    D(0,0) = E;
}

MaterialCategory CBarMaterial::GetCateogory() const {
    return MaterialCategory::Mechaincal1D;
}
