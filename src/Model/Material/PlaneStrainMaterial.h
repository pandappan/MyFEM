//
// Created by Administrator on 2026/6/10.
//
#pragma once

#include <iostream>
#include <fstream>
#include "Material.h"
template<class T>
class DenseMatrix;

class CPlaneStrainMaterial : public CMaterial {
public:
    double thk;  // 厚度（通常取1.0，用于单位厚度问题）
    CPlaneStrainMaterial();
    void Write(std::ostream& output) const override;
    unsigned int GetNumStressComponents() const override {return 3;};
    void ComputeElasticMatrix(DenseMatrix<double>& D) const override;
    double GetThickness() const override {return thk;};
    MaterialCategory GetCateogory() const override;
};
