//
// Created by Administrator on 2026/6/10.
//
#pragma once
#include <iostream>
#include "Material.h"

class Solid3DMaterial : public Material {
public:
    Solid3DMaterial();
    void Write(std::ostream& Output) const override;

    unsigned int GetNumStressComponents() const override {return 6;};

    void ComputeElasticMatrix(DenseMatrix<double>& D) const override;
    MaterialCategory GetCategory() const override;
};
