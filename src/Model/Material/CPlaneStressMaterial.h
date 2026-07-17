//
// Created by Administrator on 2026/6/10.
//

#pragma once

#include "Material.h"


class CPlaneStressMaterial : public CMaterial{
public:
    double thk;
    CPlaneStressMaterial();
    void Write(std::ostream &Output) const override;
    unsigned int GetNumStressComponents() const override {return 3;};
    void ComputeElasticMatrix(DenseMatrix<double>& D) const override;
    double GetThickness() const override {return thk;};
};
