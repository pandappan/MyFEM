//
// Created by Administrator on 2026/6/10.
//
#pragma once
#include <iostream>
#include "Material.h"

class CSolid3DMaterial : public CMaterial {
public:
    bool Read(std::ifstream& Input) override;

    void Write(std::ostream& Output) const override;

    unsigned int GetNumStressComponents() const override {return 6;};

    void ComputeElasticMatrix(DenseMatrix<double>& D) const override;

};
