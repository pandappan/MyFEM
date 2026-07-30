#pragma once
#include "Material.h"

class AxisymMaterial : public Material {
public:
    AxisymMaterial();
    void Write(std::ostream& output) const override;
    unsigned int GetNumStressComponents() const override { return 4; }
    void ComputeElasticMatrix(DenseMatrix<double>& D) const override;
    MaterialCategory GetCategory() const override {
        return MaterialCategory::MechanicalAxisym;
    }
    // 轴对称厚度概念不适用，返回 1
    double GetThickness() const override { return 1.0; }
};