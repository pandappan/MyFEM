//
// Created by Administrator on 2026/7/22.
//
#pragma once
#include "Q4.h"

class CQ4_AX : public CQ4{
protected:
    void ComputeBMatrix(unsigned int ip, DenseMatrix<double> &B) const override;
    double GetIntegrationVolumeFactor(unsigned int ip) const override;
public:
    std::string ElementTypeName() const override {return "CQ4_AX";}
    bool CalculateSurfaceLoad(unsigned int faceID, unsigned int dof, double value) override;
    MaterialCategory GetRequiredMaterial() const override;
};