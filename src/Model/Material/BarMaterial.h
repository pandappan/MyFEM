//
// Created by Administrator on 2026/6/10.
//

#pragma once

#include "Material.h"

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

    double Area;	//!< Sectional area of a bar element

public:

    //!	Read material data from stream Input
    bool Read(std::ifstream& Input) override;
    //!	Write material data to Stream
    void Write(std::ostream& output) const override;
    //  应力分量个数
    unsigned int GetNumStressComponents() const override;
    //  弹性矩阵
    void ComputeElasticMatrix(DenseMatrix<double>& D) const override;
};
