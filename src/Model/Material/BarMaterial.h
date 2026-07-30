#pragma once

#include "Material.h"

//!	Material class for bar element
class BarMaterial : public Material
{
public:

    double Area;	//!< Sectional area of a bar element
public:
    BarMaterial();
    //!	Write material data to Stream
    void Write(std::ostream& output) const override;
    //  应力分量个数
    unsigned int GetNumStressComponents() const override;
    //  弹性矩阵
    void ComputeElasticMatrix(DenseMatrix<double>& D) const override;
    MaterialCategory GetCategory() const override;
};
