#pragma once
#include <fstream>
#include <vector>
#include "Types.h"

template<class T>
class DenseMatrix;


//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class Material
{
public:

	unsigned int nset;	//!< Number of set
	MaterialCategory matType;
	double E;  //!< Young's modulus
	double nu;
	double rho;

public:
	Material() : nset(0), matType(MaterialCategory::UNDEFINED), E(0.0), nu(0.0), rho(0.0){}
//! Virtual deconstructor
    virtual ~Material() = default;

//!	Write material data to Stream
    virtual void Write(std::ostream& output) const = 0;
	// 获取应力分量数目
	virtual unsigned int GetNumStressComponents() const = 0;
	// 计算弹性矩阵
	virtual void ComputeElasticMatrix(DenseMatrix<double>& D) const = 0;
	// 计算应力
	virtual void ComputeStress(const std::vector<double>& strain, std::vector<double>& stress) const;
	// 获取单元厚度
	inline virtual double GetThickness() const { return 1.0; } // 默认厚度1
	// 材料类型
	virtual MaterialCategory GetCategory() const = 0;
};

