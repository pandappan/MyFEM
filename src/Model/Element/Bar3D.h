#pragma once
#include <ostream>
#include <fstream>
#include "Element.h"

//! Bar element class
class Bar3D : public Element
{
private:
	// 静态成员变量，属于类本身的特性，不单独属于某个类，静态常量可以在类内定义，静态常量数组需要在类外定义
	static const DOFIndex ActiveDOFs[3];
	static const unsigned int NumActiveDOFsPerNode;
public:
//!	Constructor
	Bar3D();
//!	Write element data to stream
	void Write(std::ostream& output) const override;
	void WriteElementStress(std::ostream& output) const override;
	void ElementStiffness(DenseMatrix<double>& K) const override;
	double ElementStress() const;
	unsigned int GetNumActiveDOFsPerNode() const override {return NumActiveDOFsPerNode;};
	const DOFIndex* GetActiveDOFs() const override {return ActiveDOFs;};
	void CalculateBodyForce(const double* bodyForce) override;
//! 输出单元类型
	std::string ElementTypeName() const override {return "CBar3D";};
//! 返回节点初始坐标
	virtual void GetVisualizationNodes(DenseMatrix<double>& coords) const;
//! 返回节点变形坐标
	virtual void GetVisualizationDeformeNodes(DenseMatrix<double>& deformeCoords) const;
//! 杆单元应力
	virtual double GetRepresentativeStress() const override;
	// 获取单元所述的材料类型
	MaterialCategory GetRequiredMaterial() const override;
};
