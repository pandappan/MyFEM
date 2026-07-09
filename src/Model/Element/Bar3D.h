/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once
#include <ostream>
#include <fstream>
#include "Element.h"

//! Bar element class
class CBar3D : public CElement
{
private:
	// 静态成员变量，属于类本身的特性，不单独属于某个类，静态常量可以在类内定义，静态常量数组需要在类外定义
	static const DOFIndex ActiveDOFs[3];
	static const unsigned int NumActiveDOFsPerNode;
public:
//!	Constructor
	CBar3D();
//!	Read element data from stream Input
	bool Read(std::ifstream& Input, CElementGroup& group , std::vector<CNode>& nodelist) override;
//!	Write element data to stream
	void Write(std::ostream& output) const override;
	void WriteElementStress(std::ostream& output) const override;
	void ElementStiffness(DenseMatrix<double>& K) override;
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
};
