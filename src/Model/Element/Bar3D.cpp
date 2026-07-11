/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>
#include <cmath>

#include "ElementGroup.h"
#include "Bar3D.h"
#include "../Node.h"
#include "../Material/BarMaterial.h"
#include "../../Core/DenseMatrix.h"

const DOFIndex CBar3D::ActiveDOFs[3] = {UX, UY, UZ};
const unsigned int CBar3D::NumActiveDOFsPerNode = 3;

//	Constructor
CBar3D::CBar3D()
{
	AllocateStorage(3,2,6);
}

//	Read element data from stream Input
bool CBar3D::Read(std::ifstream& Input,
				  CElementGroup& group,
				  std::vector<CNode>& nodelist) {
	unsigned int N1, N2, MSet;
	Input >> N1 >> N2 >> MSet;
	ElementMaterial_ = &group.GetMaterial(MSet - 1); // 获取基类指针
	nodes_[0] = &nodelist[N1 - 1];
	nodes_[1] = &nodelist[N2 - 1];
	return true;
}

//	Write element data to stream
void CBar3D::Write(std::ostream& output) const
{
	output << std::setw(6) << nodes_[0]->NodeNumber
	<< std::setw(6) << nodes_[1]->NodeNumber
	<< std::setw(6) << ElementMaterial_->nset << std::endl;
}

void CBar3D::WriteElementStress(std::ostream& out) const {
	auto* mat = static_cast<CBarMaterial*>(ElementMaterial_);
	double DX[3];
	double L2 = 0.0;
	for (unsigned int i = 0; i < 3; ++i) {
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 += DX[i] * DX[i];
	}
	double S[6];
	for (unsigned int i = 0; i < 3; ++i) {
		S[i]     = -DX[i] * mat->E / L2;
		S[i + 3] = -S[i];
	}
	double stress = 0.0;
	unsigned int idx = 0;
	for (unsigned int i = 0; i < NEN_; ++i)
		for (unsigned int d = 0; d < 3; ++d)
			stress += S[idx++] * nodes_[i]->Displacement[d];
	out << std::setw(16)  << ElementNumber_
		<< std::setw(16) << stress * mat->Area   // Force
		<< std::setw(16) << stress               // Stress
		<< std::endl;
}

//	Calculate element stiffness matrix
void CBar3D::ElementStiffness(DenseMatrix<double>& K)
{
	K.SetZero();
	// 双指针的nodes什么意思？
	double dx = nodes_[1]->XYZ[0] - nodes_[0]->XYZ[0];
	double dy = nodes_[1]->XYZ[1] - nodes_[0]->XYZ[1];
	double dz = nodes_[1]->XYZ[2] - nodes_[0]->XYZ[2];
	double l2 = dx * dx + dy * dy + dz * dz;
	double l = sqrt(l2);
	CBarMaterial* mat = static_cast<CBarMaterial*>(ElementMaterial_);
	double k = mat->E * mat->Area / (l * l2);
	// 刚度系数
	double cxx = k * dx * dx;
	double cyy = k * dy * dy;
	double czz = k * dz * dz;
	double cxy = k * dx * dy;
	double cxz = k * dx * dz;
	double cyz = k * dy * dz;
	// 节点1 的 XYZ 自由度索引 0,1,2；节点2 为 3,4,5
	// 左上块 (0-2, 0-2)
	K(0,0) = cxx;  K(0,1) = cxy;  K(0,2) = cxz;
	K(1,0) = cxy;  K(1,1) = cyy;  K(1,2) = cyz;
	K(2,0) = cxz;  K(2,1) = cyz;  K(2,2) = czz;

	// 右上块 (0-2, 3-5) = -左上
	K(0,3) = -cxx; K(0,4) = -cxy; K(0,5) = -cxz;
	K(1,3) = -cxy; K(1,4) = -cyy; K(1,5) = -cyz;
	K(2,3) = -cxz; K(2,4) = -cyz; K(2,5) = -czz;

	// 左下块 (3-5, 0-2) = -左上
	K(3,0) = -cxx; K(3,1) = -cxy; K(3,2) = -cxz;
	K(4,0) = -cxy; K(4,1) = -cyy; K(4,2) = -cyz;
	K(5,0) = -cxz; K(5,1) = -cyz; K(5,2) = -czz;

	// 右下块 (3-5, 3-5) = +左上
	K(3,3) = cxx;  K(3,4) = cxy;  K(3,5) = cxz;
	K(4,3) = cxy;  K(4,4) = cyy;  K(4,5) = cyz;
	K(5,3) = cxz;  K(5,4) = cyz;  K(5,5) = czz;
}

//	Calculate element stress 
double CBar3D::ElementStress() const
{
	CBarMaterial* mat = static_cast<CBarMaterial*>(ElementMaterial_);	// Pointer to material of the element
	double DX[3];
	double L2 = 0.0;
	for (unsigned int i = 0; i < 3; ++i) {
		DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
		L2 += DX[i] * DX[i];
	}
	double S[6];
	for (unsigned int i = 0; i < 3; ++i) {
		S[i]     = -DX[i] * mat->E / L2;
		S[i + 3] = -S[i];
	}
	double stress = 0.0;
	unsigned int idx = 0;
	for (unsigned int i = 0; i < NEN_; ++i)
		for (unsigned int d = 0; d < 3; ++d)
			stress += S[idx++] * nodes_[i]->Displacement[d];
	return stress;
}

void CBar3D::CalculateBodyForce(const double *bodyForce) {
	// 单元长度
	double dx = nodes_[1]->XYZ[0] - nodes_[0]->XYZ[0];
	double dy = nodes_[1]->XYZ[1] - nodes_[0]->XYZ[1];
	double dz = nodes_[1]->XYZ[2] - nodes_[0]->XYZ[2];
	double l2 = dx * dx + dy * dy + dz * dz;
	double len = sqrt(l2);
	// 单元截面参数
	CBarMaterial* mat = static_cast<CBarMaterial*>(ElementMaterial_);
	// 将体积转化为等效节点力，依次写入节点力中
	double factor = 0.5 * len * mat->Area * mat->rho;
	for (unsigned int d = 0; d < NDim_; d++) {
		double eqforce = bodyForce[d] * factor;
		nodes_[0]->AddForce(d,eqforce);
		nodes_[1]->AddForce(d,eqforce);
	}
}

void CBar3D::GetVisualizationNodes(DenseMatrix<double>& coords) const {
	coords.Resize(3, NEN_);
	for (unsigned int i = 0; i < NEN_; i++)
		for (unsigned int d = 0; d < 3; d++)
			coords(d, i) = nodes_[i]->XYZ[d];
}

void CBar3D::GetVisualizationDeformeNodes(DenseMatrix<double>& coords) const {
	coords.Resize(3, NEN_);
	for (unsigned int i = 0; i < NEN_; i++)
		for (unsigned int d = 0; d < 3; d++)
			coords(d, i) = nodes_[i]->XYZ[d] + nodes_[i]->Displacement[d];
}

double CBar3D::GetRepresentativeStress() const {
	return ElementStress();
}
