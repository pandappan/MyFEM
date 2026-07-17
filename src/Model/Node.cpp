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
#include <cassert>
#include <vector>
#include "Node.h"
#include "../Core/Types.h"

CNode::CNode(double X, double Y, double Z)
{
    XYZ[0] = X;		// Coordinates of the node，创建时默认为0.0，读取节点数据时会被赋值
    XYZ[1] = Y;
    XYZ[2] = Z;
    
    bcode[UX] = 0;	// Boundary codes，创建默认全部自由，读取节点数据时会被约束
    bcode[UY] = 0;
    bcode[UZ] = 0;

	eqn[UX] = 0;	// 节点全局方程号，默认不进入方程中，后续编号时会确定具体方程号
	eqn[UY] = 0;
	eqn[UZ] = 0;

	Displacement[UX] = 0.0;	// 节点位移值，默认全部为0，施加指定位移约束或者求解结果回代会填充值
	Displacement[UY] = 0.0;
	Displacement[UZ] = 0.0;

	NodeForce[UX] = 0.0; // 节点力，点载荷，面载，线载荷经过转换形成节点力
	NodeForce[UY] = 0.0;
	NodeForce[UZ] = 0.0;

	NodeBCForce[UX] = 0.0; // 节点约束力
	NodeBCForce[UY] = 0.0;
	NodeBCForce[UZ] = 0.0;

};

void CNode::SetGeom(unsigned int nodeId_0, std::vector<double>& xyz) {
	Index      = nodeId_0;
	XYZ[0] = xyz[0]; XYZ[1] = xyz[1]; XYZ[2] = xyz[2];
}
void CNode::SetDimConstraints(unsigned int dim) {
	if (dim == 2) bcode[UZ] = 1;
}
void CNode::SetFixConstraints(const std::vector<unsigned int>& dofs_0) {
	for (unsigned int d : dofs_0) bcode[d] = 1;
}
void CNode::SetPreDispConstraints(unsigned int dof_0, double value) {
	bcode[dof_0]        = 2;
	Displacement[dof_0] = value;
}

// 将节点约束代码转化为全局方程号，并且返回当前全局最大方程号的引用
// 自由自由度bcode=0形成方程号，而约束自由度bcode!=0不形成方程号
void CNode::GenerateNodeEquation(unsigned int &NEQ) {
	for (unsigned int i = 0; i < NDF; i++) {
		if (bcode[i] == 0) {
			NEQ += 1;
			eqn[i] = NEQ;
		} else {
			eqn[i] = 0;
		}
	}
}

// 将总的位移根据自由度关系写入节点中，便于后续输出和形成约束力
void CNode::UpdataNodeDisplacement(const std::vector<double> &displacement) {
	for (unsigned int i = 0; i < NDF; i++) {
		if (bcode[i] == 0) { // 只对自由度自由度更新节点自由度值
			Displacement[i] = displacement[eqn[i] - 1];
		}
	}
}

// 对节点的预定义位移自由度设置值
bool CNode::SetPreDisp(unsigned int dof, double value) {
	// 判定本自由度约束情况是否为预定义位移
	if (bcode[dof] != 2) return false;
	// 设置值
	Displacement[dof] = value;
	return true;
}