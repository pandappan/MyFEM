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
};
void CNode::SetGeom(unsigned int nodeId_0, std::vector<double>& xyz) {
	Index      = nodeId_0;
	XYZ[0] = xyz[0]; XYZ[1] = xyz[1]; XYZ[2] = xyz[2];
}
void CNode::SetFixConstraints(const std::vector<unsigned int>& dofs_0) {
	for (unsigned int d : dofs_0) {
		assert(d < NDF_MAX);
		bcode[d] = 1;
	}
}
void CNode::SetPreDispConstraints(unsigned int dof_0, double value) {
	assert(dof_0 < NDF_MAX);
	bcode[dof_0]        = 2;
	displacement[dof_0] = value;
}

// 将节点约束代码转化为全局方程号，并且返回当前全局最大方程号的引用
// 自由自由度bcode=0形成方程号，而约束自由度bcode!=0不形成方程号（固定，指定位移，从自由度）
// 只有激活的自由自由度才编号
void CNode::GenerateNodeEquation(unsigned int &NEQ) {
	for (unsigned int i = 0; i < NDF_MAX; i++) {
		if (IsDofActive(i) && bcode[i] == 0) {
			NEQ += 1;
			eqn[i] = NEQ;
		} else {
			eqn[i] = 0;
		}
	}
}

// 将总的位移根据自由度关系写入节点中，便于后续输出和形成约束力
// 只有激活的自由自由度才写入结果
void CNode::UpdateNodeDisplacement(const std::vector<double> &disp) {
	for (unsigned int i = 0; i < NDF_MAX; i++) {
		if (IsDofActive(i) && bcode[i] == 0) {
			displacement[i] = disp[eqn[i] - 1];
		}
	}
}