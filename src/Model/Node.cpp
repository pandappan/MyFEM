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

	eqn[UX] = 0;	// 节点全局方程号，默认不进入方程中，后续编号时会确定
	eqn[UY] = 0;
	eqn[UZ] = 0;

	Displacement[UX] = 0.0;	// 节点自由度值，默认全部为0，后续施加指定位移约束或者求解结果回代会填充值
	Displacement[UY] = 0.0;
	Displacement[UZ] = 0.0;

	NodeForce[UX] = 0.0;
	NodeForce[UY] = 0.0;
	NodeForce[UZ] = 0.0;

};

//	Read element data from stream Input
bool CNode::Read(std::ifstream& Input, unsigned int dimension)
{
	assert(dimension == 2 || dimension == 3);
	Input >> NodeNumber;	// node number
	if (dimension == 2) {   // 读取2D的自由度约束代码，将3D剩余的自由度约束代码全部置为约束
		Input >> bcode[UX] >> bcode[UY]
		      >> XYZ[0] >> XYZ[1];
		bcode[UZ]    = 1;
		XYZ[2]       = 0.0;
	} else { // 读取3D自由度约束代码
		Input >> bcode[UX] >> bcode[UY] >> bcode[UZ]
	          >> XYZ[0] >> XYZ[1] >> XYZ[2];
	}

	return true;
}

// 将节点约束代码转化为全局方程号，并且返回全局前最大方程号
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

void CNode::UpdataNodeDisplacement(const std::vector<double> &displacement) {
	for (unsigned int i = 0; i < NDF; i++) {
		if (bcode[i] == 0) { // 只对自由度自由度更新节点自由度值
			Displacement[i] = displacement[eqn[i] - 1];
		}
	}
}
