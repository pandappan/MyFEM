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
#include <cassert>
#include <iomanip>
#include <fstream>
#include <vector>
#include "../Core/Types.h"

class Writer;
//!	Node class
class CNode
{
public:

//!	Maximum number of degrees of freedom per node
/*!	For 3D bar and solid elements, NDF = 3. For 3D beam or shell elements, NDF = 6 or 6 */
	const static unsigned int NDF = 3;

//!
	unsigned int Index = 0;

//!	x, y and z coordinates of the node
	double XYZ[3];

//! 节点自由度约束代码其中， 0自由 1固定约束, 2指定位移约束（默认位移为0）,3表示从属自由度
//! 只有自由自由度才形成方程号，其余为0
//! 初始默认全部采用固定约束，读入节点约束代码后设置实际约束情况
	unsigned int bcode[NDF];
	unsigned int eqn[NDF];   // 节点全局方程号

//! 节点位移值
	double Displacement[NDF];

//! 节点外载荷，点载，面载荷，体载贡献
	double NodeForce[NDF];

//! 节点约束力，边界条件贡献
	double NodeBCForce[NDF];

//! 抹平后的节点应力
	std::vector<double> stress;
//! 抹平时的累计权重
	double stressWieghts = 0.0;

//!	Constructor
	CNode(double X = 0.0, double Y = 0.0, double Z = 0.0);

//! 网格几何基础信息设置
// 设置几何信息
	void SetGeom(unsigned int nodeId_0, std::vector<double>& XYZ);
//! 维度约束信息设置
// 根据维度信息设置2D情况的约束信息
	void SetDimConstraints(unsigned int dim);
//! 固定约束信息设置
// 输入固定约束的自由度编号数组，0基，设置自由度为1
	void SetFixConstraints(const std::vector<unsigned int>& dofs_0);
//! 指定位移约束信息设置
// 输入指定位移约束的节点自由度，0基，设置自由度为2，并指定位移值
	void SetPreDispConstraints(unsigned int dof_0, double value);

//! 累加方式设置节点力
	inline void AddForce(unsigned int dof, double value) {
		NodeForce[dof] += value;
	}
//! 输出对应自由度的节点力
    inline double GetForce(unsigned int dof) const {
		return NodeForce[dof];
	}

//! 累加方式设置节点约束力
	inline void AddBcForce(unsigned int dof, double value) {
		NodeBCForce[dof] += value;
	}

//! 设置指定位移约束
	bool SetPreDisp(unsigned int dof, double value);

//!	Output nodal point data to stream
	template <class Stream>
	void Write(Stream& output, unsigned int dimension) const;

//!	Output equation numbers of nodal point to stream OutputFile
	template <class Stream>
	void WriteEquationNo(Stream& output, unsigned int dimension) const;

//!	Write nodal displacement
	template <class Stream>
	void WriteNodalDisplacement(Stream& output, unsigned int dimension) const;

//! Write node force
	template <class Stream>
	void WriteNodeForces(Stream& output, unsigned int dimension) const;

//! Write node BCforce
	template <class Stream>
	void WriteNodeBCForces(Stream& output, unsigned int dimension) const;

//! 将节点约束代码转换为全局方程号
	void GenerateNodeEquation(unsigned int& NEQ);

//! 将求解得到位移回代，更新节点自由度值
	void UpdataNodeDisplacement(const std::vector<double>& displacement);
};

//	Output nodal point data to stream
template <class Stream>
void CNode::Write(Stream& output, unsigned int dimension) const
{
	assert(dimension == 2 || dimension == 3);
	if (dimension == 2) {
		output << std::setw(6) << Index
		<< std::setw(6) << bcode[UX] << std::setw(6) << bcode[UY]
		<< std::setw(6) << XYZ[0] << std::setw(6) << XYZ[1]<< std::endl;
	} else {
		output << std::setw(6) << Index
		<< std::setw(6) << bcode[UX] << std::setw(6) << bcode[UY] << std::setw(6) << bcode[UZ]
	    << std::setw(6) << XYZ[0] << std::setw(6) << XYZ[1] << std::setw(6) << XYZ[2] << std::endl;
	}
}

//	Output equation numbers of nodal point to stream
template <class Stream>
void CNode::WriteEquationNo(Stream& output, unsigned int dimension) const
{
	assert(dimension == 2 || dimension == 3);
	output << std::setw(9) << Index << "       ";
    if (dimension == 2) {
	    output << std::setw(6) << eqn[UX] << std::setw(6) << eqn[UY];
    } else {
    	output << std::setw(6) << eqn[UX] << std::setw(6) << eqn[UY] << std::setw(6) << eqn[UZ];
    }
	output << std::endl;
}

//	Write nodal displacement
template <class Stream>
void CNode::WriteNodalDisplacement(Stream& output, unsigned int dimension) const
{
	assert(dimension == 2 || dimension == 3);
	output << std::setw(6) << Index << "        ";
    if (dimension == 2) { // 2D问题输出前2个自由度
		if (bcode[UX] == 0) { // 自由则直接输出位移
			output << std::setw(5) << Displacement[UX];
		} else { // 固定约束则输出0.0
			output << std::setw(6) << 0.0;
		}
    	if (bcode[UY] == 0) {
    		output << std::setw(6) << Displacement[UY];
    	} else {
    		output << std::setw(6) << 0.0;
    	}
    } else { // 3D问题输出6个自由度
    	for (unsigned int j = 0; j < NDF; j++)
    	{
    		if (bcode[j] == 0) // 自由则直接输出位移
    		{
    			output << std::setw(6) << Displacement[j];
    		}
    		else // 固定约束则输出0.0
    		{
    			output << std::setw(6) << 0.0;
    		}
    	}
    }
	output << std::endl;
}

template <class Stream>
void CNode::WriteNodeForces(Stream &output, unsigned int dimension) const
{
	assert(dimension == 2 || dimension == 3);
	output << std::setw(6) << Index << "        ";
	if (dimension == 2) {
		output << std::setw(6) << NodeForce[UX] << std::setw(6) << NodeForce[UY];
	} else {
		output << std::setw(6) << NodeForce[UX]
			   << std::setw(6) << NodeForce[UY]
			   << std::setw(6) << NodeForce[UZ];
	}
	output << std::endl;
}

template <class Stream>
void CNode::WriteNodeBCForces(Stream &output, unsigned int dimension) const
{
	assert(dimension == 2 || dimension == 3);
	output << std::setw(6) << Index << "        ";
	if (dimension == 2) {
		output << std::setw(6) << NodeBCForce[UX] << std::setw(6) << NodeBCForce[UY];
	} else {
		output << std::setw(6) << NodeBCForce[UX]
			   << std::setw(6) << NodeBCForce[UY]
			   << std::setw(6) << NodeBCForce[UZ];
	}
	output << std::endl;
}