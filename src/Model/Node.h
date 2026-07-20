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
#include <array>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <fstream>
#include <vector>
#include "../Core/Types.h"

class Writer;
//!	Node class
class CNode
{
public:

//!
	unsigned int Index = 0;

//!	x, y and z coordinates of the node
	double XYZ[3];

	std::array<unsigned int, NDF_MAX> bcode {};
	std::array<unsigned int, NDF_MAX> eqn {};
	std::array<double, NDF_MAX> displacement {};
	std::array<double, NDF_MAX> nodeForce {};
	std::array<double, NDF_MAX> nodeBCForce {};

    std::uint16_t activeMask = 0;

	// 自由度激活
	inline void ActivateDof(unsigned int dof) {
		assert(dof < NDF_MAX);
		activeMask |= (1u << dof);
	}
	inline bool IsDofActive(unsigned int dof) const {
		return (activeMask & (1u << dof)) != 0;
	}
	inline bool HasAnyActiveDof() const {return activeMask != 0;}

	// 应力抹平
	std::vector<double> stress;
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
		nodeForce[dof] += value;
	}
//! 输出对应自由度的节点力
    inline double GetForce(unsigned int dof) const {
		return nodeForce[dof];
	}
//! 累加方式设置节点约束力
	inline void AddBcForce(unsigned int dof, double value) {
		nodeBCForce[dof] += value;
	}
//! 将节点约束代码转换为全局方程号
	void GenerateNodeEquation(unsigned int& NEQ);

//! 将求解得到位移回代，更新节点自由度值
	void UpdateNodeDisplacement(const std::vector<double>& disp);

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
			output << std::setw(5) << displacement[UX];
		} else { // 固定约束则输出0.0
			output << std::setw(6) << 0.0;
		}
    	if (bcode[UY] == 0) {
    		output << std::setw(6) << displacement[UY];
    	} else {
    		output << std::setw(6) << 0.0;
    	}
    } else { // 3D问题输出6个自由度
    	for (unsigned int j = 0; j < 3; j++)
    	{
    		if (bcode[j] == 0) // 自由则直接输出位移
    		{
    			output << std::setw(6) << displacement[j];
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
		output << std::setw(6) << nodeForce[UX] << std::setw(6) << nodeForce[UY];
	} else {
		output << std::setw(6) << nodeForce[UX]
			   << std::setw(6) << nodeForce[UY]
			   << std::setw(6) << nodeForce[UZ];
	}
	output << std::endl;
}

template <class Stream>
void CNode::WriteNodeBCForces(Stream &output, unsigned int dimension) const
{
	assert(dimension == 2 || dimension == 3);
	output << std::setw(6) << Index << "        ";
	if (dimension == 2) {
		output << std::setw(6) << nodeBCForce[UX] << std::setw(6) << nodeBCForce[UY];
	} else {
		output << std::setw(6) << nodeBCForce[UX]
			   << std::setw(6) << nodeBCForce[UY]
			   << std::setw(6) << nodeBCForce[UZ];
	}
	output << std::endl;
}