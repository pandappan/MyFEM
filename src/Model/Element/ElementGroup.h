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

#include <fstream>
#include <memory>
#include <vector>
#include "../../Core/Types.h"

class CNode;
class CElement;
class CMaterial;

//! Element group class
class CElementGroup {
private:
    ElementTypes type_ = ElementTypes::UNDEFINED;
    std::vector<std::unique_ptr<CElement>> elements_;
    std::vector<std::unique_ptr<CMaterial>> materials_;
public:
    CElementGroup();
    ~CElementGroup();
    // 禁止拷贝
    CElementGroup(const CElementGroup&) = delete;
    CElementGroup& operator=(const CElementGroup&) = delete;
    // 允许移动？
    CElementGroup(CElementGroup&&);
    CElementGroup& operator=(CElementGroup&&);
    bool Read(std::ifstream& Input, std::vector<CNode>& nodelist);
    CElement& GetElement(unsigned int i) { return *elements_[i]; }
    const CElement& GetElement(unsigned int i) const { return *elements_[i]; }
    CMaterial& GetMaterial(unsigned int i) { return *materials_[i]; }
    const CMaterial& GetMaterial(unsigned int i) const { return *materials_[i]; }
    ElementTypes GetElementType() const { return type_; }
    unsigned int GetNUME() const {return elements_.size();}
    unsigned int GetNUMAT() const {return materials_.size();}
    ElementTypes GetElementType() {return type_;}
};