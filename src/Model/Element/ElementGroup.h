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
#include "json.hpp"
#include "../../Core/Types.h"

class Node;
class Element;
class Material;
class Model;
using json = nlohmann::json;


//! Element group class
class CElementGroup {
private:
    ElementTypes type_ = ElementTypes::UNDEFINED;
    std::vector<std::unique_ptr<Element>> elements_;
public:
    CElementGroup();
    ~CElementGroup();
    // 禁止拷贝
    CElementGroup(const CElementGroup&) = delete;
    CElementGroup& operator=(const CElementGroup&) = delete;
    // 允许移动？
    CElementGroup(CElementGroup&&);
    CElementGroup& operator=(CElementGroup&&);
    // 初始化group信息
    void SetGroupsInfo(ElementTypes elemType, unsigned int num);
    // 加入group的单元
    void AddElement(std::unique_ptr<Element> elem);
    Element& GetElement(unsigned int i) { return *elements_[i]; }
    const Element& GetElement(unsigned int i) const { return *elements_[i]; }
    ElementTypes GetElementType() const { return type_; }
    unsigned int GetNUME() const {return elements_.size();}
    // 测试入口
    void SetTypeForTesting(ElementTypes t) { type_ = t; }
    void AddElementForTesting(std::unique_ptr<Element> elem);
};