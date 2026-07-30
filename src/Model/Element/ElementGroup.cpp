/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "ElementGroup.h"
#include "Element.h"
#include "../Node.h"
#include "../Material/Material.h"
#include "../ElementFactory.h"
#include <iostream>

#include "Model.h"

CElementGroup::CElementGroup() = default;
CElementGroup::~CElementGroup() = default;
CElementGroup::CElementGroup(CElementGroup&&) = default;
CElementGroup& CElementGroup::operator=(CElementGroup&&) = default;

void CElementGroup::SetGroupsInfo(ElementTypes elemType,unsigned int num) {
    type_ = elemType;
    elements_.clear();
    elements_.reserve(num);
}

void CElementGroup::AddElement(std::unique_ptr<Element> elem) {
    elements_.push_back(std::move(elem));
}

void CElementGroup::AddElementForTesting(std::unique_ptr<Element> elem) {
    elements_.push_back(std::move(elem));
}
