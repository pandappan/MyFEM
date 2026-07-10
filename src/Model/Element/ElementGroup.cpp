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

CElementGroup::CElementGroup() = default;
CElementGroup::~CElementGroup() = default;
CElementGroup::CElementGroup(CElementGroup&&) = default;
CElementGroup& CElementGroup::operator=(CElementGroup&&) = default;
//! Read element group data from stream Input
bool CElementGroup::Read(std::ifstream& Input, std::vector<CNode>& nodelist)
{
    unsigned int etype, NUME, NUMAT;
    Input >> etype >> NUME >> NUMAT;

    type_ = static_cast<ElementTypes>(etype);

    materials_.reserve(NUMAT);
    
//  Loop over for all material property sets in this element group
    for (unsigned int mset = 0; mset < NUMAT; mset++)
    {
        std::unique_ptr<CMaterial> mat = CreateMaterial(type_);
        if (!mat) {
            std::cerr << "Unkonw element type: " << etype << std::endl;
            return false;
        }
        if (!mat->Read(Input)) return false;
        if (mat->nset != mset + 1)
        {
            std::cerr << "*** Error *** Material sets must be inputted in order !" << std::endl
            << "    Expected set : " << mset + 1 << std::endl
            << "    Provided set : " << mat->nset << std::endl;
            return false;
        }
        materials_.push_back(std::move(mat));
    }

//  Read element data lines
    elements_.reserve(NUME);
    
//  Loop over for all elements in this element group
    for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        unsigned int N;
        Input >> N;    // element number
        if (N != Ele + 1)
        {
            std::cerr << "*** Error *** Elements must be inputted in order !" << std::endl
            << "    Expected element : " << Ele + 1 << std::endl
            << "    Provided element : " << N << std::endl;
            return false;
        }
        std::unique_ptr<CElement> element = CreateElement(type_);
        if (!element) {
            std::cerr << "Unkonw element type: " << etype << std::endl;
            return false;
        }
        element->SetElementType(type_);
        element->SetElementNumber(N);
        if (!element->Read(Input, *this, nodelist)) return false;
        elements_.push_back(std::move(element));
    }
    return true;
}

void CElementGroup::AddMaterialForTesting(std::unique_ptr<CMaterial> mat) {
    materials_.push_back(std::move(mat));
}

void CElementGroup::AddElementForTesting(std::unique_ptr<CElement> elem) {
    elements_.push_back(std::move(elem));
}
