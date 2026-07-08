//
// Created by Administrator on 2026/7/5.
//

#pragma once
#include <memory>
#include "../Core/Types.h"
class CElement;
class CMaterial;
std::unique_ptr<CElement> CreateElement(ElementTypes type);
std::unique_ptr<CMaterial> CreateMaterial(ElementTypes type);
