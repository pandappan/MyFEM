/**
 *单元工厂
 *根据不同的单元类型，形成对应的单元指针类型
 **/
#pragma once
#include <memory>
#include "../Core/Types.h"
class Element;
class Material;
std::unique_ptr<Material> CreateMaterialByString(const std::string& type);
std::unique_ptr<Element> CreateElementByString(const std::string& type);
ElementTypes StringToElementType(const std::string& type);
