//
// Created by Administrator on 2026/7/5.
//

#ifndef MYFEM_READER_H
#define MYFEM_READER_H

#include <fstream>
#include <string>
class Model;
class Reader {
public:
    bool Read(const std::string& filename, Model& model);
private:
    bool ReadHeader(std::ifstream& Input, Model& model);
    bool ReadNodes(std::ifstream& Input, Model& model);
    bool ReadGroups(std::ifstream& input, Model& model);
    bool ReadLoads(std::ifstream& Input, Model& model);
};


#endif //MYFEM_READER_H
