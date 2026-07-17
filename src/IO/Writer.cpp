//
// Created by Administrator on 2026/7/6.
//

#include "Writer.h"
#include "../Model/Model.h"
#include "../Model/Node.h"
#include "../Model/Element/Element.h"
#include "../Model/Element/ElementGroup.h"
#include "../Model/Material/Material.h"
#include "../Core/SkylineMatrix.h"
#include <ctime>
#include <cstdlib>

Writer::Writer(const std::string& fileName) : outputFile_(fileName) {
    if (!outputFile_.is_open()) {
        std::cerr << "Error opening output file." << std::endl;
        std::exit(3);
    }
}

void Writer::PrintTime() {
    time_t rawtime;
    time(&rawtime);
    tm* t = localtime(&rawtime);
    static const char* wd[] = {"Sunday","Monday","Tuesday","Wednesday",
                               "Thursday","Friday","Saturday"};
    static const char* mo[] = {"January","February","March","April","May","June",
                               "July","August","September","October","November","December"};
    *this << "        ("
          << t->tm_hour << ":" << t->tm_min << ":" << t->tm_sec
          << " on " << mo[t->tm_mon] << " " << t->tm_mday << ", "
          << t->tm_year + 1900 << ", " << wd[t->tm_wday] << ")\n\n";
}

void Writer::OutputHeading(const Model &model) {
    *this << "Title: " << model.title << std::endl;
    PrintTime();
}

void Writer::OutputNodeInfo(const Model& model) {
    *this << "C O N T R O L   I N F O R M A T I O N\n\n";
    *this << "  DIMENSION = " << model.dimension << "\n";
    *this << "  NUMNP     = " << model.nodes.size() << "\n";
    *this << "  NUMEG     = " << model.groups.size() << "\n\n";
    Tee([&](std::ostream& os) {
        for (auto& node : model.nodes) node.Write(os, model.dimension);
    });
    *this << "\n";
}

void Writer::OutputEquationNumber(const Model& model) {
    *this << " EQUATION NUMBERS\n\n"
          << "   NODE   Degrees of freedom\n";
    if (model.dimension == 2)
        *this << "     N     X    Y\n";
    else
        *this << "     N     X    Y    Z\n";
    Tee([&](std::ostream& os) {
        for (auto& n : model.nodes) n.WriteEquationNo(os, model.dimension);
    });
    *this << "\n";
}

void Writer::OutputElementInfo(const Model& model) {
    *this << " M A T E R I A L   D E F I N I T I O N\n";
    Tee([&](std::ostream& os) {
        for (std::size_t m = 0; m < model.materials.size(); ++m) {
            os << std::setw(5) << m + 1;
            model.materials[m]->Write(os);
        }
    });
    *this << "\n E L E M E N T   G R O U P   D A T A\n\n";
    unsigned int gi = 0;
    for (const auto& group : model.groups) {
        *this << " GROUP " << gi++ << "  TYPE = "
              << ElementTypeName(group.GetElementType())
              << "  NUME = " << group.GetNUME() << "\n";
        Tee([&](std::ostream& os) {
            for (std::size_t e = 0; e < group.GetNUME(); ++e) {
                os << std::setw(5) << e + 1;
                group.GetElement(e).Write(os);
            }
        });
        *this << "\n";
    }
}

void Writer::OutputTotalSystemData(const Model& model) {
    if (!model.K) return;
    *this << "\nTOTAL SYSTEM DATA\n"
          << "  NEQ = " << model.neq << "\n"
          << "  NWK = " << model.K->size() << "\n"
          << "  MK  = " << model.K->GetMaximumHalfBandwidth() << "\n";
}
