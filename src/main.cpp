#include <io.hpp>
#include <example.hpp>
#include "measures.hpp"

#include <iostream>
#include <string>


int main(int argc, char *argv[]){

    std::cout << "Hello !" << std::endl;

    if(argc < 2){
        throw std::invalid_argument("This program expects at least 1 argument (path to a mesh).");
    }

    const std::string meshPath = std::string{argv[1]};
    
    geomAlgoLib::Polyhedron myMesh;

    geomAlgoLib::readOFF(meshPath, myMesh);

    auto genus = geomAlgoLib::computeGenus(myMesh);
    std::cout << "The Genus of [" << meshPath << "] is = " << std::to_string(genus) << std::endl;

    

    geomAlgoLib::Facet_double_map tmp = geomAlgoLib::perimeter(myMesh);
    geomAlgoLib::writeCOFF(myMesh,tmp,"output.off");

    std::cout << "la valeur max est " << geomAlgoLib::getMaxValue(tmp) << std::endl;

    std::cout << "The end..." << std::endl;
    return 0;
}