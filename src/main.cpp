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
    geomAlgoLib::Polyhedron newMesh = myMesh;
    auto genus = geomAlgoLib::computeGenus(myMesh);
    std::cout << "The Genus of [" << meshPath << "] is = " << std::to_string(genus) << std::endl;

    //geomAlgoLib::writeOFF(myMesh,"output_laplacien.off");
//TP1
    // geomAlgoLib::Facet_double_map tmp = geomAlgoLib::perimeter(myMesh);
    // geomAlgoLib::writeCOFF(myMesh,tmp,"output.off");
    
    // tmp = geomAlgoLib::norm_angle(myMesh);
    // geomAlgoLib::Facet_string_map tmpStringMap = geomAlgoLib::segment_mesh(tmp,45);
    
    // geomAlgoLib::writeNormOFF(myMesh,tmpStringMap,"outputnorm.off");
    // tmp = geomAlgoLib::angle(myMesh);
    // tmpStringMap = geomAlgoLib::segment_mesh(tmp,10);
    // geomAlgoLib::writeAngleOFF(myMesh,tmpStringMap,"outputangle.off");
    //std::cout << "la valeur max est " << geomAlgoLib::getMaxValue(tmp) << std::endl;

    //TP4
    geomAlgoLib::laplacien(myMesh,newMesh);
    geomAlgoLib::writeOFF(newMesh,"output_laplacien.off");
    std::cout << "The end..." << std::endl;
    return 0;
}