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
    // for(int i = 0; i<10; i++){
    //     geomAlgoLib::laplacien(myMesh,newMesh);
    //     myMesh = newMesh;
    // }
    // geomAlgoLib::writeOFF(newMesh,"output_laplacien.off");
    //for(int i = 0;i<10; i++){
    //    geomAlgoLib::gaussien(myMesh,newMesh, 0.33);
    //    myMesh = newMesh;
    //}
    //geomAlgoLib::writeOFF(newMesh,"output_gaussien.off");
    
    /*newMesh = geomAlgoLib::taubin(myMesh,newMesh,0.33,-0.34,100);
    geomAlgoLib::writeOFF(newMesh,"output_taubin.off");

    for(int i = 0;i<1000; i++){
        geomAlgoLib::gaussien(myMesh,newMesh, 0.33);
        myMesh = newMesh;
    }*/
    

    //std::cout << "The end..." << std::endl;


    std::array<geomAlgoLib::Point3, 8> vertices = geomAlgoLib::calculateBoundingBoxVertices(geomAlgoLib::box(myMesh));

    geomAlgoLib::Polyhedron Bounding_Mesh = geomAlgoLib::createMeshFromBoundingBoxVertices(vertices);

    std::map<geomAlgoLib::vertex_const_handle,geomAlgoLib::Vertex_double_map> influence_map = geomAlgoLib::calculate_influence_map(myMesh,Bounding_Mesh);

    CGAL::Vector_3<CGAL::Exact_predicates_inexact_constructions_kernel> vect(5,0,0);
    geomAlgoLib::translate_free_form(myMesh,vect,Bounding_Mesh.vertices_begin(),influence_map);

    geomAlgoLib::writeOFF(myMesh,"output_free_form.off");

    // Maintenant, vous pouvez accéder aux coordonnées des sommets comme suit
    /*for (int i = 0; i < 8; ++i) {
        std::cout << "Sommet " << i+1 << ": " << vertices[i] << std::endl;
    }*/

    return 0;
}