#pragma once
#include "types.hpp"

namespace geomAlgoLib
{
    Facet_double_map perimeter(const Polyhedron & myMesh);
    double getMaxValue(Facet_double_map map);
    double getMinValue(Facet_double_map map);
    Facet_double_map angle(const Polyhedron & myMesh);
    Facet_double_map norm_angle(const Polyhedron & myMesh);
    Facet_string_map segment_mesh(const Facet_double_map& measures, double threshold);
    // std::vector<Kernel::Vector_3>
    Polyhedron laplacien(Polyhedron & myMesh,Polyhedron & newMesh);
    Polyhedron gaussien(Polyhedron & myMesh,Polyhedron & newMesh, double lambda);
    Polyhedron taubin(Polyhedron & myMesh,Polyhedron & newMesh, double lambda, double mu,int nb_ite);
    std::vector<Kernel::Vector_3> findNeighbors(const Polyhedron& mesh, const vertex_const_handle& vertex);
}