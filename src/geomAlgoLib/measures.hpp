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
    std::map<vertex_const_handle,Vertex_double_map> calculate_influence_map(const Polyhedron & myMesh,const Polyhedron & AABB_Box);
    CGAL::Bbox_3 box(const Polyhedron& mesh);
    std::array<Point3, 8> calculateBoundingBoxVertices(CGAL::Bbox_3 bbox);
    Polyhedron createMeshFromBoundingBoxVertices(const std::array<Point3, 8>& vertices);
    void translate_free_form(Polyhedron & myMesh, Kernel::Vector_3 translation, vertex_const_handle bounding_vertex, std::map<vertex_const_handle,Vertex_double_map> influence_map);
}