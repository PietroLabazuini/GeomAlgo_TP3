#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Projection_traits_xy_3.h>

#include <map>

namespace geomAlgoLib{

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;
using Point3 = CGAL::Point_3<Kernel>;
using Projection_traits = CGAL::Projection_traits_xy_3<Kernel>;


using Facet_iterator = Polyhedron::Facet_const_iterator;
using Vertex_iterator = Polyhedron::Vertex_const_iterator;
using Halfedge_iterator = Polyhedron::Halfedge_const_iterator;
using Halfedge_facet_circulator = Polyhedron::Halfedge_around_facet_const_circulator;

using vertex_const_handle = Polyhedron::Vertex_const_handle;
using Halfedge_around_vertex_const_circulator = Polyhedron::Halfedge_around_vertex_const_circulator;

using Facet_double_map = std::map<Polyhedron::Facet_const_handle, double>;
using Facet_int_map = std::map<Polyhedron::Facet_const_handle, int>;
using Facet_string_map = std::map<Polyhedron::Facet_const_handle, std::string>;
using Vertex_double_map = std::map<vertex_const_handle, double>;

}