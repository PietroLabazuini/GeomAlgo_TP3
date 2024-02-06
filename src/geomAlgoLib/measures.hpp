#pragma once
#include "types.hpp"

namespace geomAlgoLib
{
    Facet_double_map perimeter(const Polyhedron & myMesh);
    double getMaxValue(Facet_double_map map);
    double getMinValue(Facet_double_map map);

}