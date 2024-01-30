#include "measures.hpp"
#include <iostream>
#include <CGAL/Simple_cartesian.h>

namespace geomAlgoLib
{

    void perimeter(const Polyhedron & myMesh){
        for (  Facet_iterator i = myMesh.facets_begin(); i != myMesh.facets_end(); ++i) {
            Halfedge_facet_circulator j = i->facet_begin();
            // Facets in polyhedral surfaces are at least triangles.
            float sum=0;
            do {
                sum += sqrt(CGAL::squared_distance(j ->vertex()->point(),(j->next())->vertex()->point()));
                j--;
            } while (j != i->facet_begin());
            std::cout << "Perimetre de la face : "<< sum << std::endl;
        }
    }

}