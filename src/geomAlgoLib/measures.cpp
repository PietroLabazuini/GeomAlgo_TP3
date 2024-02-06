#include "measures.hpp"
#include <iostream>
#include <CGAL/Simple_cartesian.h>

namespace geomAlgoLib
{

    Facet_double_map perimeter(const Polyhedron & myMesh){
        Facet_double_map tmp;
        for (  Facet_iterator i = myMesh.facets_begin(); i != myMesh.facets_end(); ++i) {
            Halfedge_facet_circulator j = i->facet_begin();
            // Facets in polyhedral surfaces are at least triangles.
            float sum=0;
            do {
                sum += sqrt(CGAL::squared_distance(j ->vertex()->point(),(j->next())->vertex()->point()));
                j--;
            } while (j != i->facet_begin());
            std::cout << "Perimetre de la face : "<< sum << std::endl;
            tmp[i]=sum;
        }
        return tmp;
    }

    double getMaxValue(Facet_double_map map){
        if (map.empty()) {
            return -1.0;
        }

        // Utilisation d'un itérateur pour parcourir la map
        auto maxElement = std::max_element(
            map.begin(), map.end(),
            [](const auto& p1, const auto& p2) {
                return p1.second < p2.second;
            }
        );

        // Retourne la valeur maximale
        return maxElement->second;
        }
    double getMinValue(Facet_double_map map){
        if (map.empty()) {
            return -1.0;
        }

        // Utilisation d'un itérateur pour parcourir la map
        auto minElement = std::min_element(
            map.begin(), map.end(),
            [](const auto& p1, const auto& p2) {
                return p1.second < p2.second;
            }
        );

        // Retourne la valeur maximale
        return minElement->second;
        }

}