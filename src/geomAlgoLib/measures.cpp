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
    
    Facet_double_map angle(const Polyhedron & myMesh) {
        // Définir un map pour stocker les angles les plus petits par face
        Facet_double_map smallest_angles;

        // Parcourir chaque face du polyèdre
        for (Facet_iterator i = myMesh.facets_begin(); i != myMesh.facets_end(); ++i) {
            Halfedge_facet_circulator j = i->facet_begin();

            // Initialiser l'angle minimum avec une valeur élevée
            double min_angle = std::numeric_limits<double>::max();

            // Parcourir chaque arête de la face
            do {
                // Calculer l'angle entre les arêtes en utilisant la fonction CGAL::approximate_angle
                double angle = CGAL::approximate_angle(j->vertex()->point(),(j->next())->vertex()->point(),(j->next()->next())->vertex()->point());
                // Mettre à jour l'angle minimum si nécessaire
                if (angle < min_angle) {
                    min_angle = angle;
                }

                // Avancer à l'arête suivante
                ++j;
            } while (j != i->facet_begin());

            // Enregistrer l'angle minimum pour cette face
            smallest_angles[i] = min_angle;
            printf( "Angle minimal de la face : %f\n",min_angle);
        }

        // Retourner le map contenant les angles les plus petits par face
        return smallest_angles;
    }

    Facet_double_map norm_angle(const Polyhedron & myMesh){
        Facet_double_map tmp;
        Kernel::Vector_3 vectx(1,0,0);
        Kernel::Vector_3 vecty(0,1,0);
        Kernel::Vector_3 vectz(0,0,1);

        for (  Facet_iterator i = myMesh.facets_begin(); i != myMesh.facets_end(); ++i) {
            Halfedge_facet_circulator j = i->facet_begin();
            CGAL::Vector_3 norm_vect = CGAL::normal(j->vertex()->point(),(j->next())->vertex()->point(),(j->next()->next())->vertex()->point());
            double angle_x = CGAL::approximate_angle(norm_vect,vectx);
            double angle_y = CGAL::approximate_angle(norm_vect,vecty);
            double angle_z = CGAL::approximate_angle(norm_vect,vectz);
            tmp[i]=angle_y;
        }
        return tmp;
    }


    Facet_string_map segment_mesh(const Facet_double_map& measures, double threshold) {
        Facet_string_map segmentation;

        // Parcourir chaque face et vérifier si la mesure dépasse le seuil
        for (const auto& entry : measures) {
            if (entry.second < threshold) {
                // Associer une étiquette à la face
                segmentation[entry.first] = "InfSeuil";
            } else {
                // Face ne dépassant pas le seuil, pas besoin d'étiquette
                segmentation[entry.first] = "";
            }
        }

        return segmentation;
    }

    Polyhedron & laplacien(Polyhedron & myMesh,Polyhedron & newMesh){
        newMesh = myMesh;
        double sumx,sumy,sumz,newx,newy,newz;
        Vertex_iterator new_vertex_iter = newMesh.vertices_begin();
        for (Vertex_iterator vertex_iter = myMesh.vertices_begin(); vertex_iter != myMesh.vertices_end(); ++vertex_iter) {
            
            std::vector<Kernel::Vector_3> tmp = findNeighbors(myMesh,vertex_iter);
            
            sumx = 0;
            sumy = 0;
            sumz = 0;

            for (const Kernel::Vector_3& vector : tmp) {
                sumx += vector.x();
                sumy += vector.y();
                sumz += vector.z();
            }

            newx = sumx/tmp.size();
            newy = sumy/tmp.size();
            newz = sumz/tmp.size();

            CGAL::Point_3<Kernel> newPoint(newx,newy,newz);
            CGAL::Point_3<Kernel> oldPoint = new_vertex_iter->point();
            oldPoint = newPoint;

            ++new_vertex_iter;
        }
        return newMesh;
    }
}