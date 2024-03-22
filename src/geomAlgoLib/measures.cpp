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

    Polyhedron laplacien(Polyhedron & myMesh,Polyhedron & newMesh){

        double sumx,sumy,sumz,newx,newy,newz;
        auto new_vertex_iter = newMesh.vertices_begin();
        for (auto vertex_iter = myMesh.vertices_begin(); vertex_iter != myMesh.vertices_end(); ++vertex_iter) {
            
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
            //printf("Old coordinates : %f %f %f\n",new_vertex_iter->point().x(),new_vertex_iter->point().y(),new_vertex_iter->point().z());
            new_vertex_iter->point()=newPoint;
            //printf("New coordinates : %f %f %f\n",new_vertex_iter->point().x(),new_vertex_iter->point().y(),new_vertex_iter->point().z());
            ++new_vertex_iter;
        }
        return newMesh;
    }
    std::vector<Kernel::Vector_3> findNeighbors(const Polyhedron& mesh, const vertex_const_handle& vertex) {
        std::vector<Kernel::Vector_3> voisins;
        Halfedge_around_vertex_const_circulator cir = vertex->vertex_begin();
        Halfedge_around_vertex_const_circulator end = cir;

        if (!mesh.is_empty()) {
            do {
                
                vertex_const_handle neighbor_vertex = cir->opposite()->vertex();
                Kernel::Vector_3 neighbor_vector(neighbor_vertex->point().x(), neighbor_vertex->point().y(), neighbor_vertex->point().z());
                voisins.push_back(neighbor_vector);
                //printf("%d\n",cir);
            } while (++cir != end);
        }
        return voisins;
    }

    Polyhedron gaussien(Polyhedron & myMesh,Polyhedron & newMesh, double lambda){

        double sumx,sumy,sumz,newx,newy,newz,actualx,actualy,actualz,deltax,deltay,deltaz;
        auto new_vertex_iter = newMesh.vertices_begin();
        for (auto vertex_iter = myMesh.vertices_begin(); vertex_iter != myMesh.vertices_end(); ++vertex_iter) {
            actualx = vertex_iter->point().x();
            actualy = vertex_iter->point().y();
            actualz = vertex_iter->point().z();
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

            deltax = newx - actualx;
            deltay = newy - actualy;
            deltaz = newz - actualz;

            newx = actualx + lambda*deltax;
            newy = actualy + lambda*deltay;
            newz = actualz + lambda*deltaz;



            CGAL::Point_3<Kernel> newPoint(newx,newy,newz);
            //printf("Old coordinates : %f %f %f\n",new_vertex_iter->point().x(),new_vertex_iter->point().y(),new_vertex_iter->point().z());
            new_vertex_iter->point()=newPoint;
            //printf("New coordinates : %f %f %f\n",new_vertex_iter->point().x(),new_vertex_iter->point().y(),new_vertex_iter->point().z());
            ++new_vertex_iter;
        }
        return newMesh;
    }

    Polyhedron taubin(Polyhedron & myMesh,Polyhedron & newMesh, double lambda, double mu,int nb_ite){

        double sumx,sumy,sumz,newx,newy,newz,actualx,actualy,actualz,deltax,deltay,deltaz;
        for (int i = 0; i < nb_ite; i++){
            auto new_vertex_iter = newMesh.vertices_begin();
            for (auto vertex_iter = myMesh.vertices_begin(); vertex_iter != myMesh.vertices_end(); ++vertex_iter) {
                actualx = vertex_iter->point().x();
                actualy = vertex_iter->point().y();
                actualz = vertex_iter->point().z();
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

                deltax = newx - actualx;
                deltay = newy - actualy;
                deltaz = newz - actualz;

                if(!(i%2)){
                    newx = actualx + lambda*deltax;
                    newy = actualy + lambda*deltay;
                    newz = actualz + lambda*deltaz;
                }
                else{
                    newx = actualx + mu*deltax;
                    newy = actualy + mu*deltay;
                    newz = actualz + mu*deltaz;
                }
                
                CGAL::Point_3<Kernel> newPoint(newx,newy,newz);
                //printf("Old coordinates : %f %f %f\n",new_vertex_iter->point().x(),new_vertex_iter->point().y(),new_vertex_iter->point().z());
                new_vertex_iter->point()=newPoint;
                //printf("New coordinates : %f %f %f\n",new_vertex_iter->point().x(),new_vertex_iter->point().y(),new_vertex_iter->point().z());
                ++new_vertex_iter;
            }
            myMesh = newMesh;
        }
        
        return newMesh;
    }

    //FREE FORM

    std::map<int,Vertex_double_map> influence_map(const Polyhedron & myMesh,const Polyhedron & AABB_Box){
        std::map<int,Vertex_double_map> influence_maps;
        int map_id = 0;

        //1 ITERATION = 1 POINT DE LA BOX
        for (auto vertex_iter = AABB_Box.vertices_begin(); vertex_iter != AABB_Box.vertices_end(); ++vertex_iter) {
            Halfedge_around_vertex_const_circulator halfedge = vertex_iter->vertex_begin();

            //CALCUL DES VECTEURS INCIDENTS AU POINT DE LA BOX
            Kernel::Vector_3 vecs[3];
            int i =0;
            do{
                vertex_const_handle opposite_vertex = halfedge->opposite()->vertex();
                Kernel::Vector_3 vector(opposite_vertex->point().x()-vertex_iter->point().x(),opposite_vertex->point().y()-vertex_iter->point().y(),opposite_vertex->point().z()-vertex_iter->point().z());
                vecs[i] = vector;
                ++i;
                ++halfedge;
            }
            while(halfedge!=vertex_iter->vertex_begin());

            //PROJECTION DES POINTS DU MESH SUR CES VECTEURS
            for (auto vertex_iter2 = myMesh.vertices_begin(); vertex_iter2 != myMesh.vertices_end(); ++vertex_iter2) {
                float influence = 1;
                for (i = 0; i < 3; i++){
                    Point3 point = vertex_iter2->point();
                    CGAL::Line_3<Kernel> line(point,vecs[i]);
                    Point3 projected_point = line.projection(point);

                    //CALCUL DE LA DISTANCE ET DE L'INFLUENCE

                    influence *= sqrt(squared_distance(projected_point,vertex_iter->point()));
                }
            }
            map_id++;
        }
        return influence_maps;
    }
    CGAL::Bbox_3 box(const Polyhedron& mesh) {
        // Créez un itérateur pour parcourir les points du maillage
        auto begin = mesh.points_begin();
        auto end = mesh.points_end();

        // Utilisez la fonction bbox_3 pour calculer la boîte de délimitation des points du maillage
        return CGAL::bbox_3(begin, end);
    }

    std::array<Point3, 8> calculateBoundingBoxVertices(CGAL::Bbox_3 bbox) {
        std::array<Point3, 8> vertices;

        // Les sommets de la boîte de délimitation peuvent être calculés à partir de ses coordonnées minimales et maximales
        double xmin = bbox.xmin();
        double ymin = bbox.ymin();
        double zmin = bbox.zmin();
        double xmax = bbox.xmax();
        double ymax = bbox.ymax();
        double zmax = bbox.zmax();

        // Les huit sommets peuvent être calculés en combinant les différentes combinaisons des coordonnées minimales et maximales
        vertices[0] = Point3(xmin, ymin, zmin);
        vertices[1] = Point3(xmax, ymin, zmin);
        vertices[2] = Point3(xmin, ymax, zmin);
        vertices[3] = Point3(xmax, ymax, zmin);
        vertices[4] = Point3(xmin, ymin, zmax);
        vertices[5] = Point3(xmax, ymin, zmax);
        vertices[6] = Point3(xmin, ymax, zmax);
        vertices[7] = Point3(xmax, ymax, zmax);

        return vertices;
    }

}