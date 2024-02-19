#include "io.hpp"

namespace geomAlgoLib
{

bool readOFF(const std::string& filePath, Polyhedron& mesh){

    std::ifstream input(filePath);

	if (!input || !(input >> mesh) || mesh.is_empty())
	{
		std::cerr << "Invalid .OFF file" << std::endl;
		return false;
	}

    return true;

}

void writeOFF(const Polyhedron& mesh, const std::string& filePath)
{
	std::ofstream in_myfile;
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color information
			  << mesh.size_of_vertices() << ' ' 
			  << mesh.size_of_facets() << " 0" << std::endl; 
			  // nb of vertices, faces and edges (the latter is optional, thus 0)

	std::copy(mesh.points_begin(), mesh.points_end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		Halfedge_facet_circulator j = i->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

		} while (++j != i->facet_begin());

		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Successfully exported at path: " << filePath << " !" << std::endl;
}


void writeNormOFF(const Polyhedron& mesh, const Facet_double_map& faces,const std::string& filePath)
{
	std::ofstream in_myfile;
	double val;
	double max = getMaxValue(faces);
	double min = getMinValue(faces);
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color information
			  << mesh.size_of_vertices() << ' ' 
			  << mesh.size_of_facets() << " 0" << std::endl; 
			  // nb of vertices, faces and edges (the latter is optional, thus 0)

	std::copy(mesh.points_begin(), mesh.points_end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		Halfedge_facet_circulator j = i->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

		} while (++j != i->facet_begin());

		
		Facet_double_map tmp = geomAlgoLib::norm_angle(mesh);
		if(tmp[i]<45.0){
			in_myfile << " 0.0 0.0 1.0";
		}
		else{
			in_myfile << " 0.0 0.0 0.0";
		}		

		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Successfully exported at path: " << filePath << " !" << std::endl;
}
	void writeCOFF(const Polyhedron& mesh, const Facet_double_map& faces,const std::string& filePath)
{
	std::ofstream in_myfile;
	double val;
	double max = getMaxValue(faces);
	double min = getMinValue(faces);
	in_myfile.open(filePath);

	CGAL::set_ascii_mode(in_myfile);

	in_myfile << "COFF" << std::endl // "COFF" makes the file support color information
			  << mesh.size_of_vertices() << ' ' 
			  << mesh.size_of_facets() << " 0" << std::endl; 
			  // nb of vertices, faces and edges (the latter is optional, thus 0)

	std::copy(mesh.points_begin(), mesh.points_end(),
			  std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

	for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
	{
		Halfedge_facet_circulator j = i->facet_begin();

		CGAL_assertion(CGAL::circulator_size(j) >= 3);

		in_myfile << CGAL::circulator_size(j) << ' ';
		do
		{
			in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

		} while (++j != i->facet_begin());

		Facet_double_map::const_iterator pos = faces.find(i);
		val = (pos->second - min) / (max - min);

		in_myfile << ' ' << val;
		if (val == 1 || val == 0){
			in_myfile << ".0 0.0 0.0";
		}
		else{
			in_myfile << " 0.0 0.0";
		} 

		in_myfile << std::endl;
	}

	in_myfile.close();

	std::cout << "Successfully exported at path: " << filePath << " !" << std::endl;
	}

	void writeAngleOFF(const Polyhedron& mesh, const Facet_double_map& faces,const std::string& filePath)
{
    std::ofstream in_myfile;
    
    in_myfile.open(filePath);

    CGAL::set_ascii_mode(in_myfile);

    in_myfile << "COFF" << std::endl // "COFF" makes the file support color information
              << mesh.size_of_vertices() << ' ' 
              << mesh.size_of_facets() << " 0" << std::endl; 
              // nb of vertices, faces and edges (the latter is optional, thus 0)

    std::copy(mesh.points_begin(), mesh.points_end(),
              std::ostream_iterator<Kernel::Point_3>(in_myfile, "\n"));

    for (Facet_iterator i = mesh.facets_begin(); i != mesh.facets_end(); ++i)
    {
        Halfedge_facet_circulator j = i->facet_begin();

        CGAL_assertion(CGAL::circulator_size(j) >= 3);

        in_myfile << CGAL::circulator_size(j) << ' ';
        do
        {
            in_myfile << ' ' << std::distance(mesh.vertices_begin(), j->vertex());

        } while (++j != i->facet_begin());

        Facet_double_map tmp = angle(mesh);
        
        if (tmp[i] < 15) {
            in_myfile << " 0.0 1.0 0.0";
        }
        else{
            in_myfile << " 0.0 0.0 0.0";
        }
        in_myfile << std::endl;
    }

    in_myfile.close();

    std::cout << "Successfully exported at path: " << filePath << " !" << std::endl;



}
}



