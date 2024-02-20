#pragma once

#include "types.hpp"
#include "measures.hpp"

#include <iostream>
#include <fstream>

namespace geomAlgoLib
{

/// Read an OFF file and store the mesh in "mesh".
/// Returns true if it was loaded successfully.
bool readOFF(const std::string& filePath, Polyhedron& mesh);

/// Write a mesh at location "filePath"
void writeOFF(const Polyhedron& mesh, const std::string& filePath);

void writeCOFF(const Polyhedron& mesh, const Facet_double_map& faces, const std::string& filePath);

void writeNormOFF(const Polyhedron& mesh, const Facet_string_map& faces,const std::string& filePath);

void writeAngleOFF(const Polyhedron& mesh, const Facet_string_map& faces,const std::string& filePath);
}