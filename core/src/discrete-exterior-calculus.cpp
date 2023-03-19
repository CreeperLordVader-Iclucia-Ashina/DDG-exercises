// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const
{
    Eigen::SparseMatrix<double> spm;
    std::vector<Eigen::Triplet<double>>tList;
    spm.resize(mesh.nVertices(), mesh.nVertices());
    for(Vertex v : mesh.vertices())
    {
        size_t idx = vertexIndices[v];
        tList.emplace_back(idx, idx, barycentricDualArea(v));
    }
    spm.setFromTriplets(tList.begin(), tList.end());
    return spm; // placeholder
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const
{
    Eigen::SparseMatrix<double> spm;
    std::vector<Eigen::Triplet<double>> tList;
    spm.resize(mesh.nEdges(), mesh.nEdges());
    for(Edge e : mesh.edges())
    {
        size_t idx = edgeIndices[e];
        Halfedge he = e.halfedge();
        tList.emplace_back(idx, idx, 0.5 * (cotan(he) + cotan(he.sibling())));
    }
    spm.setFromTriplets(tList.begin(), tList.end());
    return spm; // placeholder
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    Eigen::SparseMatrix<double> spm;
    std::vector<Eigen::Triplet<double>> tList;
    spm.resize(mesh.nFaces(), mesh.nFaces());
    for(Face f : mesh.faces())
    {
        size_t idx = faceIndices[f];
        tList.emplace_back(idx, idx, 1.0 / faceArea(f));
    }
    spm.setFromTriplets(tList.begin(), tList.end());
    return spm; // placeholder
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {
    Eigen::SparseMatrix<double> spm;
    std::vector<Eigen::Triplet<double> >tList;
    spm.resize(mesh.nEdges(), mesh.nVertices());
    for(Edge e : mesh.edges())
    {
        size_t a = vertexIndices[e.firstVertex()];
        size_t b = vertexIndices[e.secondVertex()];
        uint idx = edgeIndices[e];
        tList.emplace_back(idx, a, -1.0);
        tList.emplace_back(idx, b, 1.0);
    }
    spm.setFromTriplets(tList.begin(), tList.end());
    return spm; // placeholder
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {
    Eigen::SparseMatrix<double> spm;
    spm.resize(mesh.nFaces(), mesh.nEdges());
    std::vector<Eigen::Triplet<double>> tList;
    for(Edge e : mesh.edges())
    {
        size_t idx = edgeIndices[e];
        Face fA = e.halfedge().face();
        Face fB = e.halfedge().sibling().face();
        size_t idxA = faceIndices[fA];
        size_t idxB = faceIndices[fB];
        tList.emplace_back(idxA, idx, 1.0);
        tList.emplace_back(idxB, idx, -1.0);
    }
    spm.setFromTriplets(tList.begin(), tList.end());
    return spm; // placeholder
}

} // namespace surface
} // namespace geometrycentral