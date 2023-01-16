// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"
using namespace geometrycentral;
using namespace geometrycentral::surface;
using Real = double;
using Trip = Eigen::Triplet<size_t>;
/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    std::vector<Trip>tList;
    for(Halfedge he : mesh->halfedges())
        tList.emplace_back(geometry->edgeIndices[he.edge()], geometry->vertexIndices[he.tipVertex()], 1);
    SparseMatrix<size_t> SpM(mesh->nEdges(), mesh->nVertices());
    SpM.setFromTriplets(tList.begin(), tList.end());
    return SpM;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
    // TODO
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();
    std::vector<Trip>tList;
    for(Halfedge he : mesh->halfedges())
        tList.emplace_back(geometry->faceIndices[he.face()], geometry->edgeIndices[he.edge()], 1);
    SparseMatrix<size_t> SpM(mesh->nFaces(), mesh->nEdges());
    SpM.setFromTriplets(tList.begin(), tList.end());
    return SpM;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    geometry->requireVertexIndices();
    Vector<size_t> Vec(mesh->nVertices());
    int cnt = 0;
    for(Vertex v : mesh->vertices())
        Vec[cnt++] = subset.vertices.count(geometry->vertexIndices[v]) > 0 ? 1 : 0;
    return Vec;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {
    // TODO
    geometry->requireEdgeIndices();
    Vector<size_t> Vec(mesh->nEdges());
    int cnt = 0;
    for(Edge e : mesh->edges())
        Vec[cnt++] = subset.edges.count(geometry->edgeIndices[e]) > 0 ? 1 : 0;
    return Vec;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {
    // TODO
    geometry->requireFaceIndices();
    Vector<size_t> Vec(mesh->nFaces());
    int cnt = 0;
    for(Face f : mesh->faces())
        Vec[cnt++] = subset.faces.count(geometry->faceIndices[f]) > 0 ? 1 : 0;
    return Vec;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    // TODO
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();
    MeshSubset st = subset;
    for(size_t idx : subset.edges)
    {
        Edge e = mesh->edge(idx);
        auto adjF = e.adjacentFaces();
        for(Face f : adjF)
        {
            size_t idxF = geometry->faceIndices[f];
            if(!st.faces.count(idxF)) st.faces.insert(idxF);
        }
    }
    for(size_t idx : st.vertices)
    {
        Vertex v = mesh->vertex(idx);
        auto adjE = v.adjacentEdges();
        for(Edge e : adjE)
        {
            size_t idxE = geometry->edgeIndices[e];
            if(!st.edges.count(idxE)) st.edges.insert(idxE);
        }
        auto adjF = v.adjacentFaces();
        for(Face f : adjF)
        {
            size_t idxF = geometry->faceIndices[f];
            if(!st.faces.count(idxF)) st.faces.insert(idxF);
        }
    }

    return st;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    // TODO
    geometry->requireEdgeIndices();
    geometry->requireVertexIndices();
    MeshSubset st = subset;
    for(size_t idx : st.faces)
    {
        Face f = mesh->face(idx);
        auto adjE = f.adjacentEdges();
        for(Edge e : adjE)
        {
            size_t idxE = geometry->edgeIndices[e];
            if(!st.edges.count(idxE)) st.edges.insert(idxE);
        }
    }
    for(size_t idx : st.edges)
    {
        Edge e = mesh->edge(idx);
        if(!st.vertices.count(geometry->vertexIndices[e.firstVertex()]))
            st.vertices.insert(geometry->vertexIndices[e.firstVertex()]);
        if(!st.vertices.count(geometry->vertexIndices[e.secondVertex()]))
            st.vertices.insert(geometry->vertexIndices[e.secondVertex()]);
    }
    return st;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    // TODO
    MeshSubset st1 = closure(star(subset));
    MeshSubset st2 = star(closure(subset));
    for(size_t idx : st2.vertices)
        st1.vertices.erase(idx);
    for(size_t idx : st2.faces)
        st1.faces.erase(idx);
    for(size_t idx : st2.edges)
        st1.edges.erase(idx);
    return st1;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
    // TODO
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    for(size_t idx : subset.edges)
    {
        Edge e = mesh->edge(idx);
        if(!subset.vertices.count(geometry->vertexIndices[e.firstVertex()]) ||
            !subset.vertices.count(geometry->vertexIndices[e.secondVertex()]))
            return false;
    }
    for(size_t idx : subset.faces)
    {
        Face f = mesh->face(idx);
        auto Es = f.adjacentEdges();
        for(Edge e : Es)
            if (!subset.edges.count(geometry->edgeIndices[e]))
                return false;
    }
    return true;
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    // TODO
    if(!isComplex(subset)) return -1; // placeholder
    if(!subset.faces.empty())
    {
        std::set<Edge> Es;
        for(size_t idx : subset.edges)
            Es.insert(mesh->edge(idx));
        for(size_t idx : subset.faces)
        {
            Face f = mesh->face(idx);
            for(Edge e : f.adjacentEdges())
                if(Es.count(e))Es.erase(e);
        }
        if(!Es.empty()) return -1;
        for(size_t idx : subset.vertices)
        {
            Vertex v = mesh->vertex(idx);
            if(!v.degree()) return -1;
        }
        return 2;
    }
    if(!subset.edges.empty())
    {
        for(size_t idx : subset.vertices)
        {
            Vertex v = mesh->vertex(idx);
            if(!v.degree()) return -1;
        }
        return 1;
    }
    return 0;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    // TODO
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();
    MeshSubset cl = closure(subset);
    MeshSubset st;
    for(size_t idx : cl.vertices)
    {
        Vertex v = mesh->vertex(idx);
        size_t deg = 0;
        for(Edge e : v.adjacentEdges())
        {
            if(cl.edges.count(geometry->edgeIndices[e]))
            {
                deg++;
                if(deg > 1) break;
            }
        }
        if(deg == 1) st.addVertex(idx);
    }
    for(size_t idx : cl.edges)
    {
        Edge e = mesh->edge(idx);
        size_t deg = 0;
        for(Face f : e.adjacentFaces())
        {
            if(cl.faces.count(geometry->faceIndices[f]))
            {
                deg++;
                if(deg > 1) break;
            }
        }
        if(deg == 1) st.addEdge(idx);
    }
    return closure(st);
}