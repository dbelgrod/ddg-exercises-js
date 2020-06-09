"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		// TODO
		// var dualArea = 0;
		// var T = new Triplet;
		// for (var key in vertexIndex)
		// {
		// 	for (var f of key.adjacentFaces())
		// 	{
		// 		for (var h of f.adjacentHalfedges())
		// 		{
		// 			var v = geometry.vector(h);
		// 			dualArea+= v.dot(v) * geometry.cotan(h);
		// 		}
		// 	}
		// 	T.push_back(vertexIndex[key], vertexIndex[key], dualArea/8);
		// }
		var size = Object.keys(vertexIndex).length;
		var T = new Triplet(size,size);
		for (var key in vertexIndex)
		{	
			var pos = vertexIndex[key];
			var vertex = geometry.mesh.vertices[pos];
			T.addEntry(geometry.barycentricDualArea(vertex), pos, pos);
		}
		return SparseMatrix.fromTriplet(T); 
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		// TODO
		var size = Object.keys(edgeIndex).length;
		var T = new Triplet(size,size);

		for (var key in edgeIndex)
		{
			var pos = edgeIndex[key];
			var edge = geometry.mesh.edges[pos];

			var ratio = geometry.cotan(edge.halfedge) + geometry.cotan(edge.halfedge.twin);
			ratio *= 0.5;
			T.addEntry(ratio, pos, pos);
		}
		return SparseMatrix.fromTriplet(T);
		//return SparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		// TODO
		var size = Object.keys(faceIndex).length;
		var T = new Triplet(size,size);

		for (var key in faceIndex)
		{
			var pos = faceIndex[key];
			var face = geometry.mesh.faces[pos];

			T.addEntry(1 / geometry.area(face), pos, pos );
		}
		return SparseMatrix.fromTriplet(T);
		//return SparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		// TODO

		var E = Object.keys(edgeIndex).length;
		var V = Object.keys(vertexIndex).length;
		var T = new Triplet(E,V);

		for (var key in edgeIndex)
		{
			var edge = geometry.mesh.edges[edgeIndex[key]];
			var e = edgeIndex[key];
			T.addEntry(1, e, edge.halfedge.next.vertex.index);
			T.addEntry(-1, e, edge.halfedge.vertex.index);
		}
		//return SparseMatrix.identity(1, 1); // placeholder
		return SparseMatrix.fromTriplet(T);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
		// TODO
		var E = Object.keys(edgeIndex).length;
		var F = Object.keys(faceIndex).length;
		var T = new Triplet(F,E);

		
		var seen = new Set();
		
		// for (var key in faceIndex)
		// {	
		//var face = geometry.mesh.faces[0];
		var edge_val = -1;
		
		// function helper(f, val) 
		// {
		// 	if (seen.has(f)) return;
			
		// 	for (let e of f.adjacentEdges())
		// 	{
		// 		T.addEntry(val, f.index, e.index);
		// 	}
		// 	seen.add(f);
			
		// 	for (let fadj of f.adjacentFaces())
		// 	{
		// 		helper(fadj, -1 * val);
		// 	}
		// }
		// for (let face of geometry.mesh.faces)
		// 	helper(face, edge_val);
		// console.log(seen.size);
		// console.log(F);
		//}

		var face = geometry.mesh.faces[0];
		
		return SparseMatrix.fromTriplet(T);
	}
}
