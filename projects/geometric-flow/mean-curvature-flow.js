"use strict";

class MeanCurvatureFlow {
	/**
	 * This class performs {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf mean curvature flow} on a surface mesh.
	 * @constructor module:Projects.MeanCurvatureFlow
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);
	}

	/**
	 * Builds the mean curvature flow operator.
	 * @private
	 * @method module:Projects.MeanCurvatureFlow#buildFlowOperator
	 * @param {module:LinearAlgebra.SparseMatrix} M The mass matrix of the input mesh.
	 * @param {number} h The timestep.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	buildFlowOperator(M, h) {
		// TODO
		var L = this.geometry.laplaceMatrix(this.vertexIndex);
		//var I = SparseMatrix.identity(M.nRows(), M.nRows());
		
		//return I.plus(M.invertDiagonal().timesSparse(L).timesReal(h));
		return M.plus(L.timesReal(h));
		//return SparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * Performs mean curvature flow on the input mesh with timestep h.
	 * @method module:Projects.MeanCurvatureFlow#integrate
	 * @param {number} h The timestep.
	 */
	integrate(h) {
		// TODO
		let vertices = this.geometry.mesh.vertices;
		var M = this.geometry.massMatrix(this.vertexIndex);
		var A = this.buildFlowOperator(M, h);

		// need to moves vertices into a #V by 3 matrix
		var pos = new Triplet(M.nRows(), 3);
		
		for (let v of vertices) 
		{
			let p = this.geometry.positions[v];
			pos.addEntry(p.x, v.index, 0);
			pos.addEntry(p.y, v.index, 1);
			pos.addEntry(p.z, v.index, 2);
		}
		var b = M.timesDense(SparseMatrix.fromTriplet(pos).toDense());
		
		var llt = A.chol();
		var res = llt.solvePositiveDefinite(b);
		for (let v of vertices)  
		{
			let p = this.geometry.positions[v];
			p.x = res.get(v.index,0);
			p.y = res.get(v.index,1);
			p.z = res.get(v.index,2);
		}
		console.log(b.sum() - res.sum());
		// center mesh positions around origin
		normalize(this.geometry.positions, vertices, false);
	}
}
