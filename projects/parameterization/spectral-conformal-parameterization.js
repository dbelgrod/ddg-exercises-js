"use strict";

class SpectralConformalParameterization {
	/**
	 * This class implements the {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf spectral conformal parameterization} algorithm to flatten
	 * surface meshes with boundaries conformally.
	 * @constructor module:Projects.SpectralConformalParameterization
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);
	}

	/**
	 * Builds the complex conformal energy matrix EC = ED - A.
	 * @private
	 * @method module:Projects.SpectralConformalParameterization#buildConformalEnergy
	 * @returns {module:LinearAlgebra.ComplexSparseMatrix}
	 */
	buildConformalEnergy() {
		// TODO
		var size = Object.keys(this.vertexIndex).length;
		let ED = this.geometry.complexLaplaceMatrix(this.vertexIndex);

		let T = new ComplexTriplet(size, size);
		//console.log('size:', this.geometry.mesh.halfedges.length)
		for (let h of this.geometry.mesh.halfedges){
			if (h.onBoundary){
				let i = h.vertex.index;
				let j = h.next.vertex.index
				T.addEntry(new Complex(0, +1/4), i, j);
				//console.log(-1/4, i, j)
				T.addEntry(new Complex(0, -1/4), j, i);
			}
			
		}
		let A = ComplexSparseMatrix.fromTriplet(T);
		

		return ED.minus(A);
		
		//return ComplexSparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * Flattens the input surface mesh with 1 or more boundaries conformally.
	 * @method module:Projects.SpectralConformalParameterization#flatten
	 * @returns {Object} A dictionary mapping each vertex to a vector of planar coordinates.
	 */
	flatten() {
		// TODO
		let vertices = this.geometry.mesh.vertices;
		let flattening = this.geometry.positions; // placeholder
		
		//var EC = this.buildConformalEnergy();
		
		// normalize flattening
		normalize(flattening, vertices);

		return flattening;
	}
}
