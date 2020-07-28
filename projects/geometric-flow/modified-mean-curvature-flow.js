"use strict";

class ModifiedMeanCurvatureFlow extends MeanCurvatureFlow {
	/**
	 * This class performs a {@link http://cs.jhu.edu/~misha/MyPapers/SGP12.pdf modified version} of {@link https://www.cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf mean curvature flow} on a surface mesh.
	 * @constructor module:Projects.ModifiedMeanCurvatureFlow
	 * @augments module:Projects.MeanCurvatureFlow
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:LinearAlgebra.SparseMatrix} A The laplace matrix of the input mesh.
	 */
	constructor(geometry) {
		super(geometry); //derives from MCF
		//let vertexIndex = {...this.vertexIndex}; //ends up being useless but copies dicts
		this.A = this.geometry.laplaceMatrix(this.vertexIndex); //defined before in buildFlowOperator
		// TODO: build the laplace matrix
		//this.A = SparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * @inheritdoc
	 */
	buildFlowOperator(M, h) {
		// TODO
		return M.plus(this.A.timesReal(h));
		//return SparseMatrix.identity(1, 1); // placeholder
	}
}
