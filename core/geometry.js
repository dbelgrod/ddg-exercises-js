"use strict";

class Geometry {
	/**
	 * This class represents the geometry of a {@link module:Core.Mesh Memjsh}. This includes information such
	 * as the position of vertices as well as methods to compute edge lengths, corner
	 * angles, face area, normals, discrete curvatures etc.
	 * @constructor module:Core.Geometry
	 * @param {module:Core.Mesh} mesh The mesh this class describes the geometry of.
	 * @param {module:LinearAlgebra.Vector[]} positions An array containing the position of each vertex in a mesh.
	 * @param {boolean} normalizePositions flag to indicate whether positions should be normalized. Default value is true.
	 * @property {module:Core.Mesh} mesh The mesh this class describes the geometry of.
	 * @property {Object} positions A dictionary mapping each vertex to a normalized position.
	 */
	constructor(mesh, positions, normalizePositions = true) {
		this.mesh = mesh;
		this.positions = {};
		for (let i = 0; i < positions.length; i++) {
			let v = this.mesh.vertices[i];
			let p = positions[i];

			this.positions[v] = p;
		}

		if (normalizePositions) {
			normalize(this.positions, mesh.vertices);
		}
	}

	/**
	 * Computes the vector along a halfedge.
	 * @method module:Core.Geometry#vector
	 * @param {module:Core.Halfedge} h The halfedge along which the vector needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vector(h) {
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];

		return b.minus(a);
	}

	/**
	 * Computes the length of an edge.
	 * @method module:Core.Geometry#length
	 * @param {module:Core.Edge} e The edge whose length needs to be computed.
	 * @returns {number}
	 */
	length(e) {
		return this.vector(e.halfedge).norm();
	}

	/**
	 * Computes the midpoint of an edge.
	 * @method module:Core.Geometry#midpoint
	 * @param {module:Core.Edge} e The edge whose midpoint needs to be computed.
	 * @returns {number}
	 */
	midpoint(e) {
		let h = e.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.twin.vertex];

		return (a.plus(b)).over(2);
	}

	/**
	 * Computes the mean edge length of all the edges in a mesh.
	 * @method module:Core.Geometry#meanEdgeLength
	 * @returns {number}
	 */
	meanEdgeLength() {
		let sum = 0;
		let edges = this.mesh.edges;
		for (let e of edges) {
			sum += this.length(e);
		}

		return sum / edges.length;
	}

	/**
	 * Computes the area of a face.
	 * @method module:Core.Geometry#area
	 * @param {module:Core.Face} f The face whose area needs to be computed.
	 * @returns {number}
	 */
	area(f) {
		if (f.isBoundaryLoop()) return 0.0;

		let u = this.vector(f.halfedge);
		let v = this.vector(f.halfedge.prev).negated();

		return 0.5 * u.cross(v).norm();
	}

	/**
	 * Computes the total surface area of a mesh.
	 * @method module:Core.Geometry#totalArea
	 * @returns {number}
	 */
	totalArea() {
		let sum = 0.0;
		for (let f of this.mesh.faces) {
			sum += this.area(f);
		}

		return sum;
	}

	/**
	 * Computes the normal of a face.
	 * @method module:Core.Geometry#faceNormal
	 * @param {module:Core.Face} f The face whose normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	faceNormal(f) {
		if (f.isBoundaryLoop()) return undefined;

		let u = this.vector(f.halfedge);
		let v = this.vector(f.halfedge.prev).negated();

		return u.cross(v).unit();
	}

	/**
	 * Computes the centroid of a face.
	 * @method module:Core.Geometry#centroid
	 * @param {module:Core.Face} f The face whose centroid needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	centroid(f) {
		let h = f.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];
		let c = this.positions[h.prev.vertex];

		if (f.isBoundaryLoop()) return a.plus(b).over(2);

		return a.plus(b).plus(c).over(3);
	}

	/**
	 * Computes the circumcenter of a face.
	 * @method module:Core.Geometry#circumcenter
	 * @param {module:Core.Face} f The face whose circumcenter needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	circumcenter(f) {
		let h = f.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];
		let c = this.positions[h.prev.vertex];

		if (f.isBoundaryLoop()) return a.plus(b).over(2);

		let ac = c.minus(a);
		let ab = b.minus(a);
		let w = ab.cross(ac);

		let u = (w.cross(ab)).times(ac.norm2());
		let v = (ac.cross(w)).times(ab.norm2());
		let x = (u.plus(v)).over(2 * w.norm2());

		return x.plus(a);
	}

	/**
	 * Computes an orthonormal bases for a face.
	 * @method module:Core.Geometry#orthonormalBases
	 * @param {module:Core.Face} f The face on which the orthonormal bases needs to be computed.
	 * @returns {module:LinearAlgebra.Vector[]} An array containing two orthonormal vectors tangent to the face.
	 */
	orthonormalBases(f) {
		let e1 = this.vector(f.halfedge).unit();

		let normal = this.faceNormal(f);
		let e2 = normal.cross(e1);

		return [e1, e2];
	}

	/**
	 * Computes the angle (in radians) at a corner.
	 * @method module:Core.Geometry#angle
	 * @param {module:Core.Corner} c The corner at which the angle needs to be computed.
	 * @returns {number} The angle clamped between 0 and π.
	 */
	angle(c) {
		//halfedge is opposite to corner
		var a = this.vector(c.halfedge.next).negated();
		var b = this.vector(c.halfedge.prev);

		var x = a.dot(b);
		x /= a.norm();
		x /= b.norm();

		return Math.acos(x); //returns 0 to pi
		
		//return 0.0; // placeholder
	}

	/**
	 * Computes the cotangent of the angle opposite to a halfedge.
	 * @method module:Core.Geometry#cotan
	 * @param {module:Core.Halfedge} h The halfedge opposite to the angle whose cotangent needs to be computed.
	 * @returns {number}
	 */
	cotan(h) {
		// TODO
		var a = this.vector(h.next).negated();
		var b = this.vector(h.prev);

		return a.dot(b) / a.cross(b).norm();
	}

	/**
	 * Computes the signed angle (in radians) between two adjacent faces.
	 * @method module:Core.Geometry#dihedralAngle
	 * @param {module:Core.Halfedge} h The halfedge (shared by the two adjacent faces) on which
	 * the dihedral angle is computed.
	 * @returns {number} The dihedral angle.
	 */
	dihedralAngle(h) {
		// TODO
		var e_n = this.vector(h).unit();
		var n1 = this.faceNormal(h.face);
		var n2 = this.faceNormal(h.twin.face);

		var y = e_n.dot(n1.cross(n2));
		var x = n1.dot(n2);

		return Math.atan2(y,x);
		//return 0.0; // placeholder
	}

	/**
	 * Computes the barycentric dual area of a vertex.
	 * @method module:Core.Geometry#barycentricDualArea
	 * @param {module:Core.Vertex} v The vertex whose barycentric dual area needs to be computed.
	 * @returns {number}
	 */
	barycentricDualArea(v) {
		// TODO;
		var dualArea = 0;
		for (let f of v.adjacentFaces()) 
		{
			dualArea += this.area(f);
		}
		return dualArea / 3; // placeholder
	}

	/**
	 * Computes the circumcentric dual area of a vertex.
	 * @see {@link http://www.cs.cmu.edu/~kmcrane/Projects/Other/TriangleAreasCheatSheet.pdf}
	 * @method module:Core.Geometry#circumcentricDualArea
	 * @param {module:Core.Vertex} v The vertex whose circumcentric dual area needs to be computed.
	 * @returns {number}
	 */
	circumcentricDualArea(v) {
		// TODO
		var area = 0.0;
		for (let h of v.adjacentHalfedges())
		{
			if (h.vertex == v)
			{
				var e1 = this.vector(h);
				var e2 = this.vector(h.prev).negated();

				var alpha_cot = this.cotan(h);
				var beta_cot = this.cotan(h.prev);

				area += e1.norm2()*alpha_cot + e2.norm2()*beta_cot;
			}
		}
		return area / 8;
		//return 0.0; // placeholder
	}

	/**
	 * Computes the normal at a vertex using the "equally weighted" method.
	 * @method module:Core.Geometry#vertexNormalEquallyWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalEquallyWeighted(v) {
		let n = new Vector();
		for (let f of v.adjacentFaces()) {
			let normal = this.faceNormal(f);

			n.incrementBy(normal);
		}

		n.normalize();

		return n;
	}

	/**
	 * Computes the normal at a vertex using the "face area weights" method.
	 * @method module:Core.Geometry#vertexNormalAreaWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalAreaWeighted(v) {
		// TODO
		let n = new Vector();
		for (let f of v.adjacentFaces())
		{
			let normal = this.faceNormal(f);
			let area = this.area(f);
			n.incrementBy(normal.times(area));	
		}
		n.normalize();

		return n;
		//return new Vector(); // placeholder
	}

	/**
	 * Computes the normal at a vertex using the "tip angle weights" method.
	 * @method module:Core.Geometry#vertexNormalAngleWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalAngleWeighted(v) {
		// TODO
		let n = new Vector();
		for (let f of v.adjacentFaces())
		{
			let normal = this.faceNormal(f);
			for (let c of f.adjacentCorners())
			{
				if (c.vertex == v) n.incrementBy(normal.times(this.angle(c)));
			}
			
		}
		n.normalize();

		return n;
		//return new Vector(); // placeholder
	}

	/**
	 * Computes the normal at a vertex using the "gauss curvature" method.
	 * @method module:Core.Geometry#vertexNormalGaussCurvature
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalGaussCurvature(v) {
		// TODO
		var n = new Vector();
		for (let e of v.adjacentEdges())
		{
			var gc = this.vector(e.halfedge);
			if (e.halfedge.vertex != v) gc = gc.negated();
			gc.scaleBy(this.dihedralAngle(e.halfedge) / this.vector(e.halfedge).norm());
			n.incrementBy(gc);
		}
		//return new Vector(); // placeholder
		n.scaleBy(0.5);
		n.normalize();
		return n;
	}

	/**
	 * Computes the normal at a vertex using the "mean curvature" method (same as the "area gradient" method).
	 * @method module:Core.Geometry#vertexNormalMeanCurvature
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalMeanCurvature(v) {
		// TODO
		var n = new Vector();
		for (let e of v.adjacentEdges())
		{
			var nc = this.vector(e.halfedge);
			if (e.halfedge.vertex != v) nc = nc.negated();
			nc.scaleBy(this.cotan(e.halfedge) + this.cotan(e.halfedge.twin))
			n.incrementBy(nc);
		}
		//return new Vector(); // placeholder
		n.scaleBy(0.5);
		n.normalize();
		return n;
		//return new Vector(); // placeholder
	}

	/**
	 * Computes the normal at a vertex using the "inscribed sphere" method.
	 * @method module:Core.Geometry#vertexNormalSphereInscribed
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalSphereInscribed(v) {
		// TODO
		let n = new Vector();
		for (let h of v.adjacentHalfedges())
		{
			if (h.vertex == v)
			{
				var a = this.vector(h);
				var b = this.vector(h.prev).negated()

				var sphere = a.cross(b);
				//console.log(sphere);

				sphere.divideBy(a.norm2());
				sphere.divideBy(b.norm2());
			
				n.incrementBy(sphere);
			}

		}

		n.normalize();
		return n;
		//return new Vector(); // placeholder
	}

	/**
	 * Computes the angle defect at a vertex (= 2π minus the sum of incident angles
	 * at an interior vertex or π minus the sum of incident angles at a boundary vertex).
	 * @method module:Core.Geometry#angleDefect
	 * @param {module:Core.Vertex} v The vertex whose angle defect needs to be computed.
	 * @returns {number}
	 */
	angleDefect(v) {
		// TODO
		if (v.onBoundary()) var total = Math.PI;
		else var total = 2*Math.PI;

		let inc = 0;
		for (let c of v.adjacentCorners())
		{
			if (c.vertex == v)
			{
				inc += this.angle(c);
			}
		}
		return total - inc;
		//return 0.0; // placeholder
	}

	/**
	 * Computes the (integrated) scalar gauss curvature at a vertex.
	 * @method module:Core.Geometry#scalarGaussCurvature
	 * @param {module:Core.Vertex} v The vertex whose gauss curvature needs to be computed.
	 * @returns {number}
	 */
	scalarGaussCurvature(v) {
		return this.angleDefect(v);
	}

	/**
	 * Computes the (integrated) scalar mean curvature at a vertex.
	 * @method module:Core.Geometry#scalarMeanCurvature
	 * @param {module:Core.Vertex} v The vertex whose mean curvature needs to be computed.
	 * @returns {number}
	 */
	scalarMeanCurvature(v) {
		// TODO
		var H = 0.0;
		for (let e of v.adjacentEdges())
		{
			H += this.dihedralAngle(e.halfedge) * this.vector(e.halfedge).norm();
		}
		return 0.5* H;
		//return 0.0; // placeholder
	}

	/**
	 * Computes the total angle defect (= 2π times the euler characteristic of the mesh).
	 * @method module:Core.Geometry#totalAngleDefect
	 * @returns {number}
	 */
	totalAngleDefect() {
		// TODO
		return 2*Math.PI*this.mesh.eulerCharacteristic();
		//return 0.0; // placeholder
	}

	/**
	 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
	 * @method module:Core.Geometry#principalCurvatures
	 * @param {module:Core.Vertex} v The vertex on which the principal curvatures need to be computed.
	 * @returns {number[]} An array containing the minimum and maximum principal curvature values at a vertex.
	 */
	principalCurvatures(v) {
		// TODO
		var H =  this.scalarMeanCurvature(v);
		var root = Math.sqrt(Math.pow(this.scalarMeanCurvature, 2) - this.scalarGaussCurvature);

		return [H-root, H+root];
		//return [0.0, 0.0]; // placeholder
	}

	/**
	 * Builds a sparse laplace matrix. The laplace operator is negative semidefinite;
	 * instead we build a positive definite matrix by multiplying the entries of the
	 * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
	 * @method module:Core.Geometry#laplaceMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	laplaceMatrix(vertexIndex) {
		// TODO
		var size = Object.keys(vertexIndex).length;
		var L = new Triplet(size,size);
		
		for (let key in vertexIndex){
			var pos = vertexIndex[key];
			var i = this.mesh.vertices[pos];
			var lval = 0;
			for (let h of i.adjacentHalfedges())
			{
				if (h.vertex.index == i.index) 
				{
					var j = h.next.vertex;
					var u = -0.5*(this.cotan(h) + this.cotan(h.twin));
				}
				else
				{
					var j = h.vertex;
					var u = 0.5*(this.cotan(h) + this.cotan(h.twin));
				}
				L.addEntry(u, i.index, j.index);
				lval += u;
			}
			L.addEntry(-lval+1e-8, i.index, i.index);
		}
		return SparseMatrix.fromTriplet(L);
		//return SparseMatrix.identity(1, 1); // placeholder
	}

	/**
	 * Builds a sparse diagonal mass matrix containing the barycentric dual area of each vertex
	 * of a mesh.
	 * @method module:Core.Geometry#massMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	massMatrix(vertexIndex) {
		// TODO
		var size = Object.keys(vertexIndex).length;
		var L = new Triplet(size,size);
		
		for (let key in vertexIndex)
		{
			var pos = vertexIndex[key];
			var v = this.mesh.vertices[pos];
			L.addEntry(this.barycentricDualArea(v), v.index, v.index);
		}
		//return SparseMatrix.identity(1, 1); // placeholder
		return SparseMatrix.fromTriplet(L);
	}

	/**
	 * Builds a sparse complex laplace matrix. The laplace operator is negative semidefinite;
	 * instead we build a positive definite matrix by multiplying the entries of the
	 * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
	 * @method module:Core.Geometry#complexLaplaceMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.ComplexSparseMatrix}
	 */
	complexLaplaceMatrix(vertexIndex) {
		// TODO
		var size = Object.keys(vertexIndex).length;
		var L = new ComplexTriplet(size,size);

		for (let key in vertexIndex){
			var pos = vertexIndex[key];
			var i = this.mesh.vertices[pos];
			var lval = 0;
			for (let h of i.adjacentHalfedges())
			{
				if (h.vertex.index == i.index) 
				{
					var j = h.next.vertex;
					var u = -0.25*(this.cotan(h) + this.cotan(h.twin)); //ok
					if (h.twin.onBoundary) u = -0.25*(this.cotan(h))
					if (h.onBoundary) u = -0.25*(this.cotan(h.twin))
				}
				else
				{
					var j = h.vertex;
					var u = 0.25*(this.cotan(h) + this.cotan(h.twin));
					if (h.twin.onBoundary) u = 0.25*(this.cotan(h))
					if (h.onBoundary) u = 0.25*(this.cotan(h.twin))
				}
				L.addEntry(new Complex(u,0), i.index, j.index);
				lval += u;
			}
			L.addEntry(new Complex(-lval+1e-8,0), i.index, i.index);
			
		}
		return ComplexSparseMatrix.fromTriplet(L);
		//return ComplexSparseMatrix.identity(1, 1); // placeholder
	}
}

/**
 * Centers a mesh about the origin and rescales it to unit radius.
 * @global
 * @function module:Core.normalize
 * @param {module:LinearAlgebra.Vector[]} positions The position of each vertex in the vertices array.
 * @param {module:Core.Vertex[]} vertices The vertices of a mesh.
 * @param {boolean} rescale A flag indicating whether mesh positions should be scaled to a unit radius.
 */
function normalize(positions, vertices, rescale = true) {
	// compute center of mass
	let N = vertices.length;
	let cm = new Vector();
	for (let v of vertices) {
		let p = positions[v];

		cm.incrementBy(p);
	}
	cm.divideBy(N);

	// translate to origin and determine radius
	let radius = -1;
	for (let v of vertices) {
		let p = positions[v];

		p.decrementBy(cm);
		radius = Math.max(radius, p.norm());
	}

	// rescale to unit radius
	if (rescale) {
		for (let v of vertices) {
			let p = positions[v];

			p.divideBy(radius);
		}
	}
}