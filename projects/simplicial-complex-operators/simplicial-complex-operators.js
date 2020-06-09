"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                mesh.indexElements();
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                // using mesh.A0
                let V = mesh.vertices.length;
                let E = mesh.edges.length;
                
                let T = new Triplet(E, V);
                for (let i=0; i<V; i++) {
                        let v = mesh.vertices[i];
                        for (let e of v.adjacentEdges()) {
                                T.addEntry(1, e.index, v.index); 
                        }
                }
                return SparseMatrix.fromTriplet(T);
                
        }

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.SparseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                let F = mesh.faces.length;
                let E = mesh.edges.length;
                
                let T = new Triplet(E, F);
                for (let i=0; i<F; i++){
                        let f = mesh.faces[i];
                        for (let e of f.adjacentEdges()) {
                                T.addEntry(1, e.index, f.index);
                        }   
                }
                return SparseMatrix.fromTriplet(T);    
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexO
perators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                let V = this.mesh.vertices.length;
                let A = DenseMatrix.zeros(V, 1);
                
                for (let i=0; i<V; i++){
                        if (subset.vertices.has(i)){
                                A.set(1, i, 0);
                        }
                }
                return A;
                
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                let E = this.mesh.edges.length;
                let A = DenseMatrix.zeros(E, 1);
                
                for (let i=0; i<E; i++){
                        if (subset.edges.has(i)){
                                A.set(1, i, 0);
                        }
                }
                return A;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                let F = this.mesh.faces.length;
                let A = DenseMatrix.zeros(F, 1);
                
                for (let i=0; i<F; i++){
                        if (subset.faces.has(i)){
                                A.set(1, i, 0);
                        }
                }
                return A;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                /**
                 * What happens if you apply two unsigned adjacency matrices in sequence?
                 * How do you get all the simplices in the star?
                 */
                let subset_star = new MeshSubset();
                var FV = this.A0.transpose().timesSparse(this.A1);
                assert (FV.nRows() == this.mesh.vertices.length);
                assert (FV.nCols() == this.mesh.faces.length);
                
        
                for (let v of subset.vertices) //v is an index
                {       
                        let dmat = FV.subMatrix(v, v+1, 0, FV.nCols()).toDense();
                        
                        
                        var cmp1 = dmat.sum(); //should be double cmp1
                        //console.log("cmp1:", cmp1);

                        let cmp2 = 0;
                        let ve_row = this.A0.transpose().subMatrix(v, v+1, 0, this.A0.nRows()).toDense();
                        
                        let tmp_edges = new Set();
                        let tmp_faces = new Set();
                        //for (let e of subset.edges)
                        for (let e=0; e<ve_row.nCols(); e++)
                        {
                                if (ve_row.get(0, e) == 1) 
                                {
                                        
                                        tmp_edges.add(e); 
                                        let fe_row = this.A1.subMatrix(e, e+1, 0, this.A1.nCols()).toDense();
                                        for (let f=0; f<fe_row.nCols(); f++)
                                        //for (let f of subset.faces)
                                        {
                                                if (fe_row.get(0, f) == 1) 
                                                {
                                                        cmp2++;
                                                        tmp_faces.add(f); 
                                                }
                                        }
                                }
                                
                                //console.log("cmp2:", cmp2);
                        }

                        //if (cmp1 == cmp2)
                        //{
                                subset_star.addVertex(v);
                                subset_star.addEdges(tmp_edges);
                                subset_star.addFaces(tmp_faces);
                        //}
                }
                for (let e of subset.edges)
                {
                        subset_star.addEdge(e);
                        let fe_row = this.A1.subMatrix(e, e+1, 0, this.A1.nCols()).toDense();
                        for (let f=0; f<fe_row.nCols(); f++)
                                        //for (let f of subset.faces)
                                        {
                                                if (fe_row.get(0, f) == 1) 
                                                {
                                                
                                                        subset_star.addFace(f); 
                                                }
                                        }
                }

                for (let f of subset.faces)
                {
                        subset_star.addFace(f);
                }

                return subset_star; 
        }
        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                // TODO
                //let subset_star = this.star(subset);
                //var VF = this.A0.transpose().timesSparse(this.A1); 
                let cl = MeshSubset.deepCopy(subset);
                for (let f of cl.faces){
                //for (let f of subset_star.faces){ //for each face add all edges
                        let ef_row = this.A1.subMatrix(0, this.A1.nRows(), f, f+1).toDense();
                        for (let e=0; e<this.A1.nRows(); e++)
                        {       
                                if (ef_row.get(e,0) == 1) cl.addEdge(e);
                        }      
                }
                
                for (let e of cl.edges)
                {
                        let ev_row = this.A0.subMatrix(e,e+1, 0, this.A0.nCols()).toDense();
                        for (let v=0; v<this.A0.nCols(); v++)
                        {       
                                if (ev_row.get(0,v) == 1) cl.addVertex(v);
                        }   
                }

                // let ev_col = this.A0.subMatrix(0, this.A0.nRows(), f, f+1).toDense();
                //         for (let v=0; v<VF.nRows(); v++)
                //         {       
                //                 if (vf_col.get(v,0) == 2) subset.addVertex(v);
                //                 //vertex, face appear on 2 edges
                //         }
                return cl; // placeholder
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                // TODO
                let lk = this.closure(this.star(subset));
                lk.deleteSubset(this.star(this.closure(subset)));
                return lk;
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                //use closure method 
                // is subset a simplicial complex?
                return subset.equals(this.closure(subset));
                // let cl = this.closure(this.star(subset));
                // let euler = cl.vertices.size - cl.edges.size + cl.faces.size;
                // console.log("euler:", euler)
                // let euler_sub = subset.vertices.size - subset.edges.size + subset.faces.size;
                // console.log("sub:" ,euler_sub)
                // return euler_sub == euler;
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                /**
                 * Return degre of complex if complex else -1
                 */
                //use adjacency matries
                let K = 0;
                if (!this.isComplex(subset)) return -1;

                if (subset.edges.size > 0)
                {
                        K++;
                        for (let v of subset.vertices)
                        {
                                let v_attach = 0;
                                let ev_row = this.A0.subMatrix(0, this.A0.nRows(), v, v+1).toDense();
                                
                                for (let e=0; e< this.A0.nRows(); e++)
                                {
                                        if (ev_row.get(e,0) == 1 && subset.edges.has(e)) v_attach++;
                                }
                                if (v_attach==0) return -1
                        }
                        
                }

                if (subset.faces.size > 0)
                {
                        K++;
                        for (let e of subset.edges)
                        {
                                let e_attach = 0;
                                let ef_row = this.A1.subMatrix(e,e+1, 0, this.A1.nCols()).toDense();
                                
                                for (let f=0; f< this.A1.nCols(); f++)
                                {
                                        if (ef_row.get(0,f) == 1 && subset.faces.has(f)) e_attach++;
                                }
                                if (e_attach==0) return -1
                        }
                }
                return K;
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                // TODO
                /**
                 * Use isPure first bc boundary not defined for unpure simplicial complex
                 * Hint: think carefully about what the result of applyiing an unsigned adjacency
                 *    matrix can look like. What do you notice about simplices that should be in 
                 *    output set. Use Exercise 9
                 */
                //assert (this.isPureComplex(subset));
                let bd = new MeshSubset();
                
                if (!subset.faces.size && subset.edges.size) 
                {
                        for (let v of subset.vertices)
                        {
                                let v_check=0;
                                let ev = this.A0.subMatrix(0, this.A0.nRows(), v, v+1).toDense();
                                for (let e=0; e<this.A0.nRows(); e++)
                                {
                                        if (ev.get(e,0) == 1 && subset.edges.has(e)) v_check++;
                                }
                                console.log(v_check);
                                if (v_check==1) bd.addVertex(v);
                        }
                }
                for (let e of subset.edges)
                {

                        let ef = this.A1.subMatrix(e,e+1, 0, this.A1.nCols()).toDense();
                        let proper_face_chk = 0;
                        for (let f = 0; f< this.A1.nCols(); f++)
                        {
                                if (ef.get(0,f)== 1 && subset.faces.has(f))  proper_face_chk++;
                        }
                        if (proper_face_chk == 1 )
                        {
                                bd.addEdge(e);
                                let ev = this.A0.subMatrix(e,e+1, 0, this.A0.nCols()).toDense();
                                for (let v=0; v<this.A0.nCols(); v++)
                                {
                                        if (ev.get(0,v) == 1) bd.addVertex(v);
                                }
                        }
                        
                        // if (proper_face_chk == 0)
                        // {   
                        //         let proper_vertex_check = 0;
                        //         let ev = this.A0.subMatrix(e,e+1, 0, this.A0.nCols()).toDense();
                        //         for (let v=0; v<this.A0.nCols(); v++)
                        //         {
                        //                 if (ev.get(0,v) == 1 && subset.vertices.has(v))
                        //                 {
                        //                         proper_vertex_check++;
                        //                         bd.addVertex(v);
                        //                 }
                        //         }

                        //         if (proper_vertex_check <= 2) bd.addEdge(e);
                        // }
                                 
                }
                return bd; // placeholder
        }
}
