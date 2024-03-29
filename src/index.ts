import { view_index_triangle } from "./vis/view_indexed_triangle";
import { load_file } from "./util/load-file";
import { parse_obj } from "./util/parse-obj";
import { reset_elapsed_time } from "./util/timer";
import { init_three } from "./vis/init_three";

import { indexed_triangle_to_halfedge_mesh } from "./geo/halfedge_mesh_from_indexed_triangle"
import { halfedge_mesh_to_index_triangle } from "./geo/halfedge_mesh_to_index_triangle"

import "../geometry-processing-js/linear-algebra/linear-algebra-asm";
import "../geometry-processing-js/linear-algebra/vector";
import "../geometry-processing-js/linear-algebra/emscripten-memory-manager";
import "../geometry-processing-js/linear-algebra/dense-matrix";
import "../geometry-processing-js/linear-algebra/sparse-matrix";

import { SVD } from 'svd-js'
import { Halfedge } from "geo/atom_halfedge";
import { dot_v, magnitude_v, sub_v } from "./math/linalg_standard";
import { vertex_get_neighbors } from "./geo/vertex_get_neighbors";

import * as math from "mathjs";
import { V3 } from "./geo/atom_vertex";

// <script src="../../linear-algebra/linear-algebra-asm.js"></script>
// <script src="../../linear-algebra/vector.js"></script>
// <script src="../../linear-algebra/emscripten-memory-manager.js"></script>
// <script src="../../linear-algebra/dense-matrix.js"></script>
// <script src="../../linear-algebra/sparse-matrix.js"></script>

init_three();

// load_file('/basic-arap/meshes/armadillo.obj').then(async (obj_raw_data: string) => {
load_file('/basic-arap/meshes/tetrahedron.obj').then(async (obj_raw_data: string) => {
// load_file('/basic-arap/meshes/cube.obj').then(async (obj_raw_data: string) => {
    reset_elapsed_time();
    const obj_data = parse_obj(obj_raw_data);
    const mesh = indexed_triangle_to_halfedge_mesh(obj_data.vertices, obj_data.indices);

    let N = mesh.verts.length;
    let v_to_i: {[key: symbol]: number} = {};
    mesh.verts.map((v, i) => v_to_i[v.id] = i);

    const fixed_points = [
        {
            index: 0,
            p_fixed: mesh.verts[0].add_xyz(0, 1, 0)
        },
        // {
        //     index: 1,
        //     p_fixed: mesh.verts[1]
        // },
    ];
    const fixed_indices = fixed_points.map(fixed_data => fixed_data.index);

    let num_iters = 3;

    for(let i = 0; i < num_iters; i++){
        console.log(`starting iteration ${i}...`)

        // calculate weights
        let weights: number[][] = new Array(N).fill(0).map(_ => new Array(N).fill(0));
        mesh.halfedges.forEach(he => {
            let triangle_a = [he.vert, he.next.vert, he.next.next.vert];        
            let t_a_e1 = sub_v(triangle_a[2], triangle_a[1]);
            let t_a_e2 = sub_v(triangle_a[0], triangle_a[1]);

            let triangle_b = [he.twin.vert, he.twin.next.vert, he.twin.next.next.vert];
            let t_b_e1 = sub_v(triangle_b[2], triangle_b[1]);
            let t_b_e2 = sub_v(triangle_b[0], triangle_b[1]);

            let alpha = Math.acos(dot_v(t_a_e1, t_a_e2) / (magnitude_v(t_a_e1) * magnitude_v(t_a_e2)));
            let beta = Math.acos(dot_v(t_b_e1, t_b_e2) / (magnitude_v(t_b_e1) * magnitude_v(t_b_e2)));

            let weight = 0.5 * (1 / Math.tan(alpha) + 1 / Math.tan(beta));
            // weights[he.vert.id][he.next.vert.id] = weight;
            let index_i = v_to_i[he.vert.id];
            let index_j = v_to_i[he.next.vert.id]
            weights[index_i][index_j] = weight;
        });

        // calculate R_i
        const S_i = mesh.verts.map((i, c) => {
            let neighbors = vertex_get_neighbors(i).verts;
            let e_ij = neighbors.map(j => sub_v(j, i));
            let e_ij_p = neighbors.map(j => sub_v(j, i));
            let w_ij = neighbors.map(j => {
                let index_i = v_to_i[i.id];
                let index_j = v_to_i[j.id];
                return weights[index_i][index_j]
            });

            const P_i = math.transpose(math.matrix(e_ij.map(e => e.to_THREE())));
            const P_i_pt = math.matrix(e_ij_p.map(e => e.to_THREE()));
            const D = math.diag(w_ij);

            const S_i = math.multiply(P_i, math.multiply(D, P_i_pt)).toArray() as number[][];
            return S_i;
        }); 
        const R_i = S_i.map((s, i) => {
            const {u, v, q} = SVD(s, true, true);
            const U_it = math.transpose(math.matrix(u));
            const V_i = math.matrix(v);

            const R_i = math.multiply(V_i, U_it);
            // return math.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

            const det = math.det(R_i);
            if(det > 0){
                return R_i;
            } else {
                console.log("flipping sign of U_i....")
                const smallest_singular_index = q.indexOf(Math.min(...q));
                const i0 = [0, smallest_singular_index];
                const i1 = [1, smallest_singular_index];
                const i2 = [2, smallest_singular_index];
                U_it.set(i0, U_it.get(i0) * -1);
                U_it.set(i1, U_it.get(i1) * -1);
                U_it.set(i2, U_it.get(i2) * -1);
                const R_i = math.multiply(V_i, U_it);
                return R_i;
            }
        });

        // set up L and solve Lp'=b for p'
        const L_data = new Triplet(N, N);
        const b = DenseMatrix.zeros(N, 3);
        mesh.verts.map((i) => {
            let index_i = v_to_i[i.id];

            if(fixed_indices.includes(index_i)){
                const fixed_point = fixed_points.find(data => data.index === index_i) as {
                    index: number;
                    p_fixed: V3;
                };
                L_data.addEntry(1, index_i, index_i);
                b.set(fixed_point.p_fixed.x, index_i, 0);
                b.set(fixed_point.p_fixed.y, index_i, 1);
                b.set(fixed_point.p_fixed.z, index_i, 2);
                return;
            }

            let coef_i_on_L = 0;
            let b_sum = math.matrix([0, 0, 0]);

            let neighbors = vertex_get_neighbors(i).verts; // console.log(neighbors.map(n => n.id));
            neighbors.map(j => {
                let index_j = v_to_i[j.id];
                let w_ij = weights[index_i][index_j];

                if(fixed_indices.includes(index_j)){
                    const fixed_point = fixed_points.find(data => data.index === index_j) as {
                        index: number;
                        p_fixed: V3;
                    };
                    let c_k = math.matrix(fixed_point.p_fixed.to_THREE());
                    let result = math.multiply(w_ij, c_k);
                    b_sum = math.add(b_sum, result);
                } 

                // LHS
                coef_i_on_L += w_ij;
                if(!fixed_indices.includes(index_j)){
                    L_data.addEntry(-w_ij, index_i, index_j);
                }

                // RHS
                let sum_R = math.add(R_i[index_i], R_i[index_j]);
                let diff_p = sub_v(mesh.verts[index_i], mesh.verts[index_j]).to_THREE();
                let result = math.multiply(w_ij * 0.5, math.multiply(sum_R, diff_p));
                b_sum = math.add(b_sum, result);
            });

            L_data.addEntry(coef_i_on_L, index_i, index_i);

            b.set(b_sum.get([0]), index_i, 0);
            b.set(b_sum.get([1]), index_i, 1);
            b.set(b_sum.get([2]), index_i, 2);
        });

        const L = SparseMatrix.fromTriplet(L_data);

        // const llt = L.chol();
        // let p_prime = llt.solvePositiveDefinite(b);
        const qr = L.qr();
        let new_p_prime = qr.solve(b);
        
        // const L_dense = L.toDense();
        // const L_extracted = new Array(N).fill(0).map((_, i) => new Array(N).fill(0).map((_, j) => {return L_dense.get(i, j);}));  
        // const b_extracted = new Array(N).fill(0).map((_, i) => new Array(3).fill(0).map((_, j) => {return b.get(i, j);}));  
        // console.log('$$$$$$$$$$$$$$$$$$$$$$$')
        // console.log(L, L.nRows(), L.nCols());
        // console.log(b, b.nRows(), b.nCols());
        // console.log('! Le', L_extracted)
        // console.log('! be', b_extracted)
        // console.log(p_prime, b.nRows(), b.nCols());
        // const p_prime_extracted = new Array(N).fill(0).map((_, i) => new Array(3).fill(0).map((_, j) => {return p_prime.get(i, j);}));  
        // let p_i = mesh.verts.map(v => v.to_THREE());
        // console.log('goal', p_i)
        // console.log(math.multiply(L_extracted, p_i))
        // console.log('! pe', p_prime_extracted)
        // console.log(math.multiply(L_extracted, p_prime_extracted))
        // console.log('$$$$$$$$$$$$$$$$$$$$$$$')

        // update verex positions
        for(let i = 0; i < N; i++){
            mesh.verts[i].set_xyz(
                new_p_prime.get(i, 0),
                new_p_prime.get(i, 1),
                new_p_prime.get(i, 2),
            );
        }
    }

    const re_convert = halfedge_mesh_to_index_triangle(mesh);
    view_index_triangle(re_convert, 0x00aaff);
});