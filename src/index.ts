import { view_index_triangle, view_index_triangle_updateable } from "./vis/view_indexed_triangle";
import { load_file } from "./util/load-file";
import { parse_obj } from "./util/parse-obj";
import { get_elapsed_time, get_elapsed_time2, get_elapsed_time3, reset_elapsed_time, reset_elapsed_time2, reset_elapsed_time3 } from "./util/timer";
import { camera, controls, init_three, renderer, scene } from "./vis/init_three";

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
import * as M from 'ml-matrix';

import * as THREE from "three";
import { DragControls } from "three/examples/jsm/controls/DragControls"

// <script src="../../linear-algebra/linear-algebra-asm.js"></script>
// <script src="../../linear-algebra/vector.js"></script>
// <script src="../../linear-algebra/emscripten-memory-manager.js"></script>
// <script src="../../linear-algebra/dense-matrix.js"></script>
// <script src="../../linear-algebra/sparse-matrix.js"></script>

init_three();

// load_file('/basic-arap/meshes/armadillo.obj').then(async (obj_raw_data: string) => {
// load_file('/basic-arap/meshes/tetrahedron.obj').then(async (obj_raw_data: string) => {
// load_file('/basic-arap/meshes/cube.obj').then(async (obj_raw_data: string) => {
// load_file('/basic-arap/meshes/sphere.obj').then(async (obj_raw_data: string) => {
load_file('/basic-arap/meshes/teapot.obj').then(async (obj_raw_data: string) => {
    reset_elapsed_time();
    const obj_data = parse_obj(obj_raw_data);
    const mesh = indexed_triangle_to_halfedge_mesh(obj_data.vertices, obj_data.indices);
    console.log(`Loaded mesh with ${mesh.verts.length} vertices`)

    let re_convert = halfedge_mesh_to_index_triangle(mesh);
    let viewing_object = view_index_triangle_updateable(re_convert, 0x00aaff);
    const debug_log = false;
    
    let N = mesh.verts.length;
    let v_to_i: {[key: symbol]: number} = {};
    mesh.verts.map((v, i) => v_to_i[v.id] = i);

    let weights: number[][];
    let use_QR = true;
    let L: any, qr: any, llt: any;
    let p: V3[], p_prime: V3[];

    let fixed_points: {
        index: number;
        p_fixed: V3;
    }[] = [
        {
            index: 0,
            p_fixed: mesh.verts[0].clone().add_xyz(0, 0, 0)
        },
        // {
        //     index: 20,
        //     p_fixed: mesh.verts[20].clone()
        // },
    ];
    let fixed_indices = fixed_points.map(fixed_data => fixed_data.index); 

    const calculate_L = () => {
        let re_convert = halfedge_mesh_to_index_triangle(mesh);
        let mesh_clone = indexed_triangle_to_halfedge_mesh(re_convert.vertices.map(v => v.clone_same_id()), re_convert.indices);

        p = mesh_clone.verts;
        p_prime = p.map(v => v.clone());
        for(let i = 0; i < N; i++) p[i].cache1 = p_prime[i];
    
        // calculate weights
        weights = new Array(N).fill(0).map(_ => new Array(N).fill(0));
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
    
        // set up L
        const L_data = new Triplet(N, N);
        p.map((p_i) => { // should be p_prime except we dont' have connectivity info
            let index_i = v_to_i[p_i.id];
    
            if(fixed_indices.includes(index_i)){
                L_data.addEntry(1, index_i, index_i);
                return;
            }
    
            let coef_i_on_L = 0;
    
            let neighbors = vertex_get_neighbors(p_i).verts;
            neighbors.map(p_j => {
                let index_j = v_to_i[p_j.id];
                let w_ij = weights[index_i][index_j];
    
                // LHS
                coef_i_on_L += w_ij;
                if(!fixed_indices.includes(index_j)){
                    L_data.addEntry(-w_ij, index_i, index_j);
                }
            });
    
            L_data.addEntry(coef_i_on_L, index_i, index_i);
        });
        L = SparseMatrix.fromTriplet(L_data);

        if(use_QR) qr = L.qr();
        else llt = L.chol();
    }

    const calculate_b = () => {
        reset_elapsed_time();
        reset_elapsed_time2();

        p_prime = p.map(v => v.clone());
        for(let i = 0; i < N; i++) p[i].cache1 = p_prime[i];

        let time_R_i = 0, time_b = 0, time_p_prime = 0;
        // optimization loop
        let num_iters = 2;
        for(let i = 0; i < num_iters; i++){
            if(debug_log) console.log(`starting iteration ${i}...`)
    
            // calculate R_i
            const S_i = p.map((i, c) => {
                let p_i = p[c];
                let p_i_prime = p_prime[c];
    
                let neighbors = vertex_get_neighbors(p_i).verts;
                let neighbors_prime = neighbors.map(p_i => p_i.cache1 as V3);
                let e_ij = neighbors.map(p_j => sub_v(p_j, p_i));
                let e_ij_prime = neighbors_prime.map(p_j_prime => sub_v(p_j_prime, p_i_prime));
    
                let w_ij = neighbors.map(j => {
                    let index_i = v_to_i[i.id];
                    let index_j = v_to_i[j.id];
                    return weights[index_i][index_j]
                });
    
                const P_i = math.transpose(math.matrix(e_ij.map(e => e.to_THREE())));
                const P_i_pt = math.matrix(e_ij_prime.map(e => e.to_THREE()));
                const D = math.diag(w_ij);

                const S_i = math.multiply(P_i, math.multiply(D, P_i_pt)).toArray() as number[][];
                return S_i;
            }); 
            const R_i = S_i.map((s, i) => {
                const svd_result = new M.SingularValueDecomposition(s);
                const u = svd_result.leftSingularVectors.to2DArray();
                const v = svd_result.rightSingularVectors.to2DArray();
                const q = svd_result.diagonal;
    
                const U_it = math.transpose(math.matrix(u));
                const V_i = math.matrix(v);
                const R_i = math.multiply(V_i, U_it);
    
                const det = math.det(R_i);
                if(det > 0){
                    return R_i;
                } else {
                    const qs = q.filter(n => !isNaN(n));
                    const smallest_singular_index = qs.indexOf(Math.min(...qs));
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
            time_R_i += parseFloat(get_elapsed_time(true));
    
            // set up b
            const b = DenseMatrix.zeros(N, 3);
            p.map((i) => {
                let index_i = v_to_i[i.id];
        
                if(fixed_indices.includes(index_i)){
                    const fixed_point = fixed_points.find(data => data.index === index_i) as {
                        index: number;
                        p_fixed: V3;
                    };
                    b.set(fixed_point.p_fixed.x, index_i, 0);
                    b.set(fixed_point.p_fixed.y, index_i, 1);
                    b.set(fixed_point.p_fixed.z, index_i, 2);
                    return;
                }
    
                let b_sum = math.matrix([0, 0, 0]);
        
                let neighbors = vertex_get_neighbors(i).verts; 
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
        
                    // RHS
                    let sum_R = math.add(R_i[index_i], R_i[index_j]);
                    let diff_p = sub_v(p[index_i], p[index_j]).to_THREE();
                    let result = math.multiply(w_ij * 0.5, math.multiply(sum_R, diff_p));
                    b_sum = math.add(b_sum, result);
                });
    
                b.set(b_sum.get([0]), index_i, 0);
                b.set(b_sum.get([1]), index_i, 1);
                b.set(b_sum.get([2]), index_i, 2);
                
            });
            time_b += parseFloat(get_elapsed_time(true));
    
            let new_p_prime: any;
            if(use_QR) new_p_prime = qr.solve(b);
            else new_p_prime = llt.solvePositiveDefinite(b);
            time_p_prime += parseFloat(get_elapsed_time(true));
            
            const L_dense = L.toDense();
            const L_extracted = new Array(N).fill(0).map((_, i) => new Array(N).fill(0).map((_, j) => {return L_dense.get(i, j);}));  
            const b_extracted = new Array(N).fill(0).map((_, i) => new Array(3).fill(0).map((_, j) => {return b.get(i, j);}));  
            console.log('$$$$$$$$$$$$$$$$$$$$$$$')
            console.log(L, L.nRows(), L.nCols());
            console.log(b, b.nRows(), b.nCols());
            console.log('! Le', L_extracted)
            console.log('! be', b_extracted)
            console.log(new_p_prime, b.nRows(), b.nCols());
            const p_prime_extracted = new Array(N).fill(0).map((_, i) => new Array(3).fill(0).map((_, j) => {return new_p_prime.get(i, j);}));  
            let p_i = mesh.verts.map(v => v.to_THREE());
            console.log('goal', p_i)
            console.log(math.multiply(L_extracted, p_i))
            console.log('! pe', p_prime_extracted)
            console.log(math.multiply(L_extracted, p_prime_extracted))
            console.log('$$$$$$$$$$$$$$$$$$$$$$$')
    
            // update vertex positions
            for(let i = 0; i < N; i++){
                p_prime[i].set_xyz(
                    new_p_prime.get(i, 0),
                    new_p_prime.get(i, 1),
                    new_p_prime.get(i, 2),
                );
            }
        }
    
        if(debug_log){
            console.log(`computed new positions based on fixed indices in ${get_elapsed_time2(true)}ms`)
            console.log(`    computed R_i's in ${time_R_i}ms`)
            console.log(`    computed b's in ${time_b}ms`)
            console.log(`    computed p_prime's in ${time_p_prime}ms`)
        }
    
        // set final geometry...
        // and update handle positions
        for(let i = 0; i < N; i++){
            mesh.verts[i].set_v(p_prime[i]);
            handle_objs[i].position.set(p_prime[i].x, p_prime[i].y, p_prime[i].z);
        }

        viewing_object.geometry.dispose();
        const new_positions = new THREE.BufferAttribute(new Float32Array(
            mesh.verts.map(v => v.to_THREE()).flat()
        ), 3);
        viewing_object.geometry.setAttribute('position', new_positions);

        
        if(use_QR) {
            // @ts-ignore
            if(L.delete && qr.delete) memoryManager.deleteExcept([L, qr]);
            // @ts-ignore
            else if(L.delete) memoryManager.deleteExcept([L]);
        } else {
            // @ts-ignore
            if(L.delete && llt.delete) memoryManager.deleteExcept([L, llt]);
            // @ts-ignore
            else if(L.delete) memoryManager.deleteExcept([L]);
        }
    }

    const handle_objs = mesh.verts.map(v => {
        const handle_geo = new THREE.SphereGeometry(0.02, 5, 5);
        const handle = new THREE.Mesh(handle_geo, new THREE.MeshBasicMaterial({
            color: 0xff000
        }));
        handle.position.set(v.x, v.y, v.z);
        scene.add(handle);
        return handle;
    });

    calculate_L();
    calculate_b();

    // handle dragging handles
    fixed_points.forEach(fp => handle_objs[fp.index].material.color.set(0xff0000));
    const drag_controls = new DragControls(handle_objs, camera, renderer.domElement);
    drag_controls.addEventListener('dragstart', (event) => {
        controls.enabled = false;
    });

    let waiting = false;
    const on_drag = (event: THREE.Event, callback: any) => {
        const v_index = handle_objs.findIndex(obj => event.object == obj);
        const handle = handle_objs[v_index];

        const new_pos = new V3(handle.position.x, handle.position.y, handle.position.z);
        const previous_index = fixed_points.findIndex(fp => fp.index == v_index);
        let rebuild_L = false;
        if(previous_index !== -1){
            fixed_points[previous_index].p_fixed = new_pos;
        } else {
            fixed_points.push({
                index: v_index,
                p_fixed: new_pos
            });
            rebuild_L = true;
        }
        fixed_indices = fixed_points.map(fixed_data => fixed_data.index); 
        fixed_points.forEach(fp => handle_objs[fp.index].material.color.set(0xff0000));
        
        if(rebuild_L) calculate_L();
        calculate_b();
        callback();
    }
    drag_controls.addEventListener('drag', (event) => {
        if(waiting) return;
        setTimeout(() => {
            waiting = true;
            on_drag(event, () => waiting = false);
        }, 10);
    });

    drag_controls.addEventListener('dragend', (event) => {
        controls.enabled = true;
    });
});
