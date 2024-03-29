import { OBJData } from "../util/data-structs";
import { HalfedgeMesh } from "./halfedge_mesh";

export const halfedge_mesh_to_index_triangle = (mesh: HalfedgeMesh): OBJData => {
    const vertices = mesh.verts;
    const indices: any[] = [];

    for(let i = 0; i < vertices.length; i++){
        let v = vertices[i];
        v.cache1 = i;
    }

    for(let f of mesh.faces){
        let v0 = f.he.vert;
        let v1 = f.he.next.vert;
        let v2 = f.he.next.next.vert;

        let i0 = v0.cache1;
        let i1 = v1.cache1;
        let i2 = v2.cache1;

        if(i0 === -1 || i1 === -1 || i2 === -1) throw Error("vertex not found");
        indices.push(i0, i1, i2);
    }

    return {
        vertices: vertices,
        indices: indices
    } 
}