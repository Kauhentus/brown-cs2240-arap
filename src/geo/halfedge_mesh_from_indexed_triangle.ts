import { Face } from "./atom_face";
import { Halfedge } from "./atom_halfedge";
import { V3 } from "./atom_vertex";
import { HalfedgeMesh } from "./halfedge_mesh";

export const indexed_triangle_to_halfedge_mesh = (vertices: V3[], indices: number[]): HalfedgeMesh => {
    let halfedges: Halfedge[] = [];
    let faces: Face[] = [];

    let halfedge_map : {[key: symbol]: Halfedge[]} = {};
    vertices.forEach(v => halfedge_map[v.id] = []);

    for(let i = 0; i < indices.length; i += 3){
        let i0 = indices[i];
        let i1 = indices[i + 1];
        let i2 = indices[i + 2];

        let v0 = vertices[i0]; let v1 = vertices[i1]; let v2 = vertices[i2];
        let f = new Face();

        let he0 = new Halfedge(v0); let he1 = new Halfedge(v1); let he2 = new Halfedge(v2);
        he0.next = he1; he1.next = he2; he2.next = he0;
        he0.vert = v0; he1.vert = v1; he2.vert = v2; 
        v0.he = he0; v1.he = he1; v2.he = he2;

        halfedge_map[v0.id].push(he0); 
        halfedge_map[v1.id].push(he1); 
        halfedge_map[v2.id].push(he2); 
        
        f.he = he0;
        he0.face = f; he1.face = f; he2.face = f; 

        halfedges.push(he0, he1, he2);
        faces.push(f);
    }
    
    for(let he of halfedges){
        let next_vert = he.next.vert;
        let potential_twins = halfedge_map[next_vert.id];
        let twin = potential_twins.find(he2 => he2.next.vert === he.vert);
        if(twin === undefined) throw Error("could not find twin");
        he.twin = twin;
    }

    let halfedge_mesh = new HalfedgeMesh(halfedges);
    halfedge_mesh.verts = vertices;
    halfedge_mesh.faces = faces;

    return halfedge_mesh;
}