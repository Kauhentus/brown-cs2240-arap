import { Halfedge } from "./atom_halfedge";
import { V3 } from "./atom_vertex";

let max_vertices_per_vertex = 1024;

export type V3NeighborData = {
    verts: V3[];
    halfedges: Halfedge[];
    leaving_halfedges: Halfedge[];
    entering_halfedges: Halfedge[];
}

export const vertex_get_neighbors = (v: V3): V3NeighborData => {
    let verts: V3[] = [];
    let halfedges: Halfedge[] = [];
    let leaving_halfedges: Halfedge[] = [];
    let entering_halfedges: Halfedge[] = [];

    let start_he = v.he;
    let cur_he = start_he;
    let counter = 0;
    do {
        halfedges.push(cur_he);
        leaving_halfedges.push(cur_he);
        
        let twin = cur_he.twin;
        halfedges.push(cur_he);
        entering_halfedges.push(twin);
        verts.push(twin.vert);

        // view_vertex([twin.vert], rgb_to_hex(255 - counter * 30, 100, counter * 30))

        cur_he = twin.next;
        counter += 1;
    } while(cur_he !== start_he && counter < max_vertices_per_vertex);

    if(counter === max_vertices_per_vertex) throw Error("vertex fan exceeded 1024 edges");

    return {
        verts: verts,
        halfedges: halfedges,
        entering_halfedges: entering_halfedges,
        leaving_halfedges: leaving_halfedges
    }
}