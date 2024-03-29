import { Face } from "./atom_face";
import { Halfedge } from "./atom_halfedge";
import { V3 } from "./atom_vertex";

export class HalfedgeMesh {
    halfedges: Halfedge[];
    
    verts!: V3[];
    faces!: Face[];

    constructor(halfedges: Halfedge[]){
        this.halfedges = halfedges;
    }

    get_halfedge_snapshot(){
        return this.halfedges.slice(0);
    }

    get_vertex_snapshot(){
        return this.verts.slice(0);
    }

    set_verts(verts: V3[]){
        this.verts = verts;
    }

    set_faces(faces: Face[]){
        this.faces = faces;
    }

    cull_old_elements(){
        this.faces = this.faces.filter(f => !f.to_delete);
        this.halfedges = this.halfedges.filter(he => !he.to_delete);
        this.verts = this.verts.filter(v => !v.to_delete);
    }

    reset_halfedge_flags_n(n: number){
        if(n === 1) this.halfedges.forEach(he => he.flag1 = false);
        else if(n === 2) this.halfedges.forEach(he => he.flag2 = false);
        else if(n === 3) this.halfedges.forEach(he => he.flag3 = false);
        else if(n === 4) this.halfedges.forEach(he => he.flag4 = false);
    }


    reset_flags(){
        this.reset_vert_flags();
        this.reset_halfedge_flags();
    }

    reset_vert_flags(){
        this.verts.forEach(v => v.reset_flags());
    }

    reset_halfedge_flags(){
        this.halfedges.forEach(he => he.reset_flags());
    }
}