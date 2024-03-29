import { Face } from "./atom_face";
import { V3 } from "./atom_vertex";

export class Halfedge {
    twin!: Halfedge;
    next!: Halfedge;

    vert!: V3;
    face!: Face;

    to_delete: boolean;

    flag1: boolean = false;
    flag2: boolean = false;
    flag3: boolean = false;
    flag4: boolean = false;

    cache1: any;
    cache2: any;
    cache3: any;
    cache4: any;

    constructor(v: V3){
        this.vert = v;
        this.to_delete = false;
    }

    set_twin(twin: Halfedge){
        this.twin = twin;
    }

    set_next(next: Halfedge){
        this.next = next;
    }
    
    set_vert(v: V3){
        this.vert = v;
    }

    set_face(f: Face){
        this.face = f;
    }

    reset_flags(){
        this.flag1 = false;
        this.flag2 = false;
        this.flag3 = false;
        this.flag4 = false;
    }

    reset_cache(){
        this.cache1 = undefined;
        this.cache2 = undefined;
        this.cache3 = undefined;
        this.cache4 = undefined;
    }
}