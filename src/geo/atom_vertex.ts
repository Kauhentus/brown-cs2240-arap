import { Halfedge } from "./atom_halfedge";

let c = 0;

export class V3 {
    x!: number;
    y!: number;
    z!: number;

    he!: Halfedge;
    id: symbol;
    to_delete: boolean;

    flag1: boolean = false;
    flag2: boolean = false;
    flag3: boolean = false;
    flag4: boolean = false;

    cache1: any;
    cache2: any;
    cache3: any;
    cache4: any;


    constructor(x?: number , y?: number, z?: number){
        if(x !== undefined) this.x = x;
        if(y !== undefined) this.y = y;
        if(z !== undefined) this.z = z;

        this.id = Symbol(c++);
        this.to_delete = false;
    }

    static zero(){
        return new V3(0, 0, 0);
    }

    clone(){
        return new V3(this.x, this.y, this.z);
    }
    
    set_pos(x: number, y: number, z: number){
        this.x = x;
        this.y = y;
        this.z = z;
    }

    to_THREE(){
        return [this.x, this.y, this.z];
    }

    to_homogeneous(){
        return [this.x, this.y, this.z, 1];
    }

    add_xyz(x: number, y: number, z: number){
        this.x += x;
        this.y += y;
        this.z += z;
        return this;
    }

    add_v(v: V3){
        this.x += v.x;
        this.y += v.y;
        this.z += v.z;
        return this;
    }

    sub_xyz(x: number, y: number, z: number){
        this.x -= x;
        this.y -= y;
        this.z -= z;
        return this;
    }

    sub_v(v: V3){
        this.x -= v.x;
        this.y -= v.y;
        this.z -= v.z;
        return this;
    }

    set_xyz(x: number, y: number, z: number){
        this.x = x;
        this.y = y;
        this.z = z;
        return this;
    }

    set_v(v: V3){
        this.x = v.x;
        this.y = v.y;
        this.z = v.z;
        return this;
    }

    zero_out(){
        this.x = 0;
        this.y = 0;
        this.z = 0;
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