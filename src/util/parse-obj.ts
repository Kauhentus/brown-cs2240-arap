import { V3 } from "../geo/atom_vertex";
import { OBJData } from "./data-structs";

export const parse_obj = (obj_data: string): OBJData => {
    // first parse the obj data
    const raw_obj_lines = obj_data.split('\n');
    const object_groups = {
        vertices: [] as number[],
        indices: [] as number[]
    };

    for(let raw_line of raw_obj_lines){
        let line = raw_line.replace(/\s+/g, ' ').replace(/#.*$/, '').trim();
        if(line.length === 0) continue;
        if(line[0] === '#') continue;

        else if(line.slice(0, 2) === 'v '){
            let data = line.slice(2).trim().split(' ').map(parseFloat);
            object_groups.vertices.push(...data.slice(0, 3));
        }

        else if(line.slice(0, 2) === 'f '){
            let num_vertices = object_groups.vertices.length / 3;
            let raw_indices = line.slice(2).trim().split(' ').map(triplet => {
                let i = parseInt(triplet.split('/')[0]);
                if(i > 0) return i;
                else return num_vertices + i + 1;
            }).flat();

            if(raw_indices.length == 3){
                object_groups.indices.push(...raw_indices.map(i => i - 1))
            }  else throw Error("4+ sides encountered")
        }
    }

    const indices = object_groups.indices;
    const vertices: V3[] = [];

    for(let i = 0; i < object_groups.vertices.length; i += 3){
        vertices.push(new V3(
            object_groups.vertices[i],
            object_groups.vertices[i + 1],
            object_groups.vertices[i + 2]
        ));
    }
    const final_result: OBJData = {
        indices: indices,
        vertices: vertices
    }

    return final_result;
}