import * as THREE from "three";
import { OBJData } from "../util/data-structs";
import { scene } from "./init_three";

export const view_index_triangle = (data: OBJData, color: number) => {
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(
        data.vertices.map(v => v.to_THREE()).flat()
    ), 3));
    geometry.setIndex(data.indices);
    // geometry.computeVertexNormals();
    const material = new THREE.MeshStandardMaterial({
        color: color,
        flatShading: true
    });
    const mesh = new THREE.Mesh(geometry, material);

    scene.add(mesh);

    view_index_triangle_wireframe(data, 0xffaa00);
}

export const view_index_triangle_wireframe = (data: OBJData, color: number) => {
    const geometry = new THREE.BufferGeometry();

    let vertices: number[] = [];
    for(let i = 0; i < data.indices.length; i += 3){
        let i0 = data.indices[i];
        let i1 = data.indices[i + 1];
        let i2 = data.indices[i + 2];
        let v0 = data.vertices[i0].to_THREE();
        let v1 = data.vertices[i1].to_THREE();
        let v2 = data.vertices[i2].to_THREE();
        vertices.push(...v0, ...v1);
        vertices.push(...v1, ...v2);
        vertices.push(...v2, ...v0);
    }

    geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(vertices), 3));
    const lines = new THREE.LineSegments(geometry, new THREE.LineBasicMaterial({color: color}));
    scene.add(lines);
}