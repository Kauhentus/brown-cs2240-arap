import * as THREE from "three";
import { scene } from "./init_three";
import { V3 } from "geo/atom_vertex";

export const view_fixedtriangles = (triangles: V3[][], color: number) => {
    let vertices: number[] = triangles.map(t => [
        t[0], t[1], t[1], t[2], t[2], t[0]
    ]).flat().map(v => v.to_THREE()).flat();

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(vertices), 3));
    const lines = new THREE.LineSegments(geometry, new THREE.LineBasicMaterial({color: color}));
    scene.add(lines);

    return {
        objects: [lines]
    }
}