import * as THREE from "three";
import { scene } from "./init_three";
import { V3 } from "../geo/atom_vertex";

export const view_points = (vertices: V3[], color: number, point_size: number = 0.1) => {
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(vertices.map(v => v.to_THREE()).flat()), 3));
    const points_object = new THREE.Points(geometry, new THREE.PointsMaterial({
        color: color,
        size: point_size
    }));
    scene.add(points_object);

    return points_object;
}