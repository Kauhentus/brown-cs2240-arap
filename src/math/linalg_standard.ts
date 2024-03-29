import { V3 } from "../geo/atom_vertex";

export const zero_v = () => new V3(0, 0, 0);

export const add_v = (a: V3, b: V3) => {
    return new V3(
        a.x + b.x,
        a.y + b.y,
        a.z + b.z
    )
}

export const add_vs = (...vertices: V3[]) => {
    return vertices.reduce((a, v) => add_v(a, v), new V3(0, 0, 0));
}

export const sub_v = (a: V3, b: V3) => {
    return new V3(
        a.x - b.x,
        a.y - b.y,
        a.z - b.z
    )
}

export const smul_v = (a: V3, s: number) => {
    return new V3(
        a.x * s,
        a.y * s,
        a.z * s
    )
}

export const smul_vec3 = (a: number[], s: number): [number, number, number] => {
    return [
        a[0] * s,
        a[1] * s,
        a[2] * s
    ];
}

export const sdiv_v = (a: V3, s: number) => {
    return new V3(
        a.x / s,
        a.y / s,
        a.z / s
    )
}

export const same_v = (a: V3, b: V3) => {
    const delta = squared_distance_v(a, b);
    return delta < 1e-4;
}

export const avg_v = (...vs: V3[]) => {
    return new V3(
        vs.reduce((a, v) => a + v.x, 0) / vs.length,
        vs.reduce((a, v) => a + v.y, 0) / vs.length,
        vs.reduce((a, v) => a + v.z, 0) / vs.length,
    );
}

export const lerp_v = (t: number, a: V3, b: V3) => {
    return new V3(
        t * a.x + (1 - t) * b.x,
        t * a.y + (1 - t) * b.y,
        t * a.z + (1 - t) * b.z,
    );
}

export const min_v = (a: V3, b: V3) => {
    return new V3(
        Math.min(a.x, b.x),
        Math.min(a.y, b.y),
        Math.min(a.z, b.z)
    );
}

export const max_v = (a: V3, b: V3) => {
    return new V3(
        Math.max(a.x, b.x),
        Math.max(a.y, b.y),
        Math.max(a.z, b.z)
    );
}

export const max_component_v = (a: V3) => {
    return Math.max(a.x, a.y, a.z);
}

export const abs_v = (a: V3) => {
    return new V3(
        Math.abs(a.x),
        Math.abs(a.y),
        Math.abs(a.z)
    );
}


export const distance_v = (a: V3, b: V3) => {
    return Math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2);
}

export const squared_distance_v = (a: V3, b: V3) => {
    return (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2;
}

export const magnitude_v = (a: V3) => {
    return Math.sqrt(a.x ** 2 + a.y ** 2 + a.z ** 2);
}

export const magnitude_vec3 = (a: number[]) => {
    return Math.sqrt(a[0] ** 2 + a[1] ** 2 + a[2] ** 2);
}

export const normalize_v = (a: V3) => {
    const length = magnitude_v(a);
    return new V3(
        a.x / length,
        a.y / length,
        a.z / length
    )
}

export const dot_v = (a: V3, b: V3) => {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

export const dot_vec3 = (a: number[], b: number[]) => {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

export const dot_vec4 = (a: number[], b: number[]) => {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}

export const cross_v = (a: V3, b: V3) => {
    return new V3(
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    );
}