import { V3 } from "../geo/atom_vertex";
import { cross_v, dot_v, normalize_v, smul_v, sub_v } from "./linalg_standard";

// matrix creation

export const mat3_identity = () => {
    return [
        1, 0, 0,
        0, 1, 0,
        0, 0, 1
    ];
}

export const mat4_identity = () => {
    return [
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    ];
}

export const mat3_zeroes = () => {
    return [
        0, 0, 0,
        0, 0, 0,
        0, 0, 0
    ];
}

export const mat4_zeroes = () => {
    return [
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0
    ];
}

export const mat3_clone = (mat: number[]) => {
    return mat.slice(0);
}

export const mat4_clone = (mat: number[]) => {
    return mat.slice(0);
}

export const mat3_rot_x = (theta: number) => {
    return [
        1, 0, 0,
        0, Math.cos(theta), Math.sin(theta),
        0, -Math.sin(theta), Math.cos(theta),
    ];
}

export const mat4_rot_x = (theta: number) => {
    return [
        1, 0, 0, 0,
        0, Math.cos(theta), Math.sin(theta), 0,
        0, -Math.sin(theta), Math.cos(theta), 0,
        0, 0, 0, 1
    ];
}

export const mat3_rot_y = (theta: number) => {
    return [
        Math.cos(theta), 0, Math.sin(theta),
        0, 1, 0,
        -Math.sin(theta), 0, Math.cos(theta)
    ];
}

export const mat4_rot_y = (theta: number) => {
    return [
        Math.cos(theta), 0, Math.sin(theta), 0,
        0, 1, 0, 0,
        -Math.sin(theta), 0, Math.cos(theta), 0,
        0, 0, 0, 1
    ];
}

export const mat3_rot_z = (theta: number) => {
    return [
        Math.cos(theta), -Math.sin(theta), 0,
        Math.sin(theta), Math.cos(theta), 0,
        0, 0, 1
    ];
}

export const mat4_rot_z = (theta: number) => {
    return [
        Math.cos(theta), -Math.sin(theta), 0, 0,
        Math.sin(theta), Math.cos(theta), 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    ];
}

// matrix operations

export const add_mat3 = (a: number[], b: number[]) => {
    return [
        a[0] + b[0], a[1] + b[1], a[2] + b[2],
        a[3] + b[3], a[4] + b[4], a[5] + b[5],
        a[6] + b[6], a[7] + b[7], a[8] + b[8]
    ];
}

export const add_mat4 = (a: number[], b: number[]) => {
    return [
        a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3],
        a[4] + b[4], a[5] + b[5], a[6] + b[6], a[7] + b[7],
        a[8] + b[8], a[9] + b[9], a[10] + b[10], a[11] + b[11],
        a[12] + b[12], a[13] + b[13], a[14] + b[14], a[15] + b[15]
    ];
}

export const sub_mat3 = (a: number[], b: number[]) => {
    return [
        a[0] - b[0], a[1] - b[1], a[2] - b[2],
        a[3] - b[3], a[4] - b[4], a[5] - b[5],
        a[6] - b[6], a[7] - b[7], a[8] - b[8]
    ];
}

export const sub_mat4 = (a: number[], b: number[]) => {
    return [
        a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3],
        a[4] - b[4], a[5] - b[5], a[6] - b[6], a[7] - b[7],
        a[8] - b[8], a[9] - b[9], a[10] - b[10], a[11] - b[11],
        a[12] - b[12], a[13] - b[13], a[14] - b[14], a[15] - b[15]
    ];
}

export const smul_mat3 = (a: number[], s: number) => {
    return [
        a[0] * s, a[1] * s, a[2] * s,
        a[3] * s, a[4] * s, a[5] * s,
        a[6] * s, a[7] * s, a[8] * s
    ];
}

export const smul_mat4 = (a: number[], s: number) => {
    return [
        a[0] * s, a[1] * s, a[2] * s, a[3] * s,
        a[4] * s, a[5] * s, a[6] * s, a[7] * s,
        a[8] * s, a[9] * s, a[10] * s, a[11] * s,
        a[12] * s, a[13] * s, a[14] * s, a[15] * s
    ];
}

export const mat3_matmul = (a: number[], b: number[]) => {
    return [
        // col 1
        a[0]*b[0] + a[3]*b[1] + a[6]*b[2],
        a[1]*b[0] + a[4]*b[1] + a[7]*b[2],
        a[2]*b[0] + a[5]*b[1] + a[8]*b[2],

        // col 2
        a[0]*b[3] + a[3]*b[4] + a[6]*b[5],
        a[1]*b[3] + a[4]*b[4] + a[7]*b[5],
        a[2]*b[3] + a[5]*b[4] + a[8]*b[5],

        // col 3
        a[0]*b[6] + a[3]*b[7] + a[6]*b[8],
        a[1]*b[6] + a[4]*b[7] + a[7]*b[8],
        a[2]*b[6] + a[5]*b[7] + a[8]*b[8],
    ];
}

export const mat4_matmul = (a: number[], b: number[]) => {
    return [
        // col 1
        a[0]*b[0] + a[4]*b[1] + a[8]*b[2] + a[12]*b[3],
        a[1]*b[0] + a[5]*b[1] + a[9]*b[2] + a[13]*b[3],
        a[2]*b[0] + a[6]*b[1] + a[10]*b[2] + a[14]*b[3],
        a[3]*b[0] + a[7]*b[1] + a[11]*b[2] + a[15]*b[3],

        // col 2
        a[0]*b[4] + a[4]*b[5] + a[8]*b[6] + a[12]*b[7],
        a[1]*b[4] + a[5]*b[5] + a[9]*b[6] + a[13]*b[7],
        a[2]*b[4] + a[6]*b[5] + a[10]*b[6] + a[14]*b[7],
        a[3]*b[4] + a[7]*b[5] + a[11]*b[6] + a[15]*b[7],

        // col 3
        a[0]*b[8] + a[4]*b[9] + a[8]*b[10] + a[12]*b[11],
        a[1]*b[8] + a[5]*b[9] + a[9]*b[10] + a[13]*b[11],
        a[2]*b[8] + a[6]*b[9] + a[10]*b[10] + a[14]*b[11],
        a[3]*b[8] + a[7]*b[9] + a[11]*b[10] + a[15]*b[11],

        // col 4
        a[0]*b[12] + a[4]*b[13] + a[8]*b[14] + a[12]*b[15],
        a[1]*b[12] + a[5]*b[13] + a[9]*b[14] + a[13]*b[15],
        a[2]*b[12] + a[6]*b[13] + a[10]*b[14] + a[14]*b[15],
        a[3]*b[12] + a[7]*b[13] + a[11]*b[14] + a[15]*b[15]
    ];
}

export const mat3_vecmul = (mat: number[], vec3: number[]) => {
    return [
        mat[0]*vec3[0] + mat[3]*vec3[1] + mat[6]*vec3[2],
        mat[1]*vec3[0] + mat[4]*vec3[1] + mat[7]*vec3[2],
        mat[2]*vec3[0] + mat[5]*vec3[1] + mat[8]*vec3[2],
    ]
}

export const mat4_vecmul = (mat: number[], vec4: number[]) => {
    return [
        mat[0]*vec4[0] + mat[4]*vec4[1] + mat[8]*vec4[2] + mat[12]*vec4[3],
        mat[1]*vec4[0] + mat[5]*vec4[1] + mat[9]*vec4[2] + mat[13]*vec4[3],
        mat[2]*vec4[0] + mat[6]*vec4[1] + mat[10]*vec4[2] + mat[14]*vec4[3],
        mat[3]*vec4[0] + mat[7]*vec4[1] + mat[11]*vec4[2] + mat[15]*vec4[3],
    ]
}

export const vec4_dot = (a: number[], b: number[]) => {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

export const mat3_applyv = (mat: number[], v: V3) => {
    let r = mat3_vecmul(mat, [v.x, v.y, v.z]);
    return new V3(r[0], r[1], r[2]);
}

export const mat4_applyv3 = (mat: number[], v: V3, w_one: boolean, homogenize: boolean) => {
    let r = mat4_vecmul(mat, [v.x, v.y, v.z, w_one ? 1 : 0]);
    if(homogenize){
        return new V3(r[0] / r[3], r[1] / r[3], r[2] / r[3]);
    } else {
        return new V3(r[0], r[1], r[2]);
    }
}

export const mat3_trace = (a: number[]) => {
    return a[0] + a[4] + a[8];
}

export const mat4_trace = (a: number[]) => {
    return a[0] + a[5] + a[10] + a[15];
}

export const mat3_invert = (m: number[]) => {
    // 0  3  6
    // 1  4  7  
    // 2  5  8

    const det = 
        m[0] * (m[4] * m[8] - m[7] * m[5]) -
        m[1] * (m[3] * m[8] - m[5] * m[6]) +
        m[2] * (m[3] * m[7] - m[4] * m[6]);
    if (det == 0) console.error("non-invertible matrix encountered");
    const invdet = 1 / det;
    
    const minv = new Array(9);
    minv[0] = (m[4] * m[8] - m[7] * m[5]) * invdet;
    minv[1] = (m[2] * m[7] - m[1] * m[8]) * invdet;
    minv[2] = (m[1] * m[5] - m[2] * m[4]) * invdet;
    minv[3] = (m[5] * m[6] - m[3] * m[8]) * invdet;
    minv[4] = (m[0] * m[8] - m[2] * m[6]) * invdet;
    minv[5] = (m[3] * m[2] - m[0] * m[5]) * invdet;
    minv[6] = (m[3] * m[7] - m[6] * m[4]) * invdet;
    minv[7] = (m[6] * m[1] - m[0] * m[7]) * invdet;
    minv[8] = (m[0] * m[4] - m[3] * m[1]) * invdet;
    return minv;
}

export const mat3_det = (m: number[]) => {
    const det = 
        m[0] * (m[4] * m[8] - m[7] * m[5]) -
        m[1] * (m[3] * m[8] - m[5] * m[6]) +
        m[2] * (m[3] * m[7] - m[4] * m[6]);

    return det;
}


export const mat4_invert = (m: number[]) => {
    // 0  4  8  12
    // 1  5  9  13
    // 2  6  10 14
    // 3  7  11 15

    const inv = new Array(16);
    const invOut = new Array(16);

    inv[0] = m[5]  * m[10] * m[15] - 
                m[5]  * m[11] * m[14] - 
                m[9]  * m[6]  * m[15] + 
                m[9]  * m[7]  * m[14] +
                m[13] * m[6]  * m[11] - 
                m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
                m[4]  * m[11] * m[14] + 
                m[8]  * m[6]  * m[15] - 
                m[8]  * m[7]  * m[14] - 
                m[12] * m[6]  * m[11] + 
                m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
                m[4]  * m[11] * m[13] - 
                m[8]  * m[5] * m[15] + 
                m[8]  * m[7] * m[13] + 
                m[12] * m[5] * m[11] - 
                m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
                m[4]  * m[10] * m[13] +
                m[8]  * m[5] * m[14] - 
                m[8]  * m[6] * m[13] - 
                m[12] * m[5] * m[10] + 
                m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
                m[1]  * m[11] * m[14] + 
                m[9]  * m[2] * m[15] - 
                m[9]  * m[3] * m[14] - 
                m[13] * m[2] * m[11] + 
                m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
                m[0]  * m[11] * m[14] - 
                m[8]  * m[2] * m[15] + 
                m[8]  * m[3] * m[14] + 
                m[12] * m[2] * m[11] - 
                m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
                m[0]  * m[11] * m[13] + 
                m[8]  * m[1] * m[15] - 
                m[8]  * m[3] * m[13] - 
                m[12] * m[1] * m[11] + 
                m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
                m[0]  * m[10] * m[13] - 
                m[8]  * m[1] * m[14] + 
                m[8]  * m[2] * m[13] + 
                m[12] * m[1] * m[10] - 
                m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
                m[1]  * m[7] * m[14] - 
                m[5]  * m[2] * m[15] + 
                m[5]  * m[3] * m[14] + 
                m[13] * m[2] * m[7] - 
                m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
                m[0]  * m[7] * m[14] + 
                m[4]  * m[2] * m[15] - 
                m[4]  * m[3] * m[14] - 
                m[12] * m[2] * m[7] + 
                m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
                m[0]  * m[7] * m[13] - 
                m[4]  * m[1] * m[15] + 
                m[4]  * m[3] * m[13] + 
                m[12] * m[1] * m[7] - 
                m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
                m[0]  * m[6] * m[13] + 
                m[4]  * m[1] * m[14] - 
                m[4]  * m[2] * m[13] - 
                m[12] * m[1] * m[6] + 
                m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
                m[1] * m[7] * m[10] + 
                m[5] * m[2] * m[11] - 
                m[5] * m[3] * m[10] - 
                m[9] * m[2] * m[7] + 
                m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
                m[0] * m[7] * m[10] - 
                m[4] * m[2] * m[11] + 
                m[4] * m[3] * m[10] + 
                m[8] * m[2] * m[7] - 
                m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
                m[0] * m[7] * m[9] + 
                m[4] * m[1] * m[11] - 
                m[4] * m[3] * m[9] - 
                m[8] * m[1] * m[7] + 
                m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
                m[0] * m[6] * m[9] - 
                m[4] * m[1] * m[10] + 
                m[4] * m[2] * m[9] + 
                m[8] * m[1] * m[6] - 
                m[8] * m[2] * m[5];

    let det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
    if (det == 0) console.error("non-invertible matrix encountered");
    det = 1.0 / det;
    for (let i = 0; i < 16; i++) invOut[i] = inv[i] * det;
    return invOut;
}

export const mat4_det = (m: number[]) => {
    const inv_0 = m[5]  * m[10] * m[15] - 
                m[5]  * m[11] * m[14] - 
                m[9]  * m[6]  * m[15] + 
                m[9]  * m[7]  * m[14] +
                m[13] * m[6]  * m[11] - 
                m[13] * m[7]  * m[10];

    const inv_4 = -m[4]  * m[10] * m[15] + 
                m[4]  * m[11] * m[14] + 
                m[8]  * m[6]  * m[15] - 
                m[8]  * m[7]  * m[14] - 
                m[12] * m[6]  * m[11] + 
                m[12] * m[7]  * m[10];

    const inv_8 = m[4]  * m[9] * m[15] - 
                m[4]  * m[11] * m[13] - 
                m[8]  * m[5] * m[15] + 
                m[8]  * m[7] * m[13] + 
                m[12] * m[5] * m[11] - 
                m[12] * m[7] * m[9];

    const inv_12 = -m[4]  * m[9] * m[14] + 
                m[4]  * m[10] * m[13] +
                m[8]  * m[5] * m[14] - 
                m[8]  * m[6] * m[13] - 
                m[12] * m[5] * m[10] + 
                m[12] * m[6] * m[9];

    let det = m[0] * inv_0 + m[1] * inv_4 + m[2] * inv_8 + m[3] * inv_12;

    return det;
}


export const mat4_invert_with_check = (m: number[]) => {
    // 0  4  8  12
    // 1  5  9  13
    // 2  6  10 14
    // 3  7  11 15

    const inv = new Array(16);
    const invOut = new Array(16);

    inv[0] = m[5]  * m[10] * m[15] - 
                m[5]  * m[11] * m[14] - 
                m[9]  * m[6]  * m[15] + 
                m[9]  * m[7]  * m[14] +
                m[13] * m[6]  * m[11] - 
                m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
                m[4]  * m[11] * m[14] + 
                m[8]  * m[6]  * m[15] - 
                m[8]  * m[7]  * m[14] - 
                m[12] * m[6]  * m[11] + 
                m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
                m[4]  * m[11] * m[13] - 
                m[8]  * m[5] * m[15] + 
                m[8]  * m[7] * m[13] + 
                m[12] * m[5] * m[11] - 
                m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
                m[4]  * m[10] * m[13] +
                m[8]  * m[5] * m[14] - 
                m[8]  * m[6] * m[13] - 
                m[12] * m[5] * m[10] + 
                m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
                m[1]  * m[11] * m[14] + 
                m[9]  * m[2] * m[15] - 
                m[9]  * m[3] * m[14] - 
                m[13] * m[2] * m[11] + 
                m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
                m[0]  * m[11] * m[14] - 
                m[8]  * m[2] * m[15] + 
                m[8]  * m[3] * m[14] + 
                m[12] * m[2] * m[11] - 
                m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
                m[0]  * m[11] * m[13] + 
                m[8]  * m[1] * m[15] - 
                m[8]  * m[3] * m[13] - 
                m[12] * m[1] * m[11] + 
                m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
                m[0]  * m[10] * m[13] - 
                m[8]  * m[1] * m[14] + 
                m[8]  * m[2] * m[13] + 
                m[12] * m[1] * m[10] - 
                m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
                m[1]  * m[7] * m[14] - 
                m[5]  * m[2] * m[15] + 
                m[5]  * m[3] * m[14] + 
                m[13] * m[2] * m[7] - 
                m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
                m[0]  * m[7] * m[14] + 
                m[4]  * m[2] * m[15] - 
                m[4]  * m[3] * m[14] - 
                m[12] * m[2] * m[7] + 
                m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
                m[0]  * m[7] * m[13] - 
                m[4]  * m[1] * m[15] + 
                m[4]  * m[3] * m[13] + 
                m[12] * m[1] * m[7] - 
                m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
                m[0]  * m[6] * m[13] + 
                m[4]  * m[1] * m[14] - 
                m[4]  * m[2] * m[13] - 
                m[12] * m[1] * m[6] + 
                m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
                m[1] * m[7] * m[10] + 
                m[5] * m[2] * m[11] - 
                m[5] * m[3] * m[10] - 
                m[9] * m[2] * m[7] + 
                m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
                m[0] * m[7] * m[10] - 
                m[4] * m[2] * m[11] + 
                m[4] * m[3] * m[10] + 
                m[8] * m[2] * m[7] - 
                m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
                m[0] * m[7] * m[9] + 
                m[4] * m[1] * m[11] - 
                m[4] * m[3] * m[9] - 
                m[8] * m[1] * m[7] + 
                m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
                m[0] * m[6] * m[9] - 
                m[4] * m[1] * m[10] + 
                m[4] * m[2] * m[9] + 
                m[8] * m[1] * m[6] - 
                m[8] * m[2] * m[5];

    let det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
    if (det == 0) {
        return {
            invertible: false,
            result: mat4_zeroes()
        };
    }
    det = 1.0 / det;
    for (let i = 0; i < 16; i++) invOut[i] = inv[i] * det;
    return {
        invertible: true,
        result: invOut
    };;
}

export const mat3_tranpose = (m: number[]) => {
    return [
        m[0], m[3], m[6],
        m[1], m[4], m[7],
        m[2], m[5], m[8]
    ]
}

export const mat4_tranpose = (m: number[]) => {
    return [
        m[0], m[4], m[8], m[12],
        m[1], m[5], m[9], m[13],
        m[2], m[6], m[10], m[14],
        m[3], m[7], m[11], m[15],
    ]
}

export const mat4_translate = (x: number, y: number, z: number) => {
    return [
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0, 
        x, y, z, 1
    ]
}

export const mat3_scale = (sx: number, sy: number, sz: number) => {
    return [
        sx, 0, 0,
        0, sy, 0,
        0, 0, sz
    ]
}

export const mat4_scale = (sx: number, sy: number, sz: number) => {
    return [
        sx, 0, 0, 0,
        0, sy, 0, 0,
        0, 0, sz, 0, 
        0, 0, 0, 1
    ]
}

// misc utilities

export const mat4_to_mat3 = (m: number[]) => {
    return [
        m[0], m[1], m[2],
        m[4], m[5], m[6],
        m[8], m[9], m[10]
    ]
}

export const mat3_to_mat4 = (m: number[], w_one: boolean = true) => {
    return [
        m[0], m[1], m[2], 0,
        m[4], m[5], m[6], 0,
        m[8], m[9], m[10], 0,
        0, 0, 0, w_one ? 1 : 0,
    ]
}

export const world_to_camera = (position: V3, look: V3, up: V3) => {
    const w = smul_v(normalize_v(look), -1);
    const v = normalize_v(sub_v(up, smul_v(w, dot_v(up, w))));
    const u = cross_v(v, w);

    const mat_t = mat4_translate(-position.x, -position.y, -position.z);
    const mat_r = [
        u.x, v.x, w.x, 0,
        u.y, v.y, w.y, 0,
        u.z, v.z, w.z, 0, 
        0, 0, 0, 1
    ];

    const mat_v = mat4_matmul(mat_r, mat_t);
    return mat_v;
}