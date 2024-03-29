export const integrator_euler_step = (x0: number[][], delta_x: number[][], step: number) => {
    let result = x0.map(segment => segment.slice(0));
    let num_segments = x0.length;
    
    for(let j = 0; j < num_segments; j++){
        let segment_len = x0[j].length;
        for(let i = 0; i < segment_len; i++){
            result[j][i] = x0[j][i] + step * delta_x[j][i];
        }
    }

    return {
        new_x0: result
    }
}