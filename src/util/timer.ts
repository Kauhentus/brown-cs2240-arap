let time_keeper = window.performance.now();

export const reset_elapsed_time = () => {
    time_keeper = window.performance.now();
}

export const get_elapsed_time = (reset: boolean = false) => {
    let new_time = window.performance.now();
    let elapsed_time = new_time - time_keeper;
    if(reset) reset_elapsed_time();
    return elapsed_time.toPrecision(4);
}

let time_keeper2 = window.performance.now();

export const reset_elapsed_time2 = () => {
    time_keeper2 = window.performance.now();
}

export const get_elapsed_time2 = (reset: boolean = false) => {
    let new_time = window.performance.now();
    let elapsed_time = new_time - time_keeper2;
    if(reset) reset_elapsed_time2();
    return elapsed_time.toPrecision(4);
}

let time_keeper3 = window.performance.now();

export const reset_elapsed_time3 = () => {
    time_keeper3 = window.performance.now();
}

export const get_elapsed_time3 = (reset: boolean = false) => {
    let new_time = window.performance.now();
    let elapsed_time = new_time - time_keeper3;
    if(reset) reset_elapsed_time3();
    return parseFloat(elapsed_time.toPrecision(4));
}