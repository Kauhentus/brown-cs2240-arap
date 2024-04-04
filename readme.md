# ARAP project

## Running the program

1. Install dependencies with `npm i`
2. Build the program with `npx webpack` (may need to install `webpack` and `webpack-cli`)
3. Run the program on a local server with `http-server -c-1` (may need to install `http-server`)

To load a specific mesh, change the mesh URL in line 36 of index.ts. This could have been an HTML drop down menu but I got lazy...

The program starts with one vertex fixed already, highlighted in red. Drag on other vertices to fix them to.

To unfix a vertex, hold down space and drag on it slightly. Setup is a bit scuffed but works nevertheless.

## Features implemented

My code attempts to implement all the basic requirements. It:

* Builds the L matrix using the neighbors and cotangent weights
* Applies user constraints by zero-ing out rows and columns of L
* Precomputes the decomposition of L 

* Determines the best rotation R's through SVD
* Optimizes the best positions p_i through solving the sparse system

Note: I used the `geometry-processing-js` library which contains Eigen transpiled to `asm.js`. For the Cholesky, it precomputes the decomposition so it only needs backpropagate values. I also tried using the QR decomposition
which seemed to give more stable results and again, the decomposition is precomputed.

## Bugs

1) I messed up somewhere in either setting up the L matrix or b vector (I think it's the b vector with zero-ing out fixed points since my matrix is confirmed symmetric positive definite). My meshes when loaded have small visual glitches which compound, the larger the mesh is. This is only really obvious with the teapot mesh. The glitches also compound when I repeatedly make permanent deformations via un-anchoring previous points.

2) Compounding from 1), the mesh does't deform exactly according to as-rigidly-as-possible. For instance, the 2 points that should correspond with rotation don't actually rotate, I get some intermediate deformation instead. The fact that I only run 3 iterations of the iterative solving doesn't help, but any more and the program doesn't run in realtime.

# Video demos

#### Anchor 1 point translation

#### Anchor 2 point rotation (bug)

#### Permanent deformations (small bug)

#### Waving armadillo

#### Tetrahedron stable

#### Deform large mesh