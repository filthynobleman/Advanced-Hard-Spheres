# Advanced Hard Spheres
This is a tool which implements a generalization of the classic hard spheres
simulation method, taking advantage of OpenCL to parallelize computation.  
At this moment, the project is still in its development version, but all the
main features work correctly. See below for the list of future changes.

## Prerequisites
You will need a working C++ compiler, an IDE capable of importing Visual Studio
solution and an implementation of the OpenCL library.

## Install
Simply clone the repository and open the solution with an IDE. Be sure that OpenCL
is correctly configured and compile the solution.

## How to Use It
The correct syntax to run the tool is
```
AHSSimulation.exe INPUT_FILE [OUTPUT_FILE]
```
where `INPUT_FILE` is the path to the file containing the definition of the model
to simulate and `OUTPUT_FILE` is the path to the file where the simulation results
are saved. If the latter is not given, it will be automatically generated from the
informations about the input model.

### Input File Syntax
The input file must have the following format:
```
MODEL_NAME=<string>
NUM_PARTS=<positive integer>
STOP_TIME=<positive real>
ELASTIC_COEFF=<real between 0 and 1>
X_WALL=<real>, <real>
Y_WALL=<real>, <real>
Z_WALL=<real>, <real>
SIM_TYPE=<INELASTIC|FUSION|FISSION>
[THRESHOLD=<non-negative real>] // Only if SIM_TYPE == FUSION or SIM_TYPE == FISSION
POSITIONS=<RANDOM|GIVEN>
[<real>, <real>, <real>] // Only if POSITIONS == GIVEN. Must be repeated for NUM_PARTS rows
VELOCITIES=<RANDOM|GIVEN>
[<real>, <real>, <real>] // Only if VELOCITIES == GIVEN. Must be repeated for NUM_PARTS rows
MASSES=<RANDOM|GIVEN>
[<real>] // Only if MASSES == GIVEN. Must be repeated for NUM_PARTS rows
RADII=<RANDOM|GIVEN>
[<real>] // Only if RADII == GIVEN. Must be repeated for NUM_PARTS rows
```
The explaination of the parameters is the following:
  * `MODEL_NAME`: A string identifying the name of the model. If the output file is not given,
                  the result will be stored inside `MODEL_NAME.out`
  * `NUM_PARTS`: A strictly positive integer representing the number of particles to simulate.
  * `STOP_TIME`: A strictly positive real value representing the number of time units to simulate.
  * `ELASTIC_COEFF`: A real value between 0 and 1 representing the elastic coefficient of the
                     collisions. The elastic coefficient is equals to 1 minus the loss of kinetic
                     energy at each collision.
  * `X_WALL`, `Y_WALL`, `Z_WALL`: Three couples of real values representing the position of the walls
                                  along the *X*, *Y* and *Z* axes.
  * `SIM_TYPE`: An enumerated value representing the type of the simulation. Each simulation type will
                be explained below.
  * `THRESHOLD`: A non-negative real value. It defines the velocity threshold over which the fission,
                 or the fusion, are enabled.
  * `POSITION`: An enumerated value, possibly followed by a list of triplets of real values, representing
                the initial positions of the spheres. If this parameter is equals to `RANDOM`, then all
                the spheres are initialized in random positions. If this parameter is equals to `GIVEN`,
                then it must be followed by a list of triplet of real values, a triplet for each sphere,
                representing the initial positiions of the spheres.
  * `VELOCITIES`: Same as `POSITION`, but identifies the initial velocities of the spheres.
  * `MASSES`: Same as `POSITION`, but if `GIVEN` it is followed by a real value for each sphere, rather
              than by triplets, and represents the masses of the spheres.
  * `RADII`: Same as `MASSES`, but it represents the radii of the spheres.

### Output File Format
The output file is always binary. The *inelastic* model has its own output format. The *fission* and
*fusion* models share the same output format, different from the format used by the *inelastic* model.

#### The Inelastic Model
The output file begins with an header containing:
  * 64 bits (8 bytes): unsigned integer representing the simulation type (0 is inelastic, 1 is fusion, 2 is fission).
  * 64 bits (8 bytes): unsigned integer representing the number of particles.
  * 64 bits (8 bytes): double precision floating point value representing the elastic coefficient.
  * 64 bits (8 bytes): double precision floating point value representing the simulation time units.
  * 128 bits (16 bytes): two double precision floating point values representing the walls coordinates along the *X* axis.
  * 128 bits (16 bytes): two double precision floating point values representing the walls coordinates along the *Y* axis.
  * 128 bits (16 bytes): two double precision floating point values representing the walls coordinates along the *Z* axis.
  * 64 bits (8 bytes) for each particle: a double precision floating point value representing the radius of a particle.

The file is followed by a line for each simulation step. Each line contains a double precision floating point value
representing the first point in time (after the previous line) where a collision occurred. It is then followed by
three double precision values for each particle representing its position at the time of the collision and by
three double precision values for each particle representing its velocity immediately after the collision.

#### The Fusion/Fission Model
TODO


## Types of Model
Here follows the three possible types of model.

### The Inelastic Model
This model simulates a set of particles, each of its mass and radius, moving in a closed volume at constant velocity.
At a certain point in time, two particles will collide. A collision is resolved as a classical collision between two
rigid bodies. It is possible to define the loss rate of the kinetic energy at each collision. In this way, it is
possible to simulate perfectly elastic, perfectly inelastic as well as partially elastic models.

### The Fusion Model
TODO

### The Fission Model
TODO

## Future Changes
Here a list of the possible future changes. As I will think to other changes, I will also add them here.
  * The code needs to be reorganized and cleaned.
  * Each possible improvement in performance will be implemented.
  * Seriously considering the elimination of the fission model.
  * Considering the elimination of the fusion model.
  * Considering changes in the input format.
  * Considering changes in the output format.
  * A graphic renderer in OpenGL will be implemented to visualize the simulation results.
  * Switch to Make.
  * A complete and formal documentation of the models, the code and the input and output format will be provided.
