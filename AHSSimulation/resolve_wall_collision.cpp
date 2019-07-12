#include "shared.h"


void resolve_wall_collision(cl_double* in_pos, cl_double* in_vel,
                            size_t p, cl_double* collision_axis)
{
    // Get the position and velocity of particle p
    cl_double x = in_pos[p * 3];
    cl_double y = in_pos[p * 3 + 1];
    cl_double z = in_pos[p * 3 + 2];

    cl_double vx = in_vel[p * 3];
    cl_double vy = in_vel[p * 3 + 1];
    cl_double vz = in_vel[p * 3 + 2];

    // The collision against a wall inverts the sign of the component along the collision axis
    // Compute the dot product between velocity and collision axis
    cl_double dot_cv = vx * collision_axis[0] +
                       vy * collision_axis[1] +
                       vz * collision_axis[2];
    // Scale the collision axis
    collision_axis[0] *= dot_cv;
    collision_axis[1] *= dot_cv;
    collision_axis[2] *= dot_cv;
    // Change the particle's velocity
    vx -= 2 * collision_axis[0];
    vy -= 2 * collision_axis[1];
    vz -= 2 * collision_axis[2];

    // Set the velocity inside the array
    in_vel[p * 3] = vx;
    in_vel[p * 3 + 1] = vy;
    in_vel[p * 3 + 2] = vz;
}