#pragma once

#include <CL/cl2.hpp>

void resolve_fission_part_collision(cl_double* pos, cl_double* vel, cl_double* masses, cl_double* radii,
                                    size_t num_parts, cl_double e, size_t i, size_t j, cl_double fusion_thresh,
                                    cl_double** endpos, cl_double** endvel,
                                    cl_double** endmasses, cl_double** endradii,
                                    size_t* endparts);

void fission_simulation_loop(cl_double* pos, cl_double* vel,
                             cl_double* masses, cl_double* radii,
                             cl_double* x_wall, cl_double* y_wall, cl_double* z_wall,
                             size_t num_parts, cl_double e, cl_double max_time, cl_double fusion_thresh);