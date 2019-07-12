#pragma once

#include <CL/cl2.hpp>

void resolve_inelastic_part_collision(cl_double* pos, cl_double* vel, cl_double* masses,
                                      size_t num_parts, cl_double e,
                                      size_t i, size_t j);

void inelastic_simulation_loop(cl_double* pos, cl_double* vel,
                               cl_double* masses, cl_double* radii,
                               cl_double* x_wall, cl_double* y_wall, cl_double* z_wall,
                               size_t num_parts, cl_double e, cl_double max_time);