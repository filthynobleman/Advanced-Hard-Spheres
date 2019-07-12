#pragma once

#include <CL/cl2.hpp>

void update_positions(cl_double* in_pos, cl_double* in_vel,
                      size_t num_parts, cl_double delta_time,
                      cl_double* out_pos);

void next_wall_collision(cl_double* in_pos, cl_double* in_vel, cl_double* radii, size_t num_parts,
                         cl_double* x_wall, cl_double* y_wall, cl_double* z_wall,
                         size_t* p, cl_double* delta_time, cl_double* collision_axis);

void next_part_collision(cl_double* in_pos, cl_double* in_vel, cl_double* radii, size_t num_parts,
                         size_t* i, size_t* j, cl_double* delta_time);

void resolve_wall_collision(cl_double* in_pos, cl_double* in_vel, 
                            size_t p, cl_double* collision_axis);