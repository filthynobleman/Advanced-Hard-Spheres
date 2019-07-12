#include "inelastic.h"
#include "fusion.h"
#include "fission.h"
#include "shared.h"
#include "CLSettings.h"

#include <sstream>
#include <stdio.h>

#include <iostream>

#define MIN(x, y)       ((x) < (y) ? (x) : (y))
#define MAX(x, y)       ((x) > (y) ? (x) : (y))

void inelastic_simulation_loop(cl_double* pos, cl_double* vel,
                               cl_double* masses, cl_double* radii,
                               cl_double* x_wall, cl_double* y_wall, cl_double* z_wall,
                               size_t num_parts, cl_double e, cl_double max_time)
{
    // Open the file stream for the output
    FILE* stream;
    fopen_s(&stream, CLSettings::get_output_file().c_str(), "wb");
    if (stream == NULL)
    {
        std::stringstream ss;
        ss << "Some error occurred while opening the output file in the simulation loop." << std::endl;
        throw std::runtime_error(ss.str());
    }
    // Initialize the current time to zero
    cl_double time = 0;
    // Initialize the current and next positions and velocities
    cl_double* curpos = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* curvel = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* endpos = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* endvel = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    if (curpos == NULL || curvel == NULL || endpos == NULL || endvel == NULL)
    {
        std::stringstream ss;
        ss << "Some errors occurred while allocating memory in the simulation loop." << std::endl;
        throw std::runtime_error(ss.str());
    }
    // Copy the input values in the arrays
    std::memcpy(curpos, pos, 3 * num_parts * sizeof(cl_double));
    std::memcpy(curvel, vel, 3 * num_parts * sizeof(cl_double));
    std::memcpy(endpos, pos, 3 * num_parts * sizeof(cl_double));
    std::memcpy(endvel, vel, 3 * num_parts * sizeof(cl_double));
    // Output some informations about the system
    size_t simtype = SIMULATION_TYPE_INELSATIC;
    fwrite(&simtype, sizeof(size_t), 1, stream);            // Simulation type
    fwrite(&num_parts, sizeof(size_t), 1, stream);          // Number of particles
    fwrite(&e, sizeof(cl_double), 1, stream);               // Elasticity
    fwrite(&max_time, sizeof(cl_double), 1, stream);        // Time horizon
    fwrite(x_wall, sizeof(cl_double), 2, stream);           // X wall
    fwrite(y_wall, sizeof(cl_double), 2, stream);           // Y wall
    fwrite(z_wall, sizeof(cl_double), 2, stream);           // Z wall
    fwrite(radii, sizeof(cl_double), num_parts, stream);    // Particles' radii

    // Begin the simulation loop
    std::cout << "Simulation of a system of " << num_parts
              << " particles for " << max_time << " seconds." << std::endl;
    size_t time_idx = 0;
    while (time < max_time)
    {
        size_t p, i, j;
        cl_double delta_time, dt_wall, dt_part;
        cl_double coll_axis[3];

        // Check for the next collision
        next_wall_collision(curpos, curvel, radii, num_parts, x_wall, y_wall, z_wall, &p, &dt_wall, coll_axis);
        next_part_collision(curpos, curvel, radii, num_parts, &i, &j, &dt_part);
        delta_time = MIN(dt_wall, dt_part);

        // Update positions
        update_positions(curpos, curvel, num_parts, MAX(0, delta_time), endpos);
        
        // If a collision with a wall occurs first, resolve it
        if (dt_wall < dt_part)
            resolve_wall_collision(endpos, endvel, p, coll_axis);
        // Otherwise, resolve the collision between the particles
        else
            resolve_inelastic_part_collision(endpos, endvel, masses, num_parts, e, i, j);

        // If this step has seen an increment in time different from zero, then the system
        // has changed after a static period, so we can save the current status
        if (delta_time > 0)
        {
            //std::cout << "Saving output for time instant " << time << " (index = " << time_idx++ << ")" << std::endl;
            fwrite(&time, sizeof(cl_double), 1, stream);
            fwrite(curpos, sizeof(cl_double), 3 * num_parts, stream);
            fwrite(curvel, sizeof(cl_double), 3 * num_parts, stream);
            /*for (size_t p = 0; p < num_parts; p++)
            {
                std::cout << "X" << p << " = (" << curpos[p * 3] << ", "
                                                << curpos[p * 3 + 1] << ", "
                                                << curpos[p * 3 + 2] << ")" << std::endl;
            }*/
        }

        // Make the final state the current state and update the time
        std::memcpy(curpos, endpos, 3 * num_parts * sizeof(cl_double));
        std::memcpy(curvel, endvel, 3 * num_parts * sizeof(cl_double));
        time += MAX(0, delta_time);
    }

    // Close the stream
    fclose(stream);

    std::cout << "Simulation terminated." << std::endl;
}

void fusion_simulation_loop(cl_double* pos, cl_double* vel, cl_double* masses, cl_double* radii, 
                            cl_double* x_wall, cl_double* y_wall, cl_double* z_wall, 
                            size_t num_parts, cl_double e, cl_double max_time, cl_double fusion_thresh)
{
    // Open the file stream for the output
    FILE* stream;
    fopen_s(&stream, CLSettings::get_output_file().c_str(), "wb");
    if (stream == NULL)
    {
        std::stringstream ss;
        ss << "Some error occurred while opening the output file in the simulation loop." << std::endl;
        throw std::runtime_error(ss.str());
    }
    // Initialize the current time to zero
    cl_double time = 0;
    // Initialize the current and next positions and velocities
    cl_double* curpos = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* curvel = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* curmass = (cl_double*)calloc(num_parts, sizeof(cl_double));
    cl_double* curradii = (cl_double*)calloc(num_parts, sizeof(cl_double));
    cl_double* midpos = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* endpos = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* endvel = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* endmass = (cl_double*)calloc(num_parts, sizeof(cl_double));
    cl_double* endradii = (cl_double*)calloc(num_parts, sizeof(cl_double));
    if (curpos == NULL || curvel == NULL || curmass == NULL || curradii == NULL
        || midpos == NULL
        || endpos == NULL || endvel == NULL || endmass == NULL || endradii == NULL)
    {
        std::stringstream ss;
        ss << "Some errors occurred while allocating memory in the simulation loop." << std::endl;
        throw std::runtime_error(ss.str());
    }
    // Copy the input values in the arrays
    std::memcpy(curpos, pos, 3 * num_parts * sizeof(cl_double));
    std::memcpy(curvel, vel, 3 * num_parts * sizeof(cl_double));
    std::memcpy(curmass, masses, num_parts * sizeof(cl_double));
    std::memcpy(curradii, radii, num_parts * sizeof(cl_double));
    std::memcpy(midpos, pos, 3 * num_parts * sizeof(cl_double));
    std::memcpy(endpos, pos, 3 * num_parts * sizeof(cl_double));
    std::memcpy(endvel, vel, 3 * num_parts * sizeof(cl_double));
    std::memcpy(endmass, masses, num_parts * sizeof(cl_double));
    std::memcpy(endradii, radii, num_parts * sizeof(cl_double));
    // Output some informations about the system
    size_t simtype = SIMULATION_TYPE_FUSION;
    fwrite(&simtype, sizeof(size_t), 1, stream);            // Simulation type
    fwrite(&num_parts, sizeof(size_t), 1, stream);          // Number of particles
    fwrite(&e, sizeof(cl_double), 1, stream);               // Elasticity
    fwrite(&max_time, sizeof(cl_double), 1, stream);        // Time horizon
    fwrite(x_wall, sizeof(cl_double), 2, stream);           // X wall
    fwrite(y_wall, sizeof(cl_double), 2, stream);           // Y wall
    fwrite(z_wall, sizeof(cl_double), 2, stream);           // Z wall

    // Begin the simulation loop
    std::cout << "Simulation of a system of " << num_parts
        << " particles for " << max_time << " seconds." << std::endl;
    size_t time_idx = 0;
    size_t cur_num_parts = num_parts;
    size_t next_num_parts = num_parts;
    while (time < max_time)
    {
        midpos = (cl_double*)calloc(3 * cur_num_parts, sizeof(cl_double));
        if (midpos == NULL)
        {
            std::stringstream ss;
            ss << "Some errors occurred while allocating memory in the simulation loop." << std::endl;
            throw std::runtime_error(ss.str());
        }

        size_t p, i, j;
        cl_double delta_time, dt_wall, dt_part;
        cl_double coll_axis[3];

        // Check for the next collision
        next_wall_collision(curpos, curvel, curradii, cur_num_parts, x_wall, y_wall, z_wall, &p, &dt_wall, coll_axis);
        next_part_collision(curpos, curvel, curradii, cur_num_parts, &i, &j, &dt_part);
        delta_time = MIN(dt_wall, dt_part);

        // Update positions
        update_positions(curpos, curvel, cur_num_parts, MAX(0, delta_time), midpos);
        /*for (size_t kk = 0; kk < cur_num_parts; kk++)
        {
            std::cout << "  X" << kk << " = (" << curpos[3 * kk] << ", "
                                               << curpos[3 * kk + 1] << ", "
                                               << curpos[3 * kk + 2] << ")" << std::endl;
            std::cout << "  V" << kk << " = (" << curvel[3 * kk] << ", "
                                               << curvel[3 * kk + 1] << ", "
                                               << curvel[3 * kk + 2] << ")" << std::endl;
        }*/

        // If this step has seen an increment in time different from zero, then the system
        // has changed after a static period, so we can save the current status
        if (true)//(delta_time > 0)
        {
            //std::cout << "Saving output for time instant " << time << " (index = " << time_idx++ << ")" << std::endl;
            fwrite(&time, sizeof(cl_double), 1, stream);
            fwrite(&cur_num_parts, sizeof(size_t), 1, stream);
            fwrite(curradii, sizeof(cl_double), cur_num_parts, stream);
            fwrite(curpos, sizeof(cl_double), 3 * cur_num_parts, stream);
            fwrite(curvel, sizeof(cl_double), 3 * cur_num_parts, stream);
        }

        // If a collision with a wall occurs first, resolve it
        if (dt_wall < dt_part)
        {
            resolve_wall_collision(midpos, curvel, p, coll_axis);
            cl_double* tmp = endpos;
            endpos = midpos;
            midpos = endpos;
            endvel = curvel;
            endmass = curmass;
            endradii = curradii;
            next_num_parts = cur_num_parts;
        }
        // Otherwise, resolve the collision between the particles
        else
            resolve_fusion_part_collision(midpos, curvel, curmass, curradii, cur_num_parts,
                e, i, j, fusion_thresh,
                &endpos, &endvel, &endmass, &endradii, &next_num_parts);

        // Make the final state the current state and update the time
        cl_double* tmp;
        
        tmp = curpos;
        curpos = endpos;
        endpos = tmp;

        tmp = curvel;
        curvel = endvel;
        endvel = tmp;

        tmp = curmass;
        curmass = endmass;
        endmass = tmp;

        tmp = curradii;
        curradii = endradii;
        endradii = tmp;

        cur_num_parts = next_num_parts;

        time += MAX(0, delta_time);
    }

    // Close the stream
    fclose(stream);

    std::cout << "Simulation terminated." << std::endl;
}



void fission_simulation_loop(cl_double* pos, cl_double* vel, 
                             cl_double* masses, cl_double* radii, 
                             cl_double* x_wall, cl_double* y_wall, cl_double* z_wall, 
                             size_t num_parts, cl_double e, cl_double max_time, cl_double fusion_thresh)
{
    // Open the file stream for the output
    FILE* stream;
    fopen_s(&stream, CLSettings::get_output_file().c_str(), "wb");
    if (stream == NULL)
    {
        std::stringstream ss;
        ss << "Some error occurred while opening the output file in the simulation loop." << std::endl;
        throw std::runtime_error(ss.str());
    }
    // Initialize the current time to zero
    cl_double time = 0;
    // Initialize the current and next positions and velocities
    cl_double* curpos = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* curvel = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* curmass = (cl_double*)calloc(num_parts, sizeof(cl_double));
    cl_double* curradii = (cl_double*)calloc(num_parts, sizeof(cl_double));
    cl_double* midpos = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* endpos = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* endvel = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    cl_double* endmass = (cl_double*)calloc(num_parts, sizeof(cl_double));
    cl_double* endradii = (cl_double*)calloc(num_parts, sizeof(cl_double));
    if (curpos == NULL || curvel == NULL || curmass == NULL || curradii == NULL
        || midpos == NULL
        || endpos == NULL || endvel == NULL || endmass == NULL || endradii == NULL)
    {
        std::stringstream ss;
        ss << "Some errors occurred while allocating memory in the simulation loop." << std::endl;
        throw std::runtime_error(ss.str());
    }
    // Copy the input values in the arrays
    std::memcpy(curpos, pos, 3 * num_parts * sizeof(cl_double));
    std::memcpy(curvel, vel, 3 * num_parts * sizeof(cl_double));
    std::memcpy(curmass, masses, num_parts * sizeof(cl_double));
    std::memcpy(curradii, radii, num_parts * sizeof(cl_double));
    std::memcpy(midpos, pos, 3 * num_parts * sizeof(cl_double));
    std::memcpy(endpos, pos, 3 * num_parts * sizeof(cl_double));
    std::memcpy(endvel, vel, 3 * num_parts * sizeof(cl_double));
    std::memcpy(endmass, masses, num_parts * sizeof(cl_double));
    std::memcpy(endradii, radii, num_parts * sizeof(cl_double));
    // Output some informations about the system
    size_t simtype = SIMULATION_TYPE_FISSION;
    fwrite(&simtype, sizeof(size_t), 1, stream);            // Simulation type
    fwrite(&num_parts, sizeof(size_t), 1, stream);          // Number of particles
    fwrite(&e, sizeof(cl_double), 1, stream);               // Elasticity
    fwrite(&max_time, sizeof(cl_double), 1, stream);        // Time horizon
    fwrite(x_wall, sizeof(cl_double), 2, stream);           // X wall
    fwrite(y_wall, sizeof(cl_double), 2, stream);           // Y wall
    fwrite(z_wall, sizeof(cl_double), 2, stream);           // Z wall

    // Begin the simulation loop
    std::cout << "Simulation of a system of " << num_parts
        << " particles for " << max_time << " seconds." << std::endl;
    size_t time_idx = 0;
    size_t cur_num_parts = num_parts;
    size_t next_num_parts = num_parts;
    while (time < max_time)
    {
        midpos = (cl_double*)calloc(3 * cur_num_parts, sizeof(cl_double));
        if (midpos == NULL)
        {
            std::stringstream ss;
            ss << "Some errors occurred while allocating memory in the simulation loop." << std::endl;
            throw std::runtime_error(ss.str());
        }

        size_t p, i, j;
        cl_double delta_time, dt_wall, dt_part;
        cl_double coll_axis[3];

        // Check for the next collision
        next_wall_collision(curpos, curvel, curradii, cur_num_parts, x_wall, y_wall, z_wall, &p, &dt_wall, coll_axis);
        next_part_collision(curpos, curvel, curradii, cur_num_parts, &i, &j, &dt_part);
        delta_time = MIN(dt_wall, dt_part);

        // Update positions
        update_positions(curpos, curvel, cur_num_parts, MAX(0, delta_time), midpos);
        /*for (size_t kk = 0; kk < cur_num_parts; kk++)
        {
            std::cout << "  X" << kk << " = (" << curpos[3 * kk] << ", "
                                               << curpos[3 * kk + 1] << ", "
                                               << curpos[3 * kk + 2] << ")" << std::endl;
        }*/

        // If this step has seen an increment in time different from zero, then the system
        // has changed after a static period, so we can save the current status
        if (delta_time > 0)
        {
            //std::cout << "Saving output for time instant " << time << " (index = " << time_idx++ << ")" << std::endl;
            fwrite(&time, sizeof(cl_double), 1, stream);
            fwrite(&cur_num_parts, sizeof(size_t), 1, stream);
            fwrite(curradii, sizeof(cl_double), cur_num_parts, stream);
            fwrite(curpos, sizeof(cl_double), 3 * cur_num_parts, stream);
            fwrite(curvel, sizeof(cl_double), 3 * cur_num_parts, stream);
        }

        // If a collision with a wall occurs first, resolve it
        if (dt_wall < dt_part)
        {
            resolve_wall_collision(midpos, curvel, p, coll_axis);
            cl_double* tmp = endpos;
            endpos = midpos;
            midpos = endpos;
            endvel = curvel;
            endmass = curmass;
            endradii = curradii;
            next_num_parts = cur_num_parts;
        }
        // Otherwise, resolve the collision between the particles
        else
            resolve_fission_part_collision(midpos, curvel, curmass, curradii, cur_num_parts,
                e, i, j, fusion_thresh,
                &endpos, &endvel, &endmass, &endradii, &next_num_parts);

        // Make the final state the current state and update the time
        cl_double* tmp;

        tmp = curpos;
        curpos = endpos;
        endpos = tmp;

        tmp = curvel;
        curvel = endvel;
        endvel = tmp;

        tmp = curmass;
        curmass = endmass;
        endmass = tmp;

        tmp = curradii;
        curradii = endradii;
        endradii = tmp;

        cur_num_parts = next_num_parts;

        time += MAX(0, delta_time);
    }

    // Close the stream
    fclose(stream);

    std::cout << "Simulation terminated." << std::endl;
}