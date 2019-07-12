#include "inelastic.h"
#include "fusion.h"
#include "fission.h"
#include <math.h>

#include <sstream>

#include <iostream>

cl_double dot_prod(cl_double* v1, cl_double* v2, size_t len)
{
    cl_double result = 0;
    for (register size_t i = 0; i < len; i++)
    {
        result += v1[i] * v2[i];
    }
    return result;
}

void cross_prod(cl_double* v1, cl_double* v2, cl_double* u)
{
    u[0] = v1[1] * v2[2] - v1[2] * v2[1];
    u[1] = v1[2] * v2[0] - v1[0] * v2[2];
    u[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

#define MAX(x, y)   ((x) < (y) ? (y) : (x))
#define MIN(x, y)   ((x) < (y) ? (x) : (y))
#define ABS(x)      ((x) < 0 ? -(x) : (x))

void resolve_inelastic_part_collision(cl_double* pos, cl_double* vel, cl_double* masses,
                                      size_t num_parts, cl_double e,
                                      size_t i, size_t j)
{
    // Positions and velocities
    cl_double* pi = pos + (3 * i);
    cl_double* vi = vel + (3 * i);
    cl_double* pj = pos + (3 * j);
    cl_double* vj = vel + (3 * j);
    // Pij and Vij
    cl_double pij[3] = { pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2] };
    cl_double vij[3] = { vi[0] - vj[0], vi[1] - vj[1], vi[2] - vj[2] };

    // Inelastic component
    cl_double inelastic[3];
    for (size_t k = 0; k < 3; k++)
        inelastic[k] = (masses[i] * vi[k] + masses[j] * vj[k]) / (masses[i] + masses[j]);

    // Elastic component
    cl_double pvij = dot_prod(pij, vij, 3);
    cl_double pij_norm2 = dot_prod(pij, pij, 3);
    cl_double elastic[3];
    for (size_t k = 0; k < 3; k++)
        elastic[k] = e * (2 * pvij * pij[k] / pij_norm2 - vij[k]) / (masses[i] + masses[j]);

    // Update the velocities
    for (size_t k = 0; k < 3; k++)
    {
        vi[k] = inelastic[k] - masses[i] * elastic[k];
        vj[k] = inelastic[k] + masses[j] * elastic[k];
    }
}



void resolve_fusion_part_collision(cl_double* pos, cl_double* vel, cl_double* masses, cl_double* radii,
                                   size_t num_parts, cl_double e, size_t i, size_t j, cl_double fusion_thresh, 
                                   cl_double** endpos, cl_double** endvel,
                                   cl_double** endmasses, cl_double** endradii,
                                   size_t* endparts)
{
    // Positions and velocities
    cl_double* pi = pos + (3 * i);
    cl_double* vi = vel + (3 * i);
    cl_double* pj = pos + (3 * j);
    cl_double* vj = vel + (3 * j);
    // Pij and Vij
    cl_double pij[3] = { pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2] };
    cl_double vij[3] = { vi[0] - vj[0], vi[1] - vj[1], vi[2] - vj[2] };

    // Inelastic component
    cl_double inelastic[3];
    for (size_t k = 0; k < 3; k++)
        inelastic[k] = (masses[i] * vi[k] + masses[j] * vj[k]) / (masses[i] + masses[j]);

    // Elastic component
    cl_double pvij = dot_prod(pij, vij, 3);
    cl_double pij_norm2 = dot_prod(pij, pij, 3);
    cl_double elastic[3];
    for (size_t k = 0; k < 3; k++)
        elastic[k] = e * (2 * pvij * pij[k] / pij_norm2 - vij[k]) / (masses[i] + masses[j]);

    // If the velocity component along the collision axis exceeds the threshold, fuse the particles
    if (ABS(pvij / sqrt(pij_norm2)) > fusion_thresh)
    {
        //std::cout << "Fusion collision between particles " << i << " and " << j << std::endl;
        // Allocate new space for new data
        *endpos = (cl_double*)calloc((num_parts - 1) * 3, sizeof(cl_double));
        *endvel = (cl_double*)calloc((num_parts - 1) * 3, sizeof(cl_double));
        *endmasses = (cl_double*)calloc(num_parts - 1, sizeof(cl_double));
        *endradii = (cl_double*)calloc(num_parts - 1, sizeof(cl_double));
        *endparts = num_parts - 1;
        if (endpos == NULL || endvel == NULL || endmasses == NULL || endradii == NULL)
        {
            std::stringstream ss;
            ss << "Error occurred while allocating memory in particle velocity update for fusion model." << std::endl;
            throw std::runtime_error(ss.str());
        }

        // Copy the content, excluding the fused element
        std::memcpy((*endpos), pos, 3 * MIN(i, j) * sizeof(cl_double));
        std::memcpy((*endpos) + 3 * MIN(i, j), pos + 3 * (MIN(i, j) + 1), (MAX(i, j) - MIN(i, j) - 1) * 3 * sizeof(cl_double));
        std::memcpy((*endpos) + 3 * (MAX(i, j) - 1), pos + 3 * (MAX(i, j) + 1), (num_parts - MAX(i, j) - 1) * 3 * sizeof(cl_double));

        std::memcpy((*endvel), vel, 3 * MIN(i, j) * sizeof(cl_double));
        std::memcpy((*endvel) + 3 * MIN(i, j), vel + 3 * (MIN(i, j) + 1), (MAX(i, j) - MIN(i, j) - 1) * 3 * sizeof(cl_double));
        std::memcpy((*endvel) + 3 * (MAX(i, j) - 1), vel + 3 * (MAX(i, j) + 1), (num_parts - MAX(i, j) - 1) * 3 * sizeof(cl_double));

        std::memcpy((*endmasses), masses, MIN(i, j) * sizeof(cl_double));
        std::memcpy((*endmasses) + MIN(i, j), masses + (MIN(i, j) + 1), (MAX(i, j) - MIN(i, j) - 1) * sizeof(cl_double));
        std::memcpy((*endmasses) + (MAX(i, j) - 1), masses + (MAX(i, j) + 1), (num_parts - MAX(i, j) - 1) * sizeof(cl_double));

        std::memcpy((*endradii), radii, MIN(i, j) * sizeof(cl_double));
        std::memcpy((*endradii) + MIN(i, j), radii + (MIN(i, j) + 1), (MAX(i, j) - MIN(i, j) - 1) * sizeof(cl_double));
        std::memcpy((*endradii) + (MAX(i, j) - 1), radii + (MAX(i, j) + 1), (num_parts - MAX(i, j) - 1) * sizeof(cl_double));

        // Finally, assign position, velocity, mass and radius to the new particle
        for (size_t k = 0; k < 3; k++)
        {
            (*endpos)[3 * ((*endparts) - 1) + k] = (pi[k] + pj[k]) / 2;
            (*endvel)[3 * ((*endparts) - 1) + k] = inelastic[k];
        }
        (*endmasses)[(*endparts) - 1] = masses[i] + masses[j];
        (*endradii)[(*endparts) - 1] = pow(pow(radii[i], 3) + pow(radii[j], 3), 1.0 / 3.0);
    }
    // Otherwise, the update is the normal inelastic collision update
    else
    {
        // Allocate new space for new data
        *endpos = (cl_double*)calloc(num_parts * 3, sizeof(cl_double));
        *endvel = (cl_double*)calloc(num_parts * 3, sizeof(cl_double));
        *endmasses = (cl_double*)calloc(num_parts, sizeof(cl_double));
        *endradii = (cl_double*)calloc(num_parts, sizeof(cl_double));
        *endparts = num_parts;
        if (endpos == NULL || endvel == NULL || endmasses == NULL || endradii == NULL)
        {
            std::stringstream ss;
            ss << "Error occurred while allocating memory in particle velocity update for fusion model." << std::endl;
            throw std::runtime_error(ss.str());
        }

        // Copy the data
        std::memcpy(*endpos, pos, 3 * num_parts * sizeof(cl_double));
        std::memcpy(*endvel, vel, 3 * num_parts * sizeof(cl_double));
        std::memcpy(*endmasses, masses, num_parts * sizeof(cl_double));
        std::memcpy(*endradii, radii, num_parts * sizeof(cl_double));

        // Update the velocities
        for (size_t k = 0; k < 3; k++)
        {
            (*endvel)[3 * i + k] = inelastic[k] - masses[i] * elastic[k];
            (*endvel)[3 * j + k] = inelastic[k] + masses[j] * elastic[k];
        }
    }
}


void resolve_fission_part_collision(cl_double* pos, cl_double* vel, 
                                    cl_double* masses, cl_double* radii, 
                                    size_t num_parts, cl_double e, size_t i, size_t j, cl_double fusion_thresh, 
                                    cl_double** endpos, cl_double** endvel, 
                                    cl_double** endmasses, cl_double** endradii, 
                                    size_t* endparts)
{
    // Positions and velocities
    cl_double* pi = pos + (3 * i);
    cl_double* vi = vel + (3 * i);
    cl_double* pj = pos + (3 * j);
    cl_double* vj = vel + (3 * j);
    // Pij and Vij
    cl_double pij[3] = { pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2] };
    cl_double vij[3] = { vi[0] - vj[0], vi[1] - vj[1], vi[2] - vj[2] };

    // Inelastic component
    cl_double inelastic[3];
    for (size_t k = 0; k < 3; k++)
        inelastic[k] = (masses[i] * vi[k] + masses[j] * vj[k]) / (masses[i] + masses[j]);

    // Elastic component
    cl_double pvij = dot_prod(pij, vij, 3);
    cl_double pij_norm2 = dot_prod(pij, pij, 3);
    cl_double elastic[3];
    for (size_t k = 0; k < 3; k++)
        elastic[k] = e * (2 * pvij * pij[k] / pij_norm2 - vij[k]) / (masses[i] + masses[j]);

    // If the velocity component along the collision axis exceeds the threshold, break one of the particles
    if (ABS(pvij / sqrt(pij_norm2)) > fusion_thresh)
    {
        //std::cout << "Fission collision between particles " << i << " and " << j << std::endl;
        // Allocate new space for new data
        *endpos = (cl_double*)calloc((num_parts + 1) * 3, sizeof(cl_double));
        *endvel = (cl_double*)calloc((num_parts + 1) * 3, sizeof(cl_double));
        *endmasses = (cl_double*)calloc(num_parts + 1, sizeof(cl_double));
        *endradii = (cl_double*)calloc(num_parts + 1, sizeof(cl_double));
        *endparts = num_parts + 1;
        if (*endpos == NULL || *endvel == NULL || *endmasses == NULL || *endradii == NULL)
        {
            std::stringstream ss;
            ss << "Error occurred while allocating memory in particle velocity update for fission model." << std::endl;
            throw std::runtime_error(ss.str());
        }

        // Copy the content to the new data
        std::memcpy(*endpos, pos, 3 * num_parts * sizeof(cl_double));
        std::memcpy(*endvel, vel, 3 * num_parts * sizeof(cl_double));
        std::memcpy(*endmasses, masses, num_parts * sizeof(cl_double));
        std::memcpy(*endradii, radii, num_parts * sizeof(cl_double));

        // Determine the broken particle. Select one particle at random, first
        size_t brok = rand() % 2;
        brok = brok * i + (1 - brok) * j;   // i.e. rand = 1 -> brok = i, rand = 0 -> brok = j
        // Select the fastest particle, along the collision axis
        if (ABS(dot_prod(pij, vi, 3)) < ABS(dot_prod(pij, vj, 3)))
            brok = i;
        else if (ABS(dot_prod(pij, vi, 3)) > ABS(dot_prod(pij, vj, 3)))
            brok = j;
        // If velocities are the same, select by lower density
        else if (pow(radii[i], 3) / masses[i] < pow(radii[j], 3) / masses[j])
            brok = i;
        else if (pow(radii[i], 3) / masses[i] > pow(radii[j], 3) / masses[j])
            brok = j;
        // If densities are the same, select higher dimensions
        else if (radii[i] > radii[j])
            brok = i;
        else if (radii[i] < radii[j])
            brok = j;
        // If all the quantities compared are equals, then the random selection is kept
        // Select the other particle
        size_t other = brok == i ? j : i;

        // Update the particles' velocity using normal collision update
        cl_double sign = brok == i ? 1 : -1;
        for (size_t k = 0; k < 3; k++)
        {
            (*endvel)[3 * brok + k] = inelastic[k] - sign * e * masses[other] * elastic[k];
            (*endvel)[3 * other + k] = inelastic[k] + sign * e * masses[brok] * elastic[k];
        }

        // Define the radii of the new particles
        (*endradii)[brok] = pow(4, 1.0 / 3.0) / 2 * radii[brok];
        (*endradii)[*endparts - 1] = pow(4, 1.0 / 3.0) / 2 * radii[brok];

        // Define the masses of the new particles
        (*endmasses)[brok] = masses[brok] / 2;
        (*endmasses)[*endparts - 1] = masses[brok] / 2;

        // Compute the normal component to the collision plane
        cl_double normal[3];
        cross_prod(vel + 3 * other, pij, normal);
        // Scale to unit. If length is zero, choose random
        cl_double normal_length = sqrt(dot_prod(normal, normal, 3));
        if (normal_length == 0)
        {
            for (size_t k = 0; k < 3; k++)
                normal[k] = ((cl_double)rand()) / RAND_MAX;
            normal_length = sqrt(dot_prod(normal, normal, 3));
        }
        for (size_t k = 0; k < 3; k++)
            normal[k] /= normal_length;

        // Compute the length of the normal component
        cl_double v_prev = dot_prod(vel + 3 * brok, vel + 3 * brok, 3);
        cl_double v_next = dot_prod(*endvel + 3 * brok, *endvel + 3 * brok, 3);
        normal_length = sqrt(ABS(v_prev - v_next));

        // Finally, assign position and velocity to the new particles
        for (size_t k = 0; k < 3; k++)
        {
            (*endpos)[3 * ((*endparts) - 1) + k] =  pos[3 * brok + k] + normal[k] * (*endradii)[brok];
            (*endpos)[3 * brok + k] =               pos[3 * brok + k] - normal[k] * (*endradii)[brok];
            (*endvel)[3 * ((*endparts) - 1) + k] =  (*endvel)[3 * brok + k] + normal[k] * normal_length;
            (*endvel)[3 * brok + k] =               (*endvel)[3 * brok + k] - normal[k] * normal_length;
        }
    }
    // Otherwise, the update is the normal inelastic collision update
    else
    {
        // Allocate new space for new data
        *endpos = (cl_double*)calloc(num_parts * 3, sizeof(cl_double));
        *endvel = (cl_double*)calloc(num_parts * 3, sizeof(cl_double));
        *endmasses = (cl_double*)calloc(num_parts, sizeof(cl_double));
        *endradii = (cl_double*)calloc(num_parts, sizeof(cl_double));
        *endparts = num_parts;
        if (*endpos == NULL || *endvel == NULL || *endmasses == NULL || *endradii == NULL)
        {
            std::stringstream ss;
            ss << "Error occurred while allocating memory in particle velocity update for fission model." << std::endl;
            throw std::runtime_error(ss.str());
        }

        // Copy the data
        std::memcpy(*endpos, pos, 3 * num_parts * sizeof(cl_double));
        std::memcpy(*endvel, vel, 3 * num_parts * sizeof(cl_double));
        std::memcpy(*endmasses, masses, num_parts * sizeof(cl_double));
        std::memcpy(*endradii, radii, num_parts * sizeof(cl_double));

        // Update the velocities
        for (size_t k = 0; k < 3; k++)
        {
            (*endvel)[3 * i + k] = inelastic[k] - masses[i] * elastic[k];
            (*endvel)[3 * j + k] = inelastic[k] + masses[j] * elastic[k];
        }
    }
}