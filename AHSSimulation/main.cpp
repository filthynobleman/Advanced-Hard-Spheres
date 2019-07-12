#include <iostream>
#include <string>
#include <stdio.h>
#include <chrono>

#include "ahs.h"

#define NAME_MAX_LEN    256

int main(int argc, char** argv)
{
    std::cout << "Executing " << argv[0] << "..." << std::endl;
    // Two arguments:
    // 1. A settings file
    // 2. An output file
    if (argc < 2)
    {
        std::cerr << "Cannot execute AHSSimulation with less than one arguments." << std::endl;
        return 1;
    }

    std::string inputfile(argv[1]);
    std::string outputfile;

    // Input file must exists
    FILE* instream; 
    fopen_s(&instream, inputfile.c_str(), "r");
    if (instream == NULL)
    {
        std::cerr << "Cannot open file " << inputfile << " for reading." << std::endl;
        return 1;
    }

    // Read the parameters
    std::cout << "Starting reading from file " << inputfile << std::endl;

    char model_name[NAME_MAX_LEN];
    int model_name_len;
    cl_double max_time;
    size_t num_parts;
    cl_double* positions;
    cl_double* velocities;
    cl_double* masses;
    cl_double* radii;
    cl_double x_wall[2];
    cl_double y_wall[2];
    cl_double z_wall[2];
    char simname[NAME_MAX_LEN];
    char vecgenmode[NAME_MAX_LEN];
    int simtype;
    cl_double e;
    cl_double threshold;

    int status;

    // Model name
    status = fscanf_s(instream, "MODEL_NAME=%s%n\n", model_name, NAME_MAX_LEN, &model_name_len);
    if (status == 0)
    {
        std::cerr << "Error while reading model name." << std::endl;
        return 1;
    }
    else
        model_name[model_name_len] = 0;

    // Number of particles
    fscanf_s(instream, "NUM_PARTS=%llu\n", &num_parts);
    if (status == 0)
    {
        std::cerr << "Error while reading the number of particles." << std::endl;
        return 1;
    }

    // Time horizon
    status = fscanf_s(instream, "STOP_TIME=%lf\n", &max_time);
    if (status == 0)
    {
        std::cerr << "Error while reading the time horizon." << std::endl;
        return 1;
    }
    if (max_time <= 0)
    {
        std::cerr << "The time horizon must be a strictly positive real value. Given value is " << max_time << std::endl;
        return 1;
    }

    // Elasticity coefficient
    status = fscanf_s(instream, "ELASTIC_COEFF=%lf\n", &e);
    if (status == 0)
    {
        std::cerr << "Error while reading the elasticity coefficient for the collisions." << std::endl;
        return 1;
    }
    else if (e < 0 || e > 1)
    {
        std::cerr << "Invalid value for the elasticity coefficient." << std::endl;
        std::cerr << "Legal values are in the interval [0, 1]. Given value is " << e << std::endl;
        return 1;
    }

    // X walls
    status = fscanf_s(instream, "X_WALL=%lf, %lf\n", x_wall, x_wall + 1);
    if (status == 0)
    {
        std::cerr << "Error while reading the values for walls along the X axis." << std::endl;
        return 1;
    }
    // Y walls
    status = fscanf_s(instream, "Y_WALL=%lf, %lf\n", y_wall, y_wall + 1);
    if (status == 0)
    {
        std::cerr << "Error while reading the values for walls along the Y axis." << std::endl;
        return 1;
    }
    // Z walls
    status = fscanf_s(instream, "Z_WALL=%lf, %lf\n", z_wall, z_wall + 1);
    if (status == 0)
    {
        std::cerr << "Error while reading the values for walls along the Z axis." << std::endl;
        return 1;
    }

    // Simulation type
    status = fscanf_s(instream, "SIM_TYPE=%s\n", simname, NAME_MAX_LEN);
    if (status == 0)
    {
        std::cerr << "Error while reading the simulation type." << std::endl;
        return 1;
    }
    if (strcmp(simname, "INELASTIC") == 0)
        simtype = 0;
    else if (strcmp(simname, "FUSION") == 0)
        simtype = 1;
    else if (strcmp(simname, "FISSION") == 0)
        simtype = 2;
    else
    {
        std::cerr << "Invalid value for the simulation type." << std::endl;
        std::cerr << "Legal values are \"INELASTIC\", \"FUSION\" and \"FISSION\". Given value is " << simname << std::endl;
        return 1;
    }

    // Possibly, the threshold
    if (simtype > 0)
    {
        status = fscanf_s(instream, "THRESHOLD=%lf\n", &threshold);
        if (status == 0)
        {
            std::cerr << "Error while reading the threshold." << std::endl;
            return 1;
        }
        if (threshold < 0)
        {
            std::cerr << "Threshold value for enabling fusion or fission must be a non-negative real value." << std::endl;
            std::cerr << "Given value is " << threshold << std::endl;
            return 1;
        }
    }

    // Positions
    status = fscanf_s(instream, "POSITIONS=%s\n", vecgenmode, NAME_MAX_LEN);
    if (status == 0)
    {
        std::cerr << "Error while reading the position input mode." << std::endl;
        return 1;
    }
    positions = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    if (positions == NULL)
    {
        std::cerr << "Error while allocating space for position vectors." << std::endl;
        return 1;
    }
    if (strcmp(vecgenmode, "RANDOM") == 0)
    {
        for (size_t i = 0; i < 3 * num_parts; i++)
            positions[i] = ((cl_double)rand()) / RAND_MAX;
    }
    else if (strcmp(vecgenmode, "GIVEN") == 0)
    {
        for (size_t i = 0; i < num_parts; i++)
        {
            status = fscanf_s(instream, "%lf, %lf, %lf\n", positions + 3 * i, positions + 3 * i + 1, positions + 3 * i + 2);
            if (status == 0)
            {
                std::cerr << "Error while reading the " << i << "-th position vector." << std::endl;
                return 1;
            }
        }
    }
    else
    {
        std::cerr << "Invalid value for the positions input mode." << std::endl;
        std::cerr << "Legal values are \"RANDOM\" and \"GIVEN\". Given value is " << vecgenmode << std::endl;
        return 1;
    }
    
    // Velocities
    status = fscanf_s(instream, "VELOCITIES=%s\n", vecgenmode, NAME_MAX_LEN);
    if (status == 0)
    {
        std::cerr << "Error while reading the velocity input mode." << std::endl;
        return 1;
    }
    velocities = (cl_double*)calloc(3 * num_parts, sizeof(cl_double));
    if (velocities == NULL)
    {
        std::cerr << "Error while allocating space for velocity vectors." << std::endl;
        return 1;
    }
    if (strcmp(vecgenmode, "RANDOM") == 0)
    {
        for (size_t i = 0; i < 3 * num_parts; i++)
            velocities[i] = ((cl_double)rand()) / RAND_MAX;
    }
    else if (strcmp(vecgenmode, "GIVEN") == 0)
    {
        for (size_t i = 0; i < num_parts; i++)
        {
            status = fscanf_s(instream, "%lf, %lf, %lf\n", velocities + 3 * i, velocities + 3 * i + 1, velocities + 3 * i + 2);
            if (status == 0)
            {
                std::cerr << "Error while reading the " << i << "-th velocity vector." << std::endl;
                return 1;
            }
        }
    }
    else
    {
        std::cerr << "Invalid value for the velocity input mode." << std::endl;
        std::cerr << "Legal values are \"RANDOM\" and \"GIVEN\". Given value is " << vecgenmode << std::endl;
        return 1;
    }
    
    // Masses
    status = fscanf_s(instream, "MASSES=%s\n", vecgenmode, NAME_MAX_LEN);
    if (status == 0)
    {
        std::cerr << "Error while reading the masses input mode." << std::endl;
        return 1;
    }
    masses = (cl_double*)calloc(num_parts, sizeof(cl_double));
    if (masses == NULL)
    {
        std::cerr << "Error while allocating space for mass values." << std::endl;
        return 1;
    }
    if (strcmp(vecgenmode, "RANDOM") == 0)
    {
        for (size_t i = 0; i < num_parts; i++)
            masses[i] = (((cl_double)rand()) / RAND_MAX) * 0.9 + 0.1;
    }
    else if (strcmp(vecgenmode, "GIVEN") == 0)
    {
        for (size_t i = 0; i < num_parts; i++)
        {
            status = fscanf_s(instream, "%lf\n", masses + i);
            if (status == 0)
            {
                std::cerr << "Error while reading the " << i << "-th mass value." << std::endl;
                return 1;
            }
            if (masses[i] <= 0)
            {
                std::cerr << "Masses must be strictly positive real values." << std::endl;
                std::cerr << "Given value for the " << i << "-th mass is " << masses[i] << std::endl;
                return 1;
            }
        }
    }
    else
    {
        std::cerr << "Invalid value for the mass input mode." << std::endl;
        std::cerr << "Legal values are \"RANDOM\" and \"GIVEN\". Given value is " << vecgenmode << std::endl;
        return 1;
    }
    
    // Radii
    status = fscanf_s(instream, "RADII=%s\n", vecgenmode, NAME_MAX_LEN);
    if (status == 0)
    {
        std::cerr << "Error while reading the radii input mode." << std::endl;
        return 1;
    }
    radii = (cl_double*)calloc(num_parts, sizeof(cl_double));
    if (radii == NULL)
    {
        std::cerr << "Error while allocating space for radius values." << std::endl;
        return 1;
    }
    if (strcmp(vecgenmode, "RANDOM") == 0)
    {
        for (size_t i = 0; i < num_parts; i++)
            radii[i] = (((cl_double)rand()) / RAND_MAX) * 0.9 + 0.1;
    }
    else if (strcmp(vecgenmode, "GIVEN") == 0)
    {
        for (size_t i = 0; i < num_parts; i++)
        {
            status = fscanf_s(instream, "%lf\n", radii + i);
            if (status == 0)
            {
                std::cerr << "Error while reading the " << i << "-th radius value." << std::endl;
                return 1;
            }
            if (radii[i] <= 0)
            {
                std::cerr << "Radii must be strictly positive real values." << std::endl;
                std::cerr << "Given value for the " << i << "-th radius is " << radii[i] << std::endl;
                return 1;
            }
        }
    }
    else
    {
        std::cerr << "Invalid value for the radius input mode." << std::endl;
        std::cerr << "Legal values are \"RANDOM\" and \"GIVEN\". Given value is " << vecgenmode << std::endl;
        return 1;
    }

    // Close the input file
    fclose(instream);

    std::cout << "Successfully readed the input file." << std::endl;


    // Get the output filename
    if (argc > 2)
        outputfile = std::string(argv[2]);
    else
        outputfile = std::string(model_name) + ".out";
    std::cout << "Results will be saved to " << outputfile << std::endl;

    CLSettings::set_output_file(outputfile);
    cl::Device device(CLSettings::select_device());
    std::string devname;
    device.getInfo(CL_DEVICE_NAME, &devname);
    std::cout << "Selected device " << devname << std::endl << std::endl;
    CLSettings::set_device(device);

    std::cout << "Starting the simulation..." << std::endl;
    std::chrono::nanoseconds start_time;
    start_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch());
    try
    {
        if (simtype == 0)
            inelastic_simulation_loop(positions, velocities, masses, radii,
                x_wall, y_wall, z_wall,
                num_parts, e, max_time);
        else if (simtype == 1)
            fusion_simulation_loop(positions, velocities, masses, radii,
                x_wall, y_wall, z_wall,
                num_parts, e, max_time, threshold);
        else if (simtype == 2)
            fission_simulation_loop(positions, velocities, masses, radii,
                x_wall, y_wall, z_wall,
                num_parts, e, max_time, threshold);
    }
    catch (std::exception e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }
    catch (std::exception* e)
    {
        std::cerr << e->what() << std::endl;
        return 1;
    }

    std::chrono::nanoseconds end_time;
    end_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch());

    size_t elaps = (end_time - start_time).count();
    std::cout << "Simulation took " << elaps / std::pow(10, 9) << " seconds." << std::endl;

    return 0;
}