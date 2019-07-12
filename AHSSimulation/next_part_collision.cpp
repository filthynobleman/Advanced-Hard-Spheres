#include "shared.h"
#include "CLSettings.h"
#include <sstream>
#include <math.h>

#include <iostream>

void next_part_collision(cl_double* in_pos, cl_double* in_vel, cl_double* radii, size_t num_parts,
                         size_t* i, size_t* j, cl_double* delta_time)
{
    // Get the device and create the list
    cl::Device& device = CLSettings::get_device();
    cl::vector<cl::Device> devices;
    devices.push_back(device);

    // Create the context
    static bool first_run = false;
    cl_int status = CL_SUCCESS;
    static cl::Context context(devices, NULL, NULL, NULL, &status);
    if (!first_run && status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "Errors occurred while creating the OpenCL context for position update." << std::endl;
        throw new std::runtime_error(ss.str());
    }

    // Create and build the program
    cl::vector<std::string> sources;
    sources.push_back(CLSettings::get_source_part_collision());
    static cl::Program program(context, sources, &status);
    if (!first_run && status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "Errors occurred while creating the OpenCL program for position update." << std::endl;
        throw new std::runtime_error(ss.str());
    }
    if (!first_run)
        status = program.build({ device });
    if (!first_run && status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "Errors occurred while building the OpenCL program for position update." << std::endl;
        std::string build_log;
        status = program.getBuildInfo(device, CL_PROGRAM_BUILD_LOG, &build_log);
        if (status != CL_SUCCESS)
        {
            ss << "Errors occurred while retrieving the build log." << std::endl;
            throw new std::runtime_error(ss.str());
        }
        ss << "********** BUILD LOG BEGIN **********" << std::endl
            << build_log
            << "**********  BUILD LOG END  **********" << std::endl;
        throw std::runtime_error(ss.str());
    }

    // Create the buffers
    // Input positions
    cl::Buffer cl_in_pos(context, CL_MEM_READ_ONLY, 3 * num_parts * sizeof(cl_double), &status);
    if (status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "Errors occurred while creating an OpenCL input buffer for position update." << std::endl;
        throw new std::runtime_error(ss.str());
    }
    // Input velocities
    cl::Buffer cl_in_vel(context, CL_MEM_READ_ONLY, 3 * num_parts * sizeof(cl_double), &status);
    if (status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "Errors occurred while creating an OpenCL input buffer for position update." << std::endl;
        throw new std::runtime_error(ss.str());
    }
    // Input radii
    cl::Buffer cl_in_radii(context, CL_MEM_READ_ONLY, num_parts * sizeof(cl_double), &status);
    if (status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "Errors occurred while creating an OpenCL input buffer for position update." << std::endl;
        throw new std::runtime_error(ss.str());
    }
    // Output delta time
    cl::Buffer cl_out_delta_time(context, CL_MEM_WRITE_ONLY, num_parts * num_parts * sizeof(cl_double), &status);
    if (status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "Errors occurred while creating an OpenCL output buffer for position update." << std::endl;
        throw new std::runtime_error(ss.str());
    }

    // Create the command queue and enqueue the buffers
    static cl::CommandQueue queue(context, CL_QUEUE_PROFILING_ENABLE, &status);
    if (!first_run && status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "Errors occurred while creating the OpenCL command queue in positions update." << std::endl;
        ss << "Error code: " << status << std::endl;
        throw std::runtime_error(ss.str());
    }
    status = queue.enqueueWriteBuffer(cl_in_pos, CL_TRUE, 0, 3 * num_parts * sizeof(cl_double), in_pos);
    status |= queue.enqueueWriteBuffer(cl_in_vel, CL_TRUE, 0, 3 * num_parts * sizeof(cl_double), in_vel);
    status |= queue.enqueueWriteBuffer(cl_in_radii, CL_TRUE, 0, num_parts * sizeof(cl_double), radii);
    if (status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "Errors occurred while writing buffers on OpenCL device memory for positions update." << std::endl;
        throw new std::runtime_error(ss.str());
    }

    // Create the kernel and set the arguments
    static cl::Kernel kernel(program, PART_COLLISION_KERNEL_NAME, &status);
    if (!first_run && status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "Errors occurred while creating the OpenCL kernel for position update." << std::endl;
        throw new std::runtime_error(ss.str());
    }

    status = kernel.setArg(0, cl_in_pos);
    status |= kernel.setArg(1, cl_in_vel);
    status |= kernel.setArg(2, cl_in_radii);
    status |= kernel.setArg(3, num_parts);
    status |= kernel.setArg(4, cl_out_delta_time);
    if (status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "errors occurred while setting the arguments of the opencl kernel for position update." << std::endl;
        throw new std::runtime_error(ss.str());
    }

    // Execute the kernel
    queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(num_parts, num_parts), cl::NullRange);

    // Retrieve the results
    cl_double* delta_times = (cl_double*)calloc(num_parts * num_parts, sizeof(cl_double));
    if (delta_times == NULL)
    {
        std::stringstream ss;
        ss << "Errors occurred during the allocation of memory." << std::endl;
        throw std::runtime_error(ss.str());
    }
    queue.enqueueReadBuffer(cl_out_delta_time, CL_TRUE, 0, num_parts * num_parts * sizeof(cl_double), delta_times);
    queue.finish();

    // Compute the minimum delta_time
    *delta_time = INFINITY;
    for (register size_t ii = 0; ii < num_parts; ii++)
    {
        for (register size_t jj = ii + 1; jj < num_parts; jj++)
        {
            size_t k = ii * num_parts + jj;
            if (delta_times[k] < *delta_time)
            {
                *delta_time = delta_times[k];
                *i = ii;
                *j = jj;
            }
        }
    }

    free(delta_times);

    first_run = true;
}