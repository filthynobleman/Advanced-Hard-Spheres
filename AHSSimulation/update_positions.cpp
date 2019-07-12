#include "shared.h"
#include "CLSettings.h"
#include <sstream>
#include <iostream>

void update_positions(cl_double* in_pos, cl_double* in_vel, size_t num_parts, 
                      cl_double delta_time, cl_double* out_pos)
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
    sources.push_back(CLSettings::get_source_position_update());
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
    // Output positions
    cl::Buffer cl_out_pos(context, CL_MEM_WRITE_ONLY, 3 * num_parts * sizeof(cl_double), &status);
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
    if (status != CL_SUCCESS)
    {
        switch (status)
        {
        case CL_INVALID_CONTEXT: std::cerr << "CL_INVALID_CONTEXT" << std::endl; break;
        case CL_INVALID_MEM_OBJECT: std::cerr << "CL_INVALID_MEM_OBJECT" << std::endl; break;
        case CL_INVALID_VALUE: std::cerr << "CL_INVALID_VALUE" << std::endl; break;
        case CL_INVALID_EVENT_WAIT_LIST: std::cerr << "CL_INVALID_EVENT_WAIT_LIST" << std::endl; break;
        case CL_MISALIGNED_SUB_BUFFER_OFFSET: std::cerr << "CL_MISALIGNED_SUB_BUFFER_OFFSET" << std::endl; break;
        case CL_MEM_OBJECT_ALLOCATION_FAILURE: std::cerr << "CL_MEM_OBJECT_ALLOCATION_FAILURE" << std::endl; break;
        case CL_OUT_OF_RESOURCES: std::cerr << "CL_OUT_OF_RESOURCES" << std::endl; break;
        case CL_OUT_OF_HOST_MEMORY: std::cerr << "CL_OUT_OF_HOST_MEMORY" << std::endl; break;

        default: std::cerr << "Unrecognized error: " << status << std::endl; break;
        }
        std::stringstream ss;
        ss << "Errors occurred while writing buffers on OpenCL device memory for positions update." << std::endl;
        throw new std::runtime_error(ss.str());
    }
    status = queue.enqueueWriteBuffer(cl_in_vel, CL_TRUE, 0, 3 * num_parts * sizeof(cl_double), in_vel);
    //status = queue.enqueueReadBuffer(cl_out_pos, CL_TRUE, 0, 3 * num_parts * sizeof(cl_double), out_pos);

    // Create the kernel and set the arguments
    static cl::Kernel kernel(program, POSITION_UPDATE_KERNEL_NAME, &status);
    if (!first_run && status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "Errors occurred while creating the OpenCL kernel for position update." << std::endl;
        throw new std::runtime_error(ss.str());
    }

    status = kernel.setArg(0, cl_in_pos);
    status |= kernel.setArg(1, cl_in_vel);
    status |= kernel.setArg(2, num_parts);
    status |= kernel.setArg(3, delta_time);
    status |= kernel.setArg(4, cl_out_pos);
    if (status != CL_SUCCESS)
    {
        std::stringstream ss;
        ss << "errors occurred while setting the arguments of the opencl kernel for position update." << std::endl;
        throw new std::runtime_error(ss.str());
    }

    // Execute the kernel
    queue.enqueueNDRangeKernel(kernel, cl::NullRange, cl::NDRange(num_parts), cl::NullRange);

    // Retrieve the results
    queue.enqueueReadBuffer(cl_out_pos, CL_TRUE, 0, 3 * num_parts * sizeof(cl_double), out_pos);
    queue.finish();

    first_run = true;
}