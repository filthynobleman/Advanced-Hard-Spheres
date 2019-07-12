#include "CLSettings.h"
#include <fstream>
#include <sstream>
#include <iostream>

#define POSITION_UPDATE_KERNEL_SOURCE   "pos_update.cl"
#define WALL_COLLISION_KERNEL_SOURCE    "wall_collision.cl"
#define PART_COLLISION_KERNEL_SOURCE    "part_collision.cl"

cl::Device* CLSettings::_device;
std::string CLSettings::_pos_update_source;
std::string CLSettings::_wall_collision_source;
std::string CLSettings::_part_collision_source;
std::string CLSettings::_output_file;

cl_device_id CLSettings::select_device()
{
    cl::vector<cl::Platform> plats;
    cl::Platform::get(&plats);
    
    cl::vector<cl::Device> devs;
    for (size_t i = 0; i < plats.size(); i++)
    {
        cl::vector<cl::Device> devs_loc;
        plats[i].getDevices(CL_DEVICE_TYPE_ALL, &devs_loc);
        devs.insert(devs.end(), devs_loc.begin(), devs_loc.end());
    }

    for (size_t i = 0; i < devs.size(); i++)
    {
        std::string devname;
        devs[i].getInfo(CL_DEVICE_NAME, &devname);
        std::cout << (i + 1) << ". " << devname << std::endl;
    }

    int d;
    while (true)
    {
        std::cout << "Select a device to show details (insert 0 to skip): ";
        std::cin >> d;
        if (d < 0 || d > devs.size())
            continue;
        if (d == 0)
            break;
        else
        {
            std::string strval;
            devs[d - 1].getInfo(CL_DEVICE_NAME, &strval);
            std::cout << "Device Name: " << strval << std::endl;
            devs[d - 1].getInfo(CL_DEVICE_VENDOR, &strval);
            std::cout << "Device Vendor: " << strval << std::endl;
            devs[d - 1].getInfo(CL_DEVICE_VERSION, &strval);
            std::cout << "OpenCL Version: " << strval << std::endl;
            cl::vector<size_t> wi_sizes;
            devs[d - 1].getInfo(CL_DEVICE_MAX_WORK_ITEM_SIZES, &wi_sizes);
            std::cout << "Maximum Sizes of Work Items:" << std::endl;
            for (size_t i = 0; i < wi_sizes.size(); i++)
                std::cout << "  Size " << i << ": " << wi_sizes[i] << std::endl;
            std::cout << std::endl;
        }
    }

    while (true)
    {
        std::cout << "Select a device to use for the simulation: ";
        std::cin >> d;
        if (d <= 0 || d > devs.size())
            continue;
        return devs[d - 1].get();
    }
}

void CLSettings::set_device(cl::Device& device)
{
    _device = &device;
}

void CLSettings::set_output_file(std::string& filename)
{
    _output_file = std::string(filename);
}

cl::Device& CLSettings::get_device()
{
    return *_device;
}


std::string CLSettings::get_source_position_update()
{
    if (_pos_update_source.empty())
    {
        std::stringstream source;
        std::ifstream stream;
        stream.open(POSITION_UPDATE_KERNEL_SOURCE);
        if (!stream)
        {
            std::stringstream ss;
            ss << "Errors occurred while opening kernel source file " << POSITION_UPDATE_KERNEL_SOURCE << "." << std::endl;
            throw new std::runtime_error(ss.str());
        }
        source << stream.rdbuf();

        _pos_update_source = source.str();
    }

    return _pos_update_source;
}

std::string CLSettings::get_source_wall_collision()
{
    if (_wall_collision_source.empty())
    {
        std::stringstream source;
        std::ifstream stream;
        stream.open(WALL_COLLISION_KERNEL_SOURCE);
        if (!stream)
        {
            std::stringstream ss;
            ss << "Errors occurred while opening kernel source file " << WALL_COLLISION_KERNEL_SOURCE << "." << std::endl;
            throw new std::runtime_error(ss.str());
        }
        source << stream.rdbuf();

        _wall_collision_source = source.str();
    }

    return _wall_collision_source;
}

std::string CLSettings::get_source_part_collision()
{
    if (_part_collision_source.empty())
    {
        std::stringstream source;
        std::ifstream stream;
        stream.open(PART_COLLISION_KERNEL_SOURCE);
        if (!stream)
        {
            std::stringstream ss;
            ss << "Errors occurred while opening kernel source file " << PART_COLLISION_KERNEL_SOURCE << "." << std::endl;
            throw new std::runtime_error(ss.str());
        }
        source << stream.rdbuf();

        _part_collision_source = source.str();
    }

    return _part_collision_source;
}

std::string CLSettings::get_output_file()
{
    return _output_file;
}
