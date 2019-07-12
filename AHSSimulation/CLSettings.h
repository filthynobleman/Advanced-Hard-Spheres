#pragma once

#include <CL/cl2.hpp>

#define POSITION_UPDATE_KERNEL_NAME "pos_update"
#define WALL_COLLISION_KERNEL_NAME  "wall_collision"
#define PART_COLLISION_KERNEL_NAME  "part_collision"

#define SIMULATION_TYPE_INELSATIC   (size_t)0;
#define SIMULATION_TYPE_FUSION      (size_t)1;
#define SIMULATION_TYPE_FISSION     (size_t)2;

class CLSettings
{
private:
    static cl::Device* _device;
    static std::string _pos_update_source;
    static std::string _wall_collision_source;
    static std::string _part_collision_source;
    static std::string _output_file;

    CLSettings() {};
    CLSettings(CLSettings& cls) {};
    ~CLSettings() {};
    void operator=(CLSettings& cls) {};

public:
    static cl_device_id select_device();
    static void set_device(cl::Device& device);
    static void set_output_file(std::string& filename);
    static cl::Device& get_device();
    static std::string get_source_position_update();
    static std::string get_source_wall_collision();
    static std::string get_source_part_collision();
    static std::string get_output_file();
};