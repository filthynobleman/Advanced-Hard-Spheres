__kernel void wall_collision(__global const double* pos,
							 __global const double* vel,
							 __global const double* radii,
							 const ulong num_parts,
							 __global const double* x_wall,
							 __global const double* y_wall,
							 __global const double* z_wall,
							 __global double* delta_time,
							 __global int* axis)
{
	int i = get_global_id(0);
	if (i < num_parts)
	{
		double x = INFINITY;
		if (vel[3 * i] > 0)
			x = x_wall[1] - radii[i];
		else if (vel[3 * i] < 0)
			x = x_wall[0] + radii[i];
			
		double y = INFINITY;
		if (vel[3 * i + 1] > 0)
			y = y_wall[1] - radii[i];
		else if (vel[3 * i + 1] < 0)
			y = y_wall[0] + radii[i];
			
		double z = INFINITY;
		if (vel[3 * i + 2] > 0)
			z = z_wall[1] - radii[i];
		else if (vel[3 * i + 2] < 0)
			z = z_wall[0] + radii[i];
			

		double delta_x = (x - pos[3 * i]) / vel[3 * i];
		double delta_y = (y - pos[3 * i + 1]) / vel[3 * i + 1];
		double delta_z = (z - pos[3 * i + 2]) / vel[3 * i + 2];

		delta_time[i] = delta_x;
		axis[i] = 1;
		if (vel[3 * i] < 0)
			axis[i] = -1;
		if (delta_y < delta_time[i])
		{
			delta_time[i] = delta_y;
			axis[i] = 2;
			if (vel[3 * i + 1] < 0)
				axis[i] = -2;
		}
		if (delta_z < delta_time[i])
		{
			delta_time[i] = delta_z;
			axis[i] = 3;
			if (vel[3 * i + 2] < 0)
				axis[i] = -3;
		}
	}
}