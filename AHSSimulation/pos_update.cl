__kernel void pos_update(__global const double* pos, 
						 __global const double* vel,
						 const ulong num_parts,
						 const double delta_time,
						 __global double* out_pos)
{
	int i = get_global_id(0);
	if (i < num_parts)
	{
		out_pos[3 * i] = pos[3 * i] + delta_time * vel[3 * i];
		out_pos[3 * i + 1] = pos[3 * i + 1] + delta_time * vel[3 * i + 1];
		out_pos[3 * i + 2] = pos[3 * i + 2] + delta_time * vel[3 * i + 2];
	}
}