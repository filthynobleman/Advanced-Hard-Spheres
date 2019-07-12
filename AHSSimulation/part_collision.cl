__kernel void part_collision(__global const double* pos,
							 __global const double* vel,
							 __global const double* radii,
							 const ulong num_parts,
							 __global double* delta_times)
{
	int i = get_global_id(0);
	int j = get_global_id(1);
	if (i < num_parts && j < num_parts)
	{
		// Velocities dot product
		double a = (vel[3 * i]     - vel[3 * j])     * (vel[3 * i]     - vel[3 * j])     +
				   (vel[3 * i + 1] - vel[3 * j + 1]) * (vel[3 * i + 1] - vel[3 * j + 1]) + 
				   (vel[3 * i + 2] - vel[3 * j + 2]) * (vel[3 * i + 2] - vel[3 * j + 2]);
		// Position-velocity dot product times two
		double b = 2 * ((pos[3 * i]     - pos[3 * j])     * (vel[3 * i]     - vel[3 * j])     +
					    (pos[3 * i + 1] - pos[3 * j + 1]) * (vel[3 * i + 1] - vel[3 * j + 1]) + 
					    (pos[3 * i + 2] - pos[3 * j + 2]) * (vel[3 * i + 2] - vel[3 * j + 2]));
		// Positions dot product, minus the square of the sum of radii
		double c = ((pos[3 * i]     - pos[3 * j])     * (pos[3 * i]     - pos[3 * j])     +
					(pos[3 * i + 1] - pos[3 * j + 1]) * (pos[3 * i + 1] - pos[3 * j + 1]) + 
					(pos[3 * i + 2] - pos[3 * j + 2]) * (pos[3 * i + 2] - pos[3 * j + 2]))
				   -
				   ((radii[i] + radii[j]) * (radii[i] + radii[j]));

		// Collision occurs only if b is negative
		if (b >= 0 || b*b < 4*a*c)
		{
			delta_times[j * num_parts + i] = INFINITY;
		}
		else
		{

			//printf("A(%d, %d) = %f\nB(%d, %d) = %f\nC(%d, %d) = %f\n", i, j, a, i, j, b, i, j, c);
			//printf("(%d, %d) ==> %f\n", i, j, (- b - sqrt(b*b - 4*a*c)) / (2 * a));

			// Solve the polynomial, if the relative difference between the centers is greater
			// than the sum of the radii
			if (c >= 0)
				delta_times[j * num_parts + i] = (- b - sqrt(b*b - 4*a*c)) / (2 * a);
			// Otherwise, give to the couples the highest priority, using a negative time
			else
				delta_times[j * num_parts + i] = -1;
		}
	}
}