#include "sph.h"
#include <iostream>
using namespace std;


const float		core_radius = 1.1f;
const float		gas_constant = 1000.0f;
const float		mu = 0.1f;
const float		rest_density = 1.2f;
const float		point_damping = 2.0f;
const float		sigma = 1.0f;

const float		timestep = 0.007f;

GridElement		*grid;
GridElement		*sleeping_grid;


inline float kernel(const Vector3f &r, const float h)
{
	if (length(r) > h)
		return (0.0f);

	return (315.0f / (64.0f * PI_FLOAT * POW9(h)) * CUBE(SQR(h) - dot(r, r)));
}

inline Vector3f gradient_kernel(const Vector3f &r, const float h)
{
	if (length(r) > h)
		return (Vector3f(0.0f, 0.0f, 0.0f));

	return (-945.0f / (32.0f * PI_FLOAT * POW9(h)) * SQR(SQR(h) - dot(r, r)) * r);
}

inline float laplacian_kernel(const Vector3f &r, const float h)
{
	if (length(r) > h)
		return (0.0f);

	return (  945.0f / (32.0f * PI_FLOAT * POW9(h))
	        * (SQR(h) - dot(r, r)) * (7.0f * dot(r, r) - 3.0f * SQR(h)));
}

inline Vector3f gradient_pressure_kernel(const Vector3f &r, const float h)
{
	if ((length(r) > h) || (length(r) < 0.001f))
		return (Vector3f(0.0f, 0.0f, 0.0f));

	return (-45.0f / (PI_FLOAT * POW6(h)) * SQR(h - length(r)) * normalize(r));
}

inline float laplacian_viscosity_kernel(const Vector3f &r, const float h)
{
	if (length(r) > h)
		return (0.0f);

	return (45.0f / (PI_FLOAT * POW6(h)) * (h - length(r)));
}

inline void add_density(Particle &particle, Particle &neighbour)
{
	Vector3f	r;
	float		common;

	r = particle.position - neighbour.position;

	common = kernel(r, core_radius);
	particle.density  += neighbour.mass * common;
///	neighbour.density += particle.mass * common;
}

void sum_density(GridElement &grid_element, Particle &particle)
{
	Vector3f			r;
	Particle			*neighbour;
	ParticleList		*plist;
	ParticleListIter	piter;

	/* Handle the particles directly stored in the grid. */
	neighbour = grid_element.particles;
	for (int x = 0; x < grid_element.particle_count; x++)
		add_density(particle, neighbour[x]);

	/* Handle the particles not stored in the grid. */
	plist = &grid_element.overflow_particles;
	for (piter = plist->begin(); piter != plist->end(); piter++)
		add_density(particle, *piter);
}

inline void sum_all_density(int i, int j, int k, Particle &particle)
{
	particle.density = 0.0f;

	for (int x = i - 1; x <= i + 1; x++)
		for (int y = j - 1; y <= j + 1; y++)
			for (int z = k - 1; z <= k + 1; z++)
			{
				if ((x < 0) || (x >= GRID_WIDTH))
					continue;
				if ((y < 0) || (y >= GRID_HEIGHT))
					continue;
				if ((z < 0) || (z >= GRID_DEPTH))
					continue;

				sum_density(GRID(x, y, z), particle);
			}
}

void update_density(int i, int j, int k)
{
	GridElement			*grid_element;
	Particle			*particles;
	ParticleList		*plist;
	ParticleListIter	piter;

	grid_element = &GRID(i, j, k);

	/* Handle the particles directly stored in the grid. */
	particles = grid_element->particles;
	for (int x = 0; x < grid_element->particle_count; x++)
		sum_all_density(i, j, k, particles[x]);

	/* Handle the particles not stored in the grid. */
	plist = &grid_element->overflow_particles;
	for (piter = plist->begin(); piter != plist->end(); piter++)
		sum_all_density(i, j, k, *piter);
}

inline void add_forces(Particle &particle, Particle &neighbour)
{
	Vector3f	r;
	Vector3f	common;
	float		value;

	if (&particle == &neighbour)
		return;

	r = particle.position - neighbour.position;

	/* Compute the pressure force. */
	common = 0.5f * gas_constant * (  (particle.density - rest_density)
	                                + (neighbour.density - rest_density))
	         * gradient_pressure_kernel(r, core_radius);
	particle.force  += -neighbour.mass / neighbour.density * common;
	particle.p_force  += -neighbour.mass / neighbour.density * common;
///	neighbour.force -= -particle.mass / particle.density * common;

	/* Compute the viscosity force. */
	common = mu * (neighbour.velocity - particle.velocity)
	         * laplacian_viscosity_kernel(r, core_radius);
	particle.force  += neighbour.mass / neighbour.density * common;
	particle.v_force  += neighbour.mass / neighbour.density * common;
///	neighbour.force -= particle.mass / particle.density * common;

	/* Compute the gradient of the color field. */
	common = gradient_kernel(r, core_radius);
	particle.color_gradient  += neighbour.mass / neighbour.density * common;
///	neighbour.color_gradient -= particle.mass / particle.density * common;

	/* Compute the gradient of the color field. */
	value = laplacian_kernel(r, core_radius);
	particle.color_laplacian  += neighbour.mass / neighbour.density * value;
///	neighbour.color_laplacian -= particle.mass / particle.density * value;
}

void sum_forces(GridElement &grid_element, Particle &particle)
{
	Vector3f			r;
	Vector3f			common;
	Particle			*neighbour;
	ParticleList		*plist;
	ParticleListIter	piter;

	/* Handle the particles directly stored in the grid. */
	neighbour = grid_element.particles;
	for (int x = 0; x < grid_element.particle_count; x++)
		add_forces(particle, neighbour[x]);

	/* Handle the particles not stored in the grid. */
	plist = &grid_element.overflow_particles;
	for (piter = plist->begin(); piter != plist->end(); piter++)
		add_forces(particle, *piter);
}

void sum_all_forces(int i, int j, int k, Particle &particle)
{
	particle.force = Vector3f(0.0f, 0.0f, 0.0f);
	particle.p_force = Vector3f(0.0f, 0.0f, 0.0f);
	particle.v_force = Vector3f(0.0f, 0.0f, 0.0f);
	particle.color_gradient = Vector3f(0.0f, 0.0f, 0.0f);
	particle.color_laplacian = 0.0f;

	for (int x = i - 1; x <= i + 1; x++)
		for (int y = j - 1; y <= j + 1; y++)
			for (int z = k - 1; z <= k + 1; z++)
			{
				if ((x < 0) || (x >= GRID_WIDTH))
					continue;
				if ((y < 0) || (y >= GRID_HEIGHT))
					continue;
				if ((z < 0) || (z >= GRID_DEPTH))
					continue;

				sum_forces(GRID(x, y, z), particle);
			}
}

void update_forces(int i, int j, int k)
{
	GridElement			*grid_element;
	Particle			*particles;
	ParticleList		*plist;
	ParticleListIter	piter;

	grid_element = &GRID(i, j, k);

	/* Handle the particles directly stored in the grid. */
	particles = grid_element->particles;
	for (int x = 0; x < grid_element->particle_count; x++)
		sum_all_forces(i, j, k, particles[x]);

	/* Handle the particles not stored in the grid. */
	plist = &grid_element->overflow_particles;
	for (piter = plist->begin(); piter != plist->end(); piter++)
		sum_all_forces(i, j, k, *piter);
}

inline void update_velocity(Particle &particle)
{
	Vector3f	acceleration;

	if (length(particle.color_gradient) > 0.001f)
	{
		particle.force +=   -sigma * particle.color_laplacian
		                  * normalize(particle.color_gradient);
		particle.t_force +=   -sigma * particle.color_laplacian
		                    * normalize(particle.color_gradient);
	}

	acceleration =   particle.force / particle.density
	               - point_damping * particle.velocity / particle.mass;

	particle.velocity += timestep * acceleration;
}

inline void update_position(Particle &particle)
{
	particle.position += timestep * particle.velocity;
}

void update_particles(int i, int j, int k)
{
	GridElement			*grid_element;
	Particle			*particles;
	ParticleList		*plist;
	ParticleListIter	piter;

	grid_element = &GRID(i, j, k);

	/* Handle the particles directly stored in the grid. */
	particles = grid_element->particles;
	for (int x = 0; x < grid_element->particle_count; x++)
	{
		update_velocity(particles[x]);
		update_position(particles[x]);
	}

	/* Handle the particles not stored in the grid. */
	plist = &grid_element->overflow_particles;
	for (piter = plist->begin(); piter != plist->end(); piter++)
	{
		update_velocity(*piter);
		update_position(*piter);
	}
}

inline void insert_into_grid(int i, int j, int k)
{
	GridElement			*grid_element;
	Particle			*particles;
	ParticleList		*plist;
	ParticleListIter	piter;

	grid_element = &GRID(i, j, k);

	/* Handle the particles directly stored in the grid. */
	particles = grid_element->particles;
	for (int x = 0; x < grid_element->particle_count; x++)
		ADD_TO_GRID(sleeping_grid, particles[x]);

	/* Handle the particles not stored in the grid. */
	plist = &grid_element->overflow_particles;
	for (piter = plist->begin(); piter != plist->end(); piter++)
		ADD_TO_GRID(sleeping_grid, *piter);
}

void update_grid()
{
	GridElement	*temp;

	for (int i = 0; i < GRID_WIDTH; i++)
		for (int j = 0; j < GRID_HEIGHT; j++)
			for (int k = 0; k < GRID_DEPTH; k++)
			{
				insert_into_grid(i, j, k);
				GRID(i, j, k).particle_count = 0;
				GRID(i, j, k).overflow_particles.clear();
			}

	/* Swap the grids. */
	temp = grid;
	grid = sleeping_grid;
	sleeping_grid = temp;
}

void update(void(*inter_hook)() = NULL, void(*post_hook)() = NULL)
{
	int		i, j, k;

	for (i = 0; i < GRID_WIDTH; i++)
		for (j = 0; j < GRID_HEIGHT; j++)
			for (k = 0; k < GRID_DEPTH; k++)
				update_density(i, j, k);

	for (i = 0; i < GRID_WIDTH; i++)
		for (j = 0; j < GRID_HEIGHT; j++)
			for (k = 0; k < GRID_DEPTH; k++)
				update_forces(i, j, k);

	/* User supplied hook, e.g. for adding custom forces (gravity, ...). */
	if (inter_hook != NULL)
		inter_hook();

	for (i = 0; i < GRID_WIDTH; i++)
		for (j = 0; j < GRID_HEIGHT; j++)
			for (k = 0; k < GRID_DEPTH; k++)
				update_particles(i, j, k);

	/* User supplied hook, e.g. for handling collisions. */
	if (post_hook != NULL)
		post_hook();

	update_grid();
}

void init_particles(Particle *particles, int count)
{
	grid = new GridElement[GRID_WIDTH * GRID_HEIGHT * GRID_DEPTH];
	sleeping_grid = new GridElement[GRID_WIDTH * GRID_HEIGHT * GRID_DEPTH];

	for (int x = 0; x < count; x++)
		ADD_TO_GRID(grid, particles[x]);
}

