#include "sph.h"

#include <iostream>

using namespace std;


const float core_radius = 1.1f;
const float gas_constant = 1000.0f;
const float mu = 0.1f;
const float rest_density = 1.2f;
const float point_damping = 2.0f;
const float sigma = 1.0f;

const float timestep = 0.02f;

GridElement *grid;
GridElement *sleeping_grid;


inline float kernel(const Vector3f &r, const float h) {
	if (length(r) > h) {
		return 0.0f;
	}

	return 315.0f / (64.0f * PI_FLOAT * POW9(h)) * CUBE(SQR(h) - dot(r, r));
}

inline Vector3f gradient_kernel(const Vector3f &r, const float h) {
	if (length(r) > h) {
		return Vector3f(0.0f, 0.0f, 0.0f);
	}

	return -945.0f / (32.0f * PI_FLOAT * POW9(h)) * SQR(SQR(h) - dot(r, r)) * r;
}

inline float laplacian_kernel(const Vector3f &r, const float h){
	if (length(r) > h) {
		return 0.0f;
	}

	return   945.0f / (32.0f * PI_FLOAT * POW9(h))
	       * (SQR(h) - dot(r, r)) * (7.0f * dot(r, r) - 3.0f * SQR(h));
}

inline Vector3f gradient_pressure_kernel(const Vector3f &r, const float h) {
	if ((length(r) > h) || (length(r) < 0.001f)) {
		return Vector3f(0.0f, 0.0f, 0.0f);
	}

	return (-45.0f / (PI_FLOAT * POW6(h)) * SQR(h - length(r)) * normalize(r));
}

inline float laplacian_viscosity_kernel(const Vector3f &r, const float h) {
	if (length(r) > h) {
		return 0.0f;
	}

	return (45.0f / (PI_FLOAT * POW6(h)) * (h - length(r)));
}

inline void add_density(Particle &particle, Particle &neighbour) {
	Vector3f r = particle.position - neighbour.position;
	particle.density += neighbour.mass * kernel(r, core_radius);
}

void sum_density(GridElement &grid_element, Particle &particle) {
	list<Particle> &plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) {
		add_density(particle, *piter);
	}
}

inline void sum_all_density(int i, int j, int k, Particle &particle) {
	particle.density = 0.0f;

	for (int x = i - 1; x <= i + 1; x++) {
		for (int y = j - 1; y <= j + 1; y++) {
			for (int z = k - 1; z <= k + 1; z++) {
				if (   (x < 0) || (x >= GRID_WIDTH)
					|| (y < 0) || (y >= GRID_HEIGHT)
					|| (z < 0) || (z >= GRID_DEPTH)) {
					continue;
				}

				sum_density(GRID(x, y, z), particle);
			}
		}
	}
}

void update_densities(int i, int j, int k) {
	GridElement &grid_element = GRID(i, j, k);

	list<Particle> &plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) {
		sum_all_density(i, j, k, *piter);
	}
}

inline void add_forces(Particle &particle, Particle &neighbour) {
	if (&particle == &neighbour) {
		return;
	}

	Vector3f r = particle.position - neighbour.position;

	/* Compute the pressure force. */
	Vector3f common = 0.5f * gas_constant * (  (particle.density - rest_density)
	                                + (neighbour.density - rest_density))
	         * gradient_pressure_kernel(r, core_radius);
	particle.force += -neighbour.mass / neighbour.density * common;

	/* Compute the viscosity force. */
	common = mu * (neighbour.velocity - particle.velocity)
	         * laplacian_viscosity_kernel(r, core_radius);
	particle.force += neighbour.mass / neighbour.density * common;

	/* Compute the gradient of the color field. */
	common = gradient_kernel(r, core_radius);
	particle.color_gradient += neighbour.mass / neighbour.density * common;

	/* Compute the gradient of the color field. */
	float value = laplacian_kernel(r, core_radius);
	particle.color_laplacian += neighbour.mass / neighbour.density * value;
}

void sum_forces(GridElement &grid_element, Particle &particle) {
	list<Particle>  &plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) {
		add_forces(particle, *piter);
	}
}

void sum_all_forces(int i, int j, int k, Particle &particle) {
	particle.force = Vector3f(0.0f, 0.0f, 0.0f);
	particle.color_gradient = Vector3f(0.0f, 0.0f, 0.0f);
	particle.color_laplacian = 0.0f;

	for (int x = i - 1; x <= i + 1; x++) {
		for (int y = j - 1; y <= j + 1; y++) {
			for (int z = k - 1; z <= k + 1; z++) {
				if (   (x < 0) || (x >= GRID_WIDTH)
					|| (y < 0) || (y >= GRID_HEIGHT)
					|| (z < 0) || (z >= GRID_DEPTH)) {
					continue;
				}

				sum_forces(GRID(x, y, z), particle);
			}
		}
	}
}

void update_forces(int i, int j, int k) {
	GridElement &grid_element = GRID(i, j, k);
	list<Particle>&plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) {
		sum_all_forces(i, j, k, *piter);
	}
}

inline void update_velocity(Particle &particle) {
	if (length(particle.color_gradient) > 0.001f) {
		particle.force +=   -sigma * particle.color_laplacian
		                  * normalize(particle.color_gradient);
	}

	Vector3f acceleration =   particle.force / particle.density
	               - point_damping * particle.velocity / particle.mass;

	particle.velocity += timestep * acceleration;
}

inline void update_position(Particle &particle) {
	particle.position += timestep * particle.velocity;
}

void update_particles(int i, int j, int k) {
	GridElement &grid_element = GRID(i, j, k);

	list<Particle> &plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) {
		update_velocity(*piter);
		update_position(*piter);
	}
}

inline void insert_into_grid(int i, int j, int k) {
	GridElement &grid_element = GRID(i, j, k);

	list<Particle> &plist = grid_element.particles;
	for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) {
		ADD_TO_GRID(sleeping_grid, *piter);
	}
}

void update_grid() {
	for (int i = 0; i < GRID_WIDTH; i++) {
		for (int j = 0; j < GRID_HEIGHT; j++) {
			for (int k = 0; k < GRID_DEPTH; k++) {
				insert_into_grid(i, j, k);
				GRID(i, j, k).particles.clear();
			}
		}
	}

	/* Swap the grids. */
	swap(grid, sleeping_grid);
}

void update_densities() {
    for (int i = 0; i < GRID_WIDTH; i++) {
		for (int j = 0; j < GRID_HEIGHT; j++) {
			for (int k = 0; k < GRID_DEPTH; k++) {
				update_densities(i, j, k);
			}
		}
	}
}

void update_forces() {
    for (int i = 0; i < GRID_WIDTH; i++) {
		for (int j = 0; j < GRID_HEIGHT; j++) {
			for (int k = 0; k < GRID_DEPTH; k++) {
				update_forces(i, j, k);
			}
		}
	}
}

void update_particles() {
    for (int i = 0; i < GRID_WIDTH; i++) {
		for (int j = 0; j < GRID_HEIGHT; j++) {
			for (int k = 0; k < GRID_DEPTH; k++) {
				update_particles(i, j, k);
			}
		}
	}
}

void update(void(*inter_hook)() = NULL, void(*post_hook)() = NULL) {
    update_densities();
    update_forces();

    /* User supplied hook, e.g. for adding custom forces (gravity, ...). */
	if (inter_hook != NULL) {
		inter_hook();
	}

    update_particles();

    /* User supplied hook, e.g. for handling collisions. */
	if (post_hook != NULL) {
		post_hook();
	}

	update_grid();
}

void init_particles(Particle *particles, int count) {
	grid = new GridElement[GRID_WIDTH * GRID_HEIGHT * GRID_DEPTH];
	sleeping_grid = new GridElement[GRID_WIDTH * GRID_HEIGHT * GRID_DEPTH];

	for (int x = 0; x < count; x++) {
		ADD_TO_GRID(grid, particles[x]);
	}
}

