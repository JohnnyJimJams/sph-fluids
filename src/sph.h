#include <list>

using namespace std;

class SphFluidSolver;
struct Particle;
struct GridElement;

#ifndef _SPH_H_
#define _SPH_H_

#include "Vector.h"


#define PI_FLOAT				3.141592653589793f

#define SQR(x)					((x) * (x))
#define CUBE(x)					((x) * (x) * (x))
#define POW6(x)					(CUBE(x) * CUBE(x))
#define POW9(x)					(POW6(x) * CUBE(x))


struct Particle {
	int id;
	float mass;
	float density;
	Vector3f position;
	Vector3f velocity;
	Vector3f force;
	Vector3f color_gradient;
	float color_laplacian;

	Particle() { mass = 1.0f; }
};

struct GridElement {
	list<Particle> particles;
};

class SphFluidSolver {

	const int grid_width;
	const int grid_height;
	const int grid_depth;

	const float core_radius;
	const float gas_constant;
	const float mu;
	const float rest_density;
	const float point_damping;
	const float sigma;
	const float timestep;

	GridElement *grid_elements;
	GridElement *sleeping_grid_elements;

public:

	SphFluidSolver(
			float core_radius,
			float gas_constant,
			float mu,
			float rest_density,
			float point_damping,
			float sigma,
			float timestep)
	: grid_width(32),
	  grid_height(32),
	  grid_depth(32),
	  core_radius(core_radius),
	  gas_constant(gas_constant),
	  mu(mu),
	  rest_density(rest_density),
	  point_damping(point_damping),
	  sigma(sigma),
	  timestep(timestep) {

	}

	void update(void(*inter_hook)() = NULL, void(*post_hook)() = NULL);

	void init_particles(Particle *particles, int count);

	template <typename Function>
	void foreach_particle(Function function) {
		for (int k = 0; k < grid_depth; k++) {
			for (int j = 0; j < grid_height; j++) {
				for (int i = 0; i < grid_width; i++) {
					GridElement &grid_element = grid_elements[grid_width * (k * grid_height + j) + i];

					list<Particle> &plist = grid_element.particles;
					for (list<Particle>::iterator piter = plist.begin(); piter != plist.end(); piter++) {
						function(*piter);
					}
				}
			}
		}
	}

private:

	float kernel(const Vector3f &r, const float h);

	Vector3f gradient_kernel(const Vector3f &r, const float h);

	float laplacian_kernel(const Vector3f &r, const float h);

	Vector3f gradient_pressure_kernel(const Vector3f &r, const float h);

	float laplacian_viscosity_kernel(const Vector3f &r, const float h);

	void add_density(Particle &particle, Particle &neighbour);

	void sum_density(GridElement &grid_element, Particle &particle);

	void sum_all_density(int i, int j, int k, Particle &particle);

	void update_densities(int i, int j, int k);

	void add_forces(Particle &particle, Particle &neighbour);

	void sum_forces(GridElement &grid_element, Particle &particle);

	void sum_all_forces(int i, int j, int k, Particle &particle);

	void update_forces(int i, int j, int k);

	void update_particle(Particle &particle);

	void update_particles(int i, int j, int k);

	void reset_particle(Particle &particle);

	void reset_particles();

	void insert_into_grid(int i, int j, int k);

	void update_grid();

	void update_densities();

	void update_forces();

	void update_particles();

	GridElement &grid(int i, int j, int k);

	GridElement &sleeping_grid(int i, int j, int k);

	int grid_index(int i, int j, int k);

	void add_to_grid(GridElement *target_grid, Particle &particle);
};

#endif /* _SPH_H_ */

