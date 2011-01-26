#include <list>

using namespace std;

class SphFluidSolver;
struct Particle;
struct GridElement;

#ifndef _SPH_H_
#define _SPH_H_

#include "Vector.h"


struct Particle {

	int id;
	float mass;
	float density;
	Vector3f position;
	Vector3f velocity;
	Vector3f force;
	Vector3f color_gradient;
	float color_laplacian;

	Vector3f viscosity_force;
	Vector3f pressure_force;

	Particle() { mass = 1.0f; }
};

struct GridElement {
	list<Particle> particles;
};

struct FluidMaterial {

	const float gas_constant;
	const float mu;
	const float rest_density;
	const float sigma;
	const float point_damping;

	FluidMaterial(
			float gas_constant,
			float mu,
			float rest_density,
			float sigma,
			float point_damping)
	  : gas_constant(gas_constant),
	    mu(mu),
	    rest_density(rest_density),
	    sigma(sigma),
	    point_damping(point_damping) {
	}
};

class SphFluidSolver {

	const int grid_width;
	const int grid_height;
	const int grid_depth;

	const float core_radius;
	const float timestep;

	const FluidMaterial material;

	GridElement *grid_elements;
	GridElement *sleeping_grid_elements;

public:

	SphFluidSolver(
			float domain_width,
			float domain_height,
			float domain_depth,
			float core_radius,
			float timestep,
			FluidMaterial material)
	: grid_width((int) (domain_width / core_radius) + 1),
	  grid_height((int) (domain_height / core_radius) + 1),
	  grid_depth((int) (domain_depth / core_radius) + 1),
	  core_radius(core_radius),
	  timestep(timestep),
	  material(material) {

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

