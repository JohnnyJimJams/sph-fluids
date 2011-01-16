#include <iostream>
#include <sys/time.h>
#include <cstdio>
#include <GL/glut.h>

#include "sph.h"

using namespace std;

#define WIDTH		15
#define HEIGHT		15
#define DEPTH		15


int wndWidth = 700, wndHeight = 700;

int oldX, oldY, rotX = 0, rotY = 0, zoomZ = 0;
int oldTransX, oldTransY, transX = 0, transY = 0;
bool zoom, trans;

GLuint sphereId;

int simulationSteps = 2;

GLfloat rotation_matrix[16];

Vector3f gravity_direction;


extern void init_particles(Particle *particles, int count);
extern void update(void(*inter_hook)() = NULL, void(*post_hook)() = NULL);

const int particle_count = 2000;


void init_liquid()
 {
	Particle *particles = new Particle[particle_count];

	int count = particle_count;
	Particle *particle_iter = particles;
	while (true) {
		for (int j = 0; j < HEIGHT; j++) {
			for (int k = 0; k < DEPTH; k++) {
				for (int i = 0; i < WIDTH; i++) {
					if (count-- == 0) {
						init_particles(particles, particle_count);
						return;
					}

					particle_iter->position.x = i;
					particle_iter->position.y = j;
					particle_iter->position.z = k;
					particle_iter++;
				}
			}
		}
	}
}

void draw_particle(Particle &particle)
{
	glTranslatef(+particle.position.x, +particle.position.y, +particle.position.z);
	glCallList(sphereId);
#if 0
	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
	glColor3ub(0, 0, 255);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glColor3ub(255, 0, 0);
	pos = normalize(particle.force);
	glVertex3fv((GLfloat *)&pos);
	glEnd();
	glEnable(GL_LIGHTING);
#endif
	glTranslatef(-particle.position.x, -particle.position.y, -particle.position.z);
}

void add_gravity_force(Particle &particle)
{
	particle.force += 100.0f * gravity_direction * particle.density / particle.mass;
}

void add_global_forces()
{
	foreach_particle(add_gravity_force);
}

void handle_particle_collision_cylinder(Particle &particle)
{
	Vector3f mid = Vector3f(WIDTH, 0.0f, DEPTH) / 2.0f;
	Vector3f distance = Vector3f(particle.position.x, 0.0f, particle.position.z) - mid;

	if (length(distance) >= WIDTH / 2)
	{
		distance = normalize(distance);

		particle.position.x = (mid + WIDTH / 2 * distance).x;
		particle.position.z = (mid + WIDTH / 2 * distance).z;

		particle.velocity -= 2.0f * dot(particle.velocity, distance) * distance;
	}

	if (particle.position.y >= HEIGHT - 1)
	{
		particle.position.y = HEIGHT - 1;
		particle.velocity.y *= -1.0f;
	}
	else if (particle.position.y < 0.0f)
	{
		particle.position.y = 0.0f;
		particle.velocity.y *= -1.0f;
	}
}

void handle_particle_collision_cube(Particle &particle)
{
	float &px = particle.position.x;
	float &py = particle.position.y;
	float &pz = particle.position.z;

	float &vx = particle.velocity.x;
	float &vy = particle.velocity.y;
	float &vz = particle.velocity.z;

	if (px < 0 || px > WIDTH) {
		px = min(max(px, 0.0f), (float) WIDTH);
		vx *= -1.0f;
	}
	if (py < 0 || py > HEIGHT) {
		py = min(max(py, 0.0f), (float) HEIGHT);
		vy *= -1.0f;
	}
	if (pz < 0 || pz > DEPTH) {
		pz = min(max(pz, 0.0f), (float) DEPTH);
		vz *= -1.0f;
	}
}

void handle_collisions()
{
#if 0
	foreach_particle(handle_particle_collision_cylinder);
#else
	foreach_particle(handle_particle_collision_cube);
#endif
}

void extract_gravity_direction()
{
	gravity_direction.x = -rotation_matrix[1];
	gravity_direction.y = -rotation_matrix[5];
	gravity_direction.z = -rotation_matrix[9];
}

void init()
{
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	glEnable(GL_DEPTH_TEST);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);

	sphereId = glGenLists(1);
	glNewList(sphereId, GL_COMPILE);
	glutSolidSphere(0.3, 20, 20);
	glEndList();
}

void reshape(int width, int height)
{
	glViewport(0, 0, width, height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//glOrtho(-0.1 * WIDTH, 1.1 * WIDTH, +0.1 * HEIGHT, -1.1 * HEIGHT, 1.0, 1000.0);
	gluPerspective(45.0, 1.0, 1.0, 1000.0);

	wndWidth = width;
	wndHeight = height;
}

void mouse(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON)
	{
		if (glutGetModifiers() & GLUT_ACTIVE_CTRL)
		{
			trans = true;
			oldTransX = x - transX;
			oldTransY = y - transY;
		}
		else
		{
			trans = false;
			oldX = x;
			oldY = y;
		}
	}
	else if (button == GLUT_RIGHT_BUTTON)
	{
		zoom = !zoom;
		oldY = y - zoomZ;
	}
}

void motion(int x, int y)
{
	if (!zoom)
	{
		if (trans)
		{
			transX = x - oldTransX;
			transY = y - oldTransY;
		}
		else
		{
			rotY = x - oldX;
			oldX = x;
			rotX = y - oldY;
			oldY = y;
		}
	}
	else
		zoomZ = y - oldY;

	glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'q':
	case 'Q':
	case 0x1bU: /* ESC */
		exit(0);
	default:
		break;
	}
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glGetFloatv(GL_MODELVIEW_MATRIX, rotation_matrix);

	/* Handle rotations and translations separately. */
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glRotatef(rotY, 0.0f, 1.0f, 0.0f);
	glRotatef(rotX, 1.0f, 0.0f, 0.0f);
	glMultMatrixf(rotation_matrix);

	/* Save the rotational component of the modelview matrix. */
	glPushMatrix();
	glLoadIdentity();

	gluLookAt(0.0, 0.0, 30.0,
	          0.0, 0.0, 0.0,
	          0.0, 1.0, 0.0);

	glTranslatef(0.0f, 0.0f, (float)zoomZ / 20.0f);
	glRotatef(rotY, 0.0f, 1.0f, 0.0f);
	glRotatef(rotX, 1.0f, 0.0f, 0.0f);
	glMultMatrixf(rotation_matrix);
	glTranslatef(-WIDTH / 2.0f, -HEIGHT / 2.0f, -DEPTH / 2.0f);

	rotX = rotY = 0;

	extract_gravity_direction();

	timeval tv1, tv2;
	gettimeofday(&tv1, NULL);

	for (int i = 0; i < simulationSteps; ++i) {
		update(add_global_forces, handle_collisions);
	}

	gettimeofday(&tv2, NULL);
	int simulationTime = 1000 * (tv2.tv_sec - tv1.tv_sec) + (tv2.tv_usec - tv1.tv_usec) / 1000;
	printf("TIME[simulation]      : %d ms\n", simulationTime);

	gettimeofday(&tv1, NULL);
	foreach_particle(draw_particle);
	gettimeofday(&tv2, NULL);
	int renderingTime = 1000 * (tv2.tv_sec - tv1.tv_sec) + (tv2.tv_usec - tv1.tv_usec) / 1000;
	printf("TIME[rendering]       : %d ms\n", renderingTime);
	printf("\n");

	glutSwapBuffers();

	/* Restore the rotational component of the modelview matrix. */
	glPopMatrix();
}

void idle()
{
	glutPostRedisplay();
}

void print_usage()
{
	cout << endl;
	cout << "KEYSTROKE       ACTION" << endl;
	cout << "=========       ======" << endl << endl;
	cout << "q, Q, <ESC>     exit the program" << endl;
	cout << endl;
}

int main(int argc, char *argv[])
{
	print_usage();

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
	glutInitWindowSize(wndWidth, wndHeight);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("SPH Fluids");

	glutReshapeFunc(reshape);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);
	glutIdleFunc(idle);

	init();
	init_liquid();

	glutMainLoop();

	return (0);
}

