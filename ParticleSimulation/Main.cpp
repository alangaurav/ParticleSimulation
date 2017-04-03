#include<stdio.h>
#include<stdlib.h>
#include<conio.h>
#include<math.h>
#include<mpi.h>
#include<iostream>

using namespace std;

#define PI 3.1415926535897932384 
#define dr 0.0001
#define dt 0.0001

/*Considering a two dimensional space*/

/* Steps to do
1. Initialize all charges and positions --done
2. Find Electric force on one particle --done
3. Find Magnetic force on one particle --done
4. Find Total force --done
5. f = ma, get acceleration --done
6. s = ut + 0.5at, get new position --done
7. v = u + at, get new velocity --done
8. Parallelize it --done
*/

const double POS_CHARGE = 1.60217662 * pow(10, -19); //in Coloumbs
const double NEG_CHARGE = -1.60217662 * pow(10, -19); //in Coloumbs
const double PROTON_MASS = 1.6726219 * pow(10, -27); //in Kilos
const double ELECTRON_MASS = 9.10938356 * pow(10, -31);// in Kilos
const double PERMEABILITY_CONSTANT = 4 * PI*pow(10.0, -7.0);
const double  COULOMB_CONSTANT = 9 * pow(10, 9);

typedef struct pos {
	double x, y;
}Position;

typedef struct vel {
	double x, y;
}Velocity;

typedef struct force {
	double x, y;
}Force;

typedef struct acc {
	double x, y;
}Acceleration;

typedef struct par {
	Position position;
	Velocity velocity;
	double mass;
	double charge;
}Particle;

Particle *list;
Force magfield;
Force elecfield;
int numberOfParticles;

void updatePosition(int index);
void updateVelocity(int index, Acceleration acceleration);
Acceleration getAcceleration(int index);
Force getTotalForce(int index);
Force getElecForce(int index);
Force getMagForce(int index);
Force getExtForce(int index);

/*Initialize the particle values*/

void initializeParticle(Particle *particle, Position p, int type) {
	particle->position.x = p.x;
	particle->position.y = p.y;
	particle->velocity.x = 0;
	particle->velocity.y = 0;
	if (type > 0) {
		particle->charge = POS_CHARGE;
		particle->mass = PROTON_MASS;
	}
	else {
		particle->charge = NEG_CHARGE;
		particle->mass = ELECTRON_MASS;
	}
}

/*Updates the position of a particle s = ut + 0.5at^2*/

void updatePosition(int index)
{
	Acceleration acceleration = getAcceleration(index);
	list[index].position.x = list[index].velocity.x * dt + 0.5 * acceleration.x * pow(dt, 2);
	list[index].position.y = list[index].velocity.y * dt + 0.5 * acceleration.y * pow(dt, 2);
	updateVelocity(index, acceleration);
}

/*Updates the velocity of a particle v = u + at*/

void updateVelocity(int index, Acceleration acceleration)
{
	list[index].velocity.x = list[index].velocity.x + acceleration.x * dt;
	list[index].velocity.y = list[index].velocity.y + acceleration.y * dt;
}

/*Returns the acceleration of the particle at any point of time*/

Acceleration getAcceleration(int index)
{
	Acceleration returnvalue;
	Force totalforce = getTotalForce(index);
	returnvalue.x = totalforce.x / list[index].mass;
	returnvalue.y = totalforce.y / list[index].mass;
	return returnvalue;
}

/*Returns the Total force acting on the particles*/

Force getTotalForce(int index)
{
	Force elecforce = getElecForce(index);
	Force magforce = getMagForce(index);
	Force extforce = getExtForce(index);
	Force returnvalue;
	returnvalue.x = elecforce.x + magforce.x + extforce.x;
	returnvalue.y = elecforce.y + magforce.y + extforce.y;
	return returnvalue;
}

/*Returns the Particle on particle force experienced by particle index*/

Force getElecForce(int index)
{
	Force returnvalue;
	returnvalue.x = 0;
	returnvalue.y = 0;
	for (int i = 0; i<numberOfParticles; ++i)
	{
		if (i == index)
			continue;

		double magnitude = PERMEABILITY_CONSTANT*list[index].charge*list[i].charge;
		double xdistance = pow(list[i].position.x - list[index].position.x, 2);
		double ydistance = pow(list[i].position.y - list[index].position.y, 2);
		double distance = xdistance + ydistance;
		distance = pow(distance, 0.5);
		magnitude = magnitude / distance;
		returnvalue.x = returnvalue.x + (xdistance / magnitude);
		returnvalue.y = returnvalue.y + (ydistance / magnitude);
	}
	return returnvalue;
}

/*Returns the Force acting on the particle index due to magnetic field*/

Force getMagForce(int index)
{
	Force returnvalue;
	returnvalue.x = list[index].velocity.x * magfield.x;
	returnvalue.y = list[index].velocity.y * magfield.y;
	return returnvalue;
}

/*Returns the Force acting on the particle index due to electric field*/

Force getExtForce(int index)
{
	Force returnvalue;
	returnvalue.x = list[index].charge * elecfield.x;
	returnvalue.y = list[index].charge * elecfield.y;
	return returnvalue;
}

/*Function to initialize magnetic field*/
void getMagneticField() {
	/*Using a generic magnetic field equation based on Biosarvats law and ampieres law*/
	magfield.x = 1.245261;
	magfield.y = 0.82791;
}

/*Function to initialize electric field*/
void getElectricField() {
	/*Using a generic electric field equation based on dipole intercations creating an isolated electric field per particle*/
	elecfield.x = 1.67182;
	elecfield.y = -1.2832;
}

int main(int argc, char *argv[]) {

	int worldRank, worldSize;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

	numberOfParticles = worldSize;
	list = (Particle *)malloc(worldSize * sizeof(Particle));

	getMagneticField();
	getElectricField();

	if (worldRank == 0) {

		/*Get values and create list of particles*/

		for (int i = 0; i<worldSize; i++)
		{
			Position p;
			int type;
			cout << "Enter the particle " << i+1 << "'s type (1/-1): \n";
			cin >> type;
			cout << "Enter the initial position of particle " << i + 1 << " (x,y)" << endl;
			cin >> p.x >> p.y;
			initializeParticle(&list[i], p, type);
		}

	}

	/*Send the initial particle properties to all processes*/

	MPI_Bcast(list, sizeof(Particle)*worldSize, MPI_BYTE, 0, MPI_COMM_WORLD);

	/*Print the updated positions*/

	cout << endl << endl;
	
	for(int i = 0; i<worldSize; ++i)
	{
		updatePosition(worldRank);
		MPI_Barrier(MPI_COMM_WORLD);
		if (worldRank == 0) {
			cout << endl << "Iteration " << i + 1 << endl;
			cout << "In time " << i*dt << endl;
			cout << "With displacement accuracy of " << dr <<endl;
			for (int j = 0; j < worldSize; j++) {
				cout << endl << "Particle " << j + 1 << endl;
				cout << "Velocity : " << "(" << list[j].velocity.x << " , " << list[j].velocity.y << ")" << endl;
				cout << "Mass : " << list[j].mass << endl;
				cout << "Charge : " << list[j].charge << endl;
				cout << "Position : " << "(" << list[j].position.x << " , " << list[j].position.y  << ")" << endl;
			}
		}
		/*Update position and velocity of current particle*/
		MPI_Allgather(&list[worldRank], sizeof(Particle), MPI_BYTE, &list[worldRank], sizeof(Particle), MPI_BYTE, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Finalize();
}
