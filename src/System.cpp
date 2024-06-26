/*
 * System.cpp
 *
 *  Created on: Nov 8 2012
 *      Author: dh4gan
 */

#include "System.h"
#include "Constants.h"
#include <iostream>
#include <stdio.h>
#include <algorithm>

System::System()
    {
    name = "System";
    bodies = vector<Body*> (0);
    bodyCount = 0;
    Vector3D zeroVector;

    totalMass = 0.0;
    initialEnergy = 0.0;
    totalEnergy = 0.0;

    timeStep = 0.0;
    timeControl = 0.00002;

    //G = 1.0;
    G = 6.67425e-11;
    softeningLength = 1.0e-5;

    initialAngularMomentum = zeroVector;
    totalAngularMomentum = zeroVector;

    deltaAngularMomentum = 0.0;
    deltaEnergy = 0.0;

    positionCOM = zeroVector;
    velocityCOM = zeroVector;
    accelerationCOM = zeroVector;

    planetaryIlluminationOn = false;


    }

System::System(string &namestring, vector<Body*> &bodyarray)
    {

    name = namestring;
    bodies = bodyarray;
    bodyCount = bodies.size();
    Vector3D zeroVector;
    totalMass = 0.0;
    for (int i = 0; i < bodyCount; i++)
	{
	totalMass = totalMass + bodies[i]->getMass();
	}

    initialEnergy = 0.0;
    totalEnergy = 0.0;

    timeStep = 0.0;
    timeControl = 0.00002;

    G = 1.0;
    softeningLength = 1.0e-5;

    initialAngularMomentum = zeroVector;
    totalAngularMomentum = zeroVector;

    deltaAngularMomentum = 0.0;
    deltaEnergy = 0.0;

    positionCOM = zeroVector;
    velocityCOM = zeroVector;
    accelerationCOM = zeroVector;

    planetaryIlluminationOn = false;
    }

System::~System()
    {

    for (int i = 0; i << bodies.size(); i++)
	{
	delete bodies[i];
	bodies[i] = 0;
	}

    }

double safeAcos(double x) {
	if (x < -1.0)
		x = -1.0;
	else if (x > 1.0)
		x = 1.0;
	return acos(x);
}

/* Methods for Controlling Population of Bodies */

void System::addBody(Body* &newBody)
    {/* Author: dh4gan
     Adds a Body object to the System object */

    bodies.push_back(newBody);
    bodyCount = bodies.size();

    totalMass = totalMass + newBody->getMass();

    }

void System::removeBody(int bodyindex)
    {
    /* Author: dh4gan
     * removes Body (bodyindex) from the system
     */

    bodies.erase(bodies.begin() + bodyindex);

    // As a body has been deleted, must recalculate initial values of Energy, Angular Momentum
    calcInitialProperties();

    }


void System::setHostBodies(vector<int> orbitCentre)
    {
    /*
     * Written 19/8/14 by dh4gan
     * Sets up Host Bodies for Worlds where appropriate
     * Catalogue how much mass is orbiting each body
     */

    for (int i=0; i<bodyCount; i++)
	{
	    if(orbitCentre[i]>0)
		{
		bodies[i]->setHostBody(bodies[orbitCentre[i]-1]);
		bodies[orbitCentre[i]-1]->setHostMass(bodies[orbitCentre[i]-1]->getHostMass() + bodies[i]->getMass());
		}
	}

    }

/* Calculation Methods */




void System::calcCOMFrame(vector<int> participants)
    {/* Author: dh4gan
     * Calculates the position and velocity of the centre of mass
     * for a limited range of participants
     */

    int i;
    double m, participantMass;
    Vector3D pos, vel, zerovector;

    positionCOM = zerovector;
    velocityCOM = zerovector;

    participantMass = 0.0;

    for (i = 0; i < bodyCount; i++)
	{

	if (participants[i] == 1)
	    {
	    pos = bodies[i]->getPosition();
	    vel = bodies[i]->getVelocity();
	    m = bodies[i]->getMass();
	    participantMass +=m;
	    positionCOM = positionCOM.addVector(pos.scaleVector(m));
	    velocityCOM = velocityCOM.addVector(vel.scaleVector(m));

	    }
	}

    if(participantMass>0.0)
    {
    participantMass = 1.0/participantMass;
    }

    positionCOM = positionCOM.scaleVector(participantMass);
    velocityCOM = velocityCOM.scaleVector(participantMass);

    }



void System::calcCOMFrame()
    {/* Author: dh4gan
     Calculates the position and velocity of the centre of mass
     for all bodies */

    vector<int> participants(bodyCount,1);

    calcCOMFrame(participants);

    }
void System::transformToCOMFrame(vector<int> participants)
    {
    /* Author: dh4gan
     * Calculates Centre of Mass Frame and transforms the system to it
     * Calls the System method calcCOMFrame first */

    int i;
    Vector3D pos, vel;

    // Calculate Centre of Mass Vectors

    calcCOMFrame(participants);

    // Subtract these from Body vectors

    for (i = 0; i < bodyCount; i++)
	{
	if(participants[i]==1)
	    {
	pos = bodies[i]->getPosition();
	vel = bodies[i]->getVelocity();

	pos = pos.subtractVector(positionCOM);
	vel = vel.subtractVector(velocityCOM);

	bodies[i]->setPosition(pos);
	bodies[i]->setVelocity(vel);
	    }

	}
    //end of module
    }

void System::transformToCOMFrame()
    {
    /* Author: dh4gan
     * Calculates Centre of Mass Frame and transforms the system to it
     * Calls the System method calcCOMFrame first */

    vector<int> participants(bodyCount,1);

    transformToCOMFrame(participants);

    }


void System::transformToBodyFrame(int bodyIndex)
    {
    /* Author: dh4gan
     * Transforms system to the frame where body bodyIndex is at rest*/


    int i;
    Vector3D pos, vel;

    // Calculate Frame's vectors
    Vector3D framepos = bodies[bodyIndex]->getPosition();
    Vector3D framevel = bodies[bodyIndex]->getVelocity();


    // Subtract these from Body vectors

    for (i = 0; i < bodyCount; i++)
	{
	pos = bodies[i]->getPosition();
	vel = bodies[i]->getVelocity();

	pos = pos.subtractVector(framepos);
	vel = vel.subtractVector(framevel);

	bodies[i]->setPosition(pos);
	bodies[i]->setVelocity(vel);

	}
    //end of module
    }


void System::calcHostCOMFrame(Body* host, Vector3D &hostCOM, Vector3D &hostvelCOM)
{
	/*
	 * Written 8/1/14 by dh4gan
	 * Transforms system so that COM of system belonging to body host is at centre
	 */

	vector<int> participants(bodyCount,0);
	double m, participantMass;
	Vector3D pos, vel, zerovector;


	// Find all bodies with this host
	for(int i=0; i<bodyCount; i++)
	{
		if(bodies[i]==host)
		{
			participants[i]=1;
		}
		if(bodies[i]->getHostBody()==host)
		{
			participants[i]=1;
		}

	}

	// Calculate COM frame of these bodies

	hostCOM = zerovector;
	hostvelCOM = zerovector;
	participantMass = 0.0;
	    for (int i = 0; i < bodyCount; i++)
		{

		if (participants[i] == 1)
		    {
		    pos = bodies[i]->getPosition();
		    vel = bodies[i]->getVelocity();
		    m = bodies[i]->getMass();
		    participantMass +=m;
		    hostCOM = hostCOM.addVector(pos.scaleVector(m));
		    hostvelCOM = hostvelCOM.addVector(vel.scaleVector(m));

		    }
		}

	    if(participantMass>0.0)
	    {
	    participantMass = 1.0/participantMass;
	    }

	    hostCOM = hostCOM.scaleVector(participantMass);
	    hostvelCOM = hostvelCOM.scaleVector(participantMass);

	    }

void System::transformToArbitraryFrame(Vector3D framePosition, Vector3D frameVelocity, vector<int> participants)
{
	/*
	 * Written 8/1/14 by dh4gan
	 * Transforms part of the system to arbitrary reference frame
	 * Bodies taking part in the transformation have non-zero entries in participant array
	 */

	// Find all bodies with this host
	for(int i=0; i<bodyCount; i++)
	{
		if(participants[i]==1)
		{
		bodies[i]->setPosition(bodies[i]->getPosition().subtractVector(framePosition));
		bodies[i]->setVelocity(bodies[i]->getVelocity().subtractVector(frameVelocity));
		}
	}

}

void System::transformToArbitraryFrame(Vector3D framePosition, Vector3D frameVelocity)
{
	/*
	 * Written 8/1/14 by dh4gan
	 * Transforms entire system to arbitrary reference frame
	 */
	vector<int> participants(bodyCount,1);

	transformToArbitraryFrame(framePosition,frameVelocity,participants);

}




void System::calcTotalEnergy()
    {
    /* Author: David the legend Harvey
     *
     * Purpose :
     * 		Calculate the total energy of the system
     *
     * Method :
     *      The total energy is the kinetic plus the gravitational potential
     *      so will need to loop through each body and work out the total of
     *      each
     *
     */

    int body;
    int other_body;
    int number_bodies;

    double r_distance;
    double body_ref_mass;
    double body_question_mass;
    double gravitational_potential = 0.0;
    double kinetic_energy = 0.0;
    double velocity;

    Vector3D r_vector;
    Vector3D body_ref_vel;
    Vector3D body_ref_pos;
    Vector3D body_question_pos;
    Vector3D r_hat;
    Vector3D force_total;

    //Should Test this as this could be faulty
    number_bodies = bodies.size();

    for (body = 0; body < number_bodies; body++)
	{

	//This is the body that is the reference for calc
	//the velocity
	body_ref_pos = bodies[body]->getPosition();
	body_ref_mass = bodies[body]->getMass();
	body_ref_vel = bodies[body]->getVelocity();

	for (other_body = 0; other_body < number_bodies; other_body++)
	    {
	    if (body != other_body)
		{
		//Need to loop through the other bodies in order to
		//find the gravitational potential

		body_question_pos = bodies[other_body]->getPosition();
		body_question_mass = bodies[other_body]->getMass();

		r_vector = body_question_pos.relativeVector(body_ref_pos);

		r_distance = r_vector.magVector();

		gravitational_potential = (-G * body_ref_mass
			* body_question_mass) / r_distance;

        //  cout << "body_ref_mass =" << body_ref_mass << endl;
        //  cout << "body_question_mass =" << body_question_mass << endl;
        //  cout << "r_distance =" << r_distance << endl;
        //  cout << "GP =" << gravitational_potential << endl;

		}
	    }
	velocity = body_ref_vel.magVector();
	kinetic_energy += 0.5 * body_ref_mass * velocity * velocity;
	}

    totalEnergy = kinetic_energy + gravitational_potential;

    deltaEnergy = (totalEnergy - initialEnergy) / initialEnergy;

    }

void System::calcTotalAngularMomentum()
    {
    /*Author: dh4gan
     * Calculates the total Angular Momentum in the System, where the rotation axis is at the origin
     * Also calculates the deviation from the initial value.
     */

    Vector3D angularMomentum;
    Vector3D pos, vel, zerovector;

    double magTotalAngularMomentum, magInitialAngularMomentum;

    totalAngularMomentum = zerovector;

    for (int i = 0; i < bodyCount; i++)
	{
	pos = bodies[i]->getPosition().scaleVector(bodies[i]->getMass()); // multiply r by mass here
	vel = bodies[i]->getVelocity();

	totalAngularMomentum = totalAngularMomentum.addVector(pos.crossProduct(
		vel));

	}

    magTotalAngularMomentum = totalAngularMomentum.magVector();
    magInitialAngularMomentum = initialAngularMomentum.magVector();

    if (magInitialAngularMomentum != 0.0)
	{
	deltaAngularMomentum = (magInitialAngularMomentum
		- magTotalAngularMomentum) / magInitialAngularMomentum;
	}

    }

void System::calcNBodyTimestep(vector<Body*> &bodyarray, double dtmax)
    {/* Author: dh4gan
     Calls the calcTimestep method for every Body object in the System object,
     and finds the minimum value */

    int i;
    vector<double> dt(bodyCount,0.0);
    double dtarraymax;

#pragma omp parallel default(none) \
	shared(dt) \
	private(i)
	{
	for (i = 0; i < bodyCount; i++)
	    {

	    bodyarray[i]->calcTimestep(timeControl);

	    dt[i] = bodyarray[i]->getTimestep();
	    }

	}
    dtarraymax = *(min_element(dt.begin(), dt.end()));

    if(dtarraymax>dtmax)
	{
	dtarraymax = dtmax;
	}
    timeStep = dtarraymax;

    }


void System::calcNBodyTimestep(double dtmax)
{/* Author: dh4gan
  overloaded Timestep method - defaults to the object's own body array
  the calcTimestep method for every Body object in the System object,
  and finds the minimum value */
    
    calcNBodyTimestep(bodies,dtmax);
    
}

void System::setupOrbits(vector<int> orbitCentre)
    {

    /*
     * Written 22/1/14 by dh4gan
     * Given a list of orbit centres (i=0: CoM; i>0: body, i<0: (0,0,0))
     *
     */

    vector <int> participants(bodyCount,0);
    vector <int> hosts(bodyCount, 0);

    Vector3D hostPosition,hostVelocity;
    Vector3D hostCOM,hostvelCOM;

    // Firstly, set up desired objects around centre of mass

    for (int b = 0; b < bodyCount; b++)
	{
	if (orbitCentre[b] == 0)
	    {
	    participants[b] = 1;

	    // Set up their orbits
	    bodies[b]->calcVectorFromOrbit(G,totalMass);
	    }
	if(orbitCentre[b]<0)
	    {

	    bodies[b]->calcVectorFromOrbit(G,totalMass-bodies[b]->getMass());
	    }


	}


    // Transform them to COM Frame
    transformToCOMFrame(participants);


    // Set up bodies with specific hosts
    // i) - place satellites in orbits assuming host is at the centre
    // ii) - shift host so that CoM of host-satellite system follows original orbit specified for host


    // r_host new =  1/m_host* (r_host*mtot - (sum_i(i/=host) m_i r_i))
    // similar for v_host

    // r_host new = r_host as mass of satellites goes to zero


    // Define arrays to compute total mass of host systems and CoM of satellites
    vector<double> totalMassAroundHost(bodyCount,0.0);
    vector<Vector3D> xCOMroundHost(bodyCount);
    vector<Vector3D> vCOMroundHost(bodyCount);

    // i) Place satellites around host at centre

    for (int b = 0; b < bodyCount; b++)
	{
        
        // Identify hosts and compute total mass belonging to each host
	if (orbitCentre[b] > 0)
	    {

	    int ihost = orbitCentre[b]-1;
	    hosts[ihost]=1;  // Record host status for later

		totalMassAroundHost[ihost] = totalMassAroundHost[ihost]+bodies[b]->getMass();

	    }
	}


        // Add host mass to totals
	for(int b=0; b< bodyCount; b++)

	{
	if(hosts[b]==1)

	{totalMassAroundHost[b] = totalMassAroundHost[b] + bodies[b]->getMass();}
	}


    // Setup bodies in orbit around each host
   for (int b=0; b<bodyCount; b++)

	{

	if(orbitCentre[b]>0)

	{
		int ihost = orbitCentre[b]-1;

	    bodies[b]->calcVectorFromOrbit(G,
		    totalMassAroundHost[ihost]);

	    Vector3D framepos =
		    bodies[ihost]->getPosition().scaleVector(-1.0);
	    Vector3D framevel =
		    bodies[ihost]->getVelocity().scaleVector(-1.0);
	    bodies[b]->changeFrame(framepos, framevel);

	    // Compute COM of host system (minus host contribution)
	    xCOMroundHost[ihost] = xCOMroundHost[ihost].addVector(bodies[b]->getPosition().scaleVector(bodies[b]->getMass()));
	    vCOMroundHost[ihost] = vCOMroundHost[ihost].addVector(bodies[b]->getVelocity().scaleVector(bodies[b]->getMass()));

	    }

	}

    // Now ensure that each host system is setup such that the CoM of the system
    // moves along the host's assigned orbit

    // Do this by moving the host and altering its velocity as specified in above comments

    for (int b=0; b<bodyCount; b++)

      {

	// If this body is a host, then alter its position and velocity
	if (hosts[b]==1)
	  {

	    Vector3D newPosition = bodies[b]->getPosition().scaleVector(totalMassAroundHost[b]).subtractVector(xCOMroundHost[b]);
	    newPosition = newPosition.scaleVector(1.0/bodies[b]->getMass());

	    Vector3D newVelocity = bodies[b]->getVelocity().scaleVector(totalMassAroundHost[b]).subtractVector(vCOMroundHost[b]);
	    	    newVelocity = newVelocity.scaleVector(1.0/bodies[b]->getMass());

	    bodies[b]->setPosition(newPosition);
	    bodies[b]->setVelocity(newVelocity);

	  }
      }


}


void System::calcForces(vector<Body*> &bodyarray)
    {
    /* Author: dh4gan
     This is a wrapper method, which abstracts the process of
     calculating accelerations, jerks, snaps and crackles of an array of bodies
     */

    int i;
    int length = bodyarray.size();
    Vector3D zeroVector;

#pragma omp parallel default(none) \
	shared(length,zeroVector) \
	private(i)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < length; i++)
	    {
	    bodyarray[i]->setAcceleration(zeroVector);
	    bodyarray[i]->setJerk(zeroVector);
	    bodyarray[i]->calcAccelJerk(G, bodyarray, softeningLength);

	    }
	}
#pragma omp parallel default(none) \
	shared(length,zeroVector) \
	private(i)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < length; i++)
	    {
	    bodyarray[i]->setSnap(zeroVector);
	    bodyarray[i]->setCrackle(zeroVector);
	    bodyarray[i]->calcSnapCrackle(G, bodyarray, softeningLength);
	    }
	}

    }

void System::calcInitialProperties()
    {
    /* Author: dh4gan
     * Calculates initial energy and angular momentum
     * and stores it for later evaluation
     *
     */

    // Update total number of bodies
    bodyCount = bodies.size();

    // Total Mass
    totalMass = 0.0;

    for (int i = 0; i < bodyCount; i++)
	{
	totalMass += bodies[i]->getMass();
	}

    // Put the system in the Centre of Mass Frame
    transformToCOMFrame();

    // Calculate initial forces
    calcForces(bodies);


    // Calculate Total Energy, and define initial value

    initialEnergy = totalEnergy;
    deltaEnergy = 0.0;


    //Calculate Total Angular Momentum, and define initial value
    calcTotalAngularMomentum();
    initialAngularMomentum = totalAngularMomentum;
    deltaAngularMomentum = 0.0;


    }

vector<double> System::checkForEclipses(int bodyIndex)
{
    /* Author: dh4gan
     This method looks from the position of the body given by bodyindex
     and checks whether any other bodies are obscured from view
     It does this by looking at impact parameters from the relative
     vector of a given body and bodyindex.  Any body which has an
     impact parameter smaller than its radius, and is in front of
     the given body, obscures it.

     The method returns a vector of doubles, describing what fraction of the star is eclipsed [0,1]
     */

	vector<double> eclipsefrac(bodyCount,0.0);
	Vector3D vector_i, vector_j;
	double mag_i, mag_j, idotj, b;
	double rad_i, rad_j, rad_i2, rad_j2;
	double angle1,angle2, area_i, area_j;

	for(int i=0; i<bodyCount; i++)
	{

		// Body cannot eclipse itself
		if(i==bodyIndex)
		{
			eclipsefrac[i]=0.0;
			continue;
		}

		//calculate relative vector between i and bodyindex

		vector_i = getBody(bodyIndex)->getPosition().relativeVector(getBody(i)->getPosition());
		mag_i = vector_i.magVector();

		// Get Radius of body i

		rad_i = getBody(i)->getRadius();

		// Now loop over other bodies to get impact parameters

		for (int j=0; j <bodyCount; j++)
		{
			// Skip for bodies i and bodyIndex

			if(j==i or j==bodyIndex)
			{
				continue;
			}

			// Calculate relative vector between body j and bodyIndex

			vector_j = getBody(bodyIndex)->getPosition().relativeVector(getBody(j)->getPosition());
			mag_j = vector_j.magVector();

			// Get Radius of body j

			rad_j = getBody(j)->getRadius();

			// impact parameter = mag(vector_j)sin alpha = mag(vector_j)*sqrt(1-(vector_i.vector_j)^2)

			if(mag_i*mag_j >0.0 and mag_i>mag_j)
			{
				idotj = vector_i.dotProduct(vector_j)/(mag_i*mag_j);
				b = mag_j*sqrt(1.0-idotj*idotj);


				// If impact parameter less than radius, and i further away than j, eclipse of i!
				// Calculate area covered during eclipse
				// (sum of two circular segment areas, one for each circle)

				if(b < (rad_i+rad_j) and idotj < 0.0)
				{

				    //b = b/rad_i;
				    //rad_i = 1.0;
				    //rad_j = rad_j/rad_i;


					rad_i2 = rad_i*rad_i;
					rad_j2 = rad_j*rad_j;
					angle1 = 2.0*safeAcos((rad_i2 + b*b -rad_j2)/(2.0*b*rad_i));
					angle2 = 2.0*safeAcos((rad_j2 + b*b -rad_i2)/(2.0*b*rad_j));

					area_i = 0.5*rad_i2*(angle1 - sin(angle1));
					area_j = 0.5*rad_j2*(angle2 - sin(angle2));

					eclipsefrac[i] = (area_i+area_j)/(pi*rad_i2);


					if(eclipsefrac[i]>1.0) {eclipsefrac[i] = 1.0;}
					if(eclipsefrac[i]<0.0) {eclipsefrac[i] = 0.0;}

				}
			}
			else
			{
				continue;
			}

		}// End loop over j to get impact parameters

	} //End loop over i to test eclipses


	return eclipsefrac;

}


void System::evolveSystem(double tbegin, double tend)
    {
    /* Author: dh4gan
     This method evolves the System of Body Objects using a 4th order Hermite integrator
     This is a predictor-corrector algorithm:
     the method stores predicted data in a vector of body objects
     */

    int i;
    double time;
    double dtmax = (tend-tbegin)/2.0;

    Vector3D pos, vel, acc, jerk, snap, crack; // Holders for the body state vectors
    Vector3D pos_p, vel_p, acc_p, jerk_p, snap_p, crack_p; // Holders for the predicted body state vectors
    Vector3D pos_c, vel_c, acc_c, jerk_c, snap_c, crack_c; // Holders for the predicted body state vectors

    Vector3D velterm, accterm, jerkterm;
    vector<Body*> predicted; // Vector to store the predicted body properties

    /* i. Calculate initial accelerations, jerks, snaps and crackles for all particles */

    calcInitialProperties();

    calcForces(bodies);

    /* ii. Calculate initial global timestep, total energy and total angular momentum */

    calcNBodyTimestep(bodies, dtmax);
    calcTotalEnergy();
    calcTotalAngularMomentum();

    /* iii Set predicted body array equal to the current body array */
    for (i=0;i<bodyCount; i++){
	predicted.push_back(bodies[i]->nBodyClone());
    }

    /* Begin loop over time */
    time = tbegin;

    while (time < tend)
	{


#pragma omp parallel default(none) \
	shared(predicted)\
	private(i,pos,vel,acc,jerk,pos_p,vel_p)
	{
#pragma omp for schedule(runtime) ordered
	/* Calculate predicted positions and velocities */
	for (i = 0; i < bodyCount; i++)
	    {

	    // Pull the body object's data //
	    pos = bodies[i]->getPosition();
	    vel = bodies[i]->getVelocity();
	    acc = bodies[i]->getAcceleration();
	    jerk = bodies[i]->getJerk();


	    // 1. Calculate predicted position and velocity //
	    pos_p = pos.addVector(vel.scaleVector(timeStep), acc.scaleVector(
		    0.5 * timeStep * timeStep), jerk.scaleVector(timeStep
		    * timeStep * timeStep / 6.0));

	    vel_p = vel.addVector(acc.scaleVector(timeStep), jerk.scaleVector(
		    0.5 * timeStep * timeStep));

	    // Update the object holding the predicted data
	    predicted[i]->setPosition(pos_p);
	    predicted[i]->setVelocity(vel_p);

	    }
	}

	/* 2. Use predicted positions and velocities to calculate
	 * predicted accelerations, jerks, snaps and crackles */

	calcForces(predicted);



#pragma omp parallel default(none) \
	shared(predicted)\
	private(i,pos,vel,acc,jerk)\
	private(pos_p,vel_p,acc_p,jerk_p)\
	private(velterm,accterm,jerkterm,vel_c,pos_c)
	{
#pragma omp for schedule(runtime) ordered
	for (i = 0; i < bodyCount; i++)
	    {

	    pos_p = predicted[i]->getPosition();
	    vel_p = predicted[i]->getVelocity();
	    acc_p = predicted[i]->getAcceleration();
	    jerk_p = predicted[i]->getJerk();

	    pos = bodies[i]->getPosition();
	    vel = bodies[i]->getVelocity();
	    acc = bodies[i]->getAcceleration();
	    jerk = bodies[i]->getJerk();

	    accterm = acc_p.addVector(acc).scaleVector(0.5 * timeStep);

	    jerkterm = jerk_p.relativeVector(jerk).scaleVector(timeStep
		    * timeStep / 12.0);

	    vel_c = vel.addVector(accterm, jerkterm);

	    accterm = acc_p.relativeVector(acc).scaleVector(timeStep * timeStep
		    / 12.0);
	    velterm = vel_c.addVector(vel).scaleVector(0.5 * timeStep);
	    pos_c = pos.addVector(velterm, accterm);

	    // update appropriate object with these new vectors

	    bodies[i]->setPosition(pos_c);
	    bodies[i]->setVelocity(vel_c);



	    }
	}

	/* 5. Calculate acceleration, jerk, snap and crackle for next step */

	calcForces(bodies);

	/* 5. Increase time by dt */
	time = time + timeStep;

	/* 6. Calculate new timestep, total energy, angular momentum and orbital data */
	calcNBodyTimestep(bodies, dtmax);
	calcTotalEnergy();
	calcTotalAngularMomentum();


	}
    // End of loop over time

    // Garbage collection on predicted pointers
    for(i=0; i< bodyCount; i++)
	{
	delete predicted[i];
	predicted[i]=0;
	}

    }

void System::evolveSystem(double dt)
    {

    double tbegin = 0.0;
    double tstop = dt;
    evolveSystem(tbegin,tstop);

    }

void System::calcPlanetaryEquilibriumTemperatures()
    {
    /*
     * Written 11/8/14 by dh4gan
     * Takes all Planets in the System, and calculates their equilibrium temperature
     * given the Stars in the System
     * Also calculates reflected starlight and thermal luminosity
     */

    double sep, temp, lum, rad, albedo;

    for (int j=0; j< bodyCount; j++)
	{
	if(bodies[j]->getType()=="Planet")
	    {

	    rad = bodies[j]->getRadius()*rsol/AU;
	    albedo = bodies[j]->getAlbedo();
	    temp = 0.0;
	    lum = 0.0;

	    for (int i=0; i<bodyCount; i++)
		{
		if(bodies[i]->getType() =="Star" and i!=j)
		    {

		    // Calculate separation between each star and planet
		    sep = bodies[j]->getPosition().subtractVector(bodies[i]->getPosition()).magVector();

		    // Add contribution to equilibrium temperature
		    temp = temp+ bodies[i]->getLuminosity()*lsol*(1.0-albedo)/(16.0*pi*sigma_SB*sep*sep*AU*AU);

		    //Calculate reflected starlight
		    lum = lum + bodies[i]->getLuminosity()*rad*rad*albedo/(sep*sep);

		    }
		}

	    temp = pow(temp,0.25);

	    // Set Planet's Equilibrium temperature
	    bodies[j]->setEquilibriumTemperature(temp);

	    // Set Planet's Reflective Luminosity
	    bodies[j]->setReflectiveLuminosity(lum);

	    // Calculate Planet's total luminosity
	    bodies[j]->calcLuminosity();

	    }

	}



    }

void System::calcLongitudesOfNoon()
{
    for (int j = 0; j < bodyCount; j++)
    {
        
        // If body is a PlanetSurface Object, loop through other bodies
        // and calculate the flux they emit onto its surface
        
        if (bodies[j]->getType() == "PlanetSurface")
        {
            
            for (int i = 0; i < bodyCount; i++)
            {
                
                if (bodies[i]->getType() == "Star")
                {
                    bodies[j]->calcLongitudeOfNoon(bodies[i], i);
                }
                
                if (bodies[i]->getType() == "Planet"
                    and planetaryIlluminationOn)
                {
                    bodies[j]->calcLongitudeOfNoon(bodies[i], i);
                }
            }
            
        }
    }
    
    
}

void System::calc2DFlux(int snapshotNumber, double &time, double &dt)
{
    /*
     * Written 11/12/14 by dh4gan
     * Method allows PlanetSurface Objects to call their flux calculation routines
     *
     */
    
    vector<double> eclipsefrac;
    
    for (int j=0; j<bodyCount; j++)
    {
        
        // If body is a PlanetSurface Object, loop through other bodies
        // and calculate the flux they emit onto its surface
        
        if(bodies[j]->getType()=="PlanetSurface")
        {
            
            bodies[j]->resetFluxTotals();
            eclipsefrac = checkForEclipses(j);
            for(int i=0; i< bodyCount; i++)
            {
                
                if(i==j) continue; // PlanetSurface can't emit flux onto itself
                
                
                if(bodies[i]->getType()=="Star")
                {
                    bodies[j]->calcFlux(i, bodies[i], eclipsefrac[i], time, dt);
                }
                
                if(bodies[i]->getType()=="Planet" && planetaryIlluminationOn)
                {
                    bodies[j]->calcFlux(i, bodies[i], eclipsefrac[i], time, dt);
                }
                
            }
            
            bodies[j]->calcIntegratedQuantities(dt);
            bodies[j]->calcAverageFlux(snapshotNumber, dt);
    
        }
        
    }
    
    
}

void System::calcPhotosynthRate(double &dt)
{
    //Written 4/17/24 by LSP7654
    //Method allows PlanetSurface Objects to call their photosythesis calculation routine

    for (int j=0; j<bodyCount; j++)
    {
        if(bodies[j]->getType()=="PlanetSurface")   
        {
            bodies[j]->calcPhotosynthRate(dt);
        }
    }

}

void System::initialise2DFluxOutput(string prefixString)
{
    /*
     * Written 12/12/14 by dh4gan
     * Initialises files for all PlanetSurface objects in the system
     *
     */
    
    for (int j = 0; j < bodyCount; j++)
    {
        
        if (bodies[j]->getType() == "PlanetSurface")
        {
            bodies[j]->initialiseOutputVariables(prefixString, bodies);
        }
        
    }
    
    calcLongitudesOfNoon();
    
}

void System::outputNBodyData(FILE *outputfile, double &time, vector<int> orbitCentre)
{
    /*
     * Written 10/1/14 by dh4gan
     * Method writes N Body information to
     * already open file pointer
     * orbitCentre vector determines where the orbits are calculated from
     */
    
    Vector3D position, velocity;
    // Transform to the Centre of Mass Frame
    transformToCOMFrame();
    //transformToBodyFrame(0);
    
    for (int j = 0; j < bodyCount; j++)
    {
        
        if(orbitCentre[j]>0)
        {
            bodies[j]->calcOrbitFromVector(G, bodies[orbitCentre[j]-1]);
        }
        
        else
        {
            bodies[j]->calcOrbitFromVector(G, totalMass);
        }
        
        position = bodies[j]->getPosition();
        velocity = bodies[j]->getVelocity();
        
        //Write out the data
        //Output format  CSV
        // mass,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z
        
        fprintf(outputfile,
                "%+.4E,%+.4E, %s,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,"
                "%+.4E,%+.4E,%+.4E,%+.4E,%+.4E,%+.4E\n", time,
                totalEnergy, bodies[j]->getName().c_str(), bodies[j]->getMass(),
                bodies[j]->getRadius(), position.elements[0],
                position.elements[1], position.elements[2],
                velocity.elements[0], velocity.elements[1],
                velocity.elements[2], bodies[j]->getSemiMajorAxis(),
                bodies[j]->getEccentricity(), bodies[j]->getInclination(),
                bodies[j]->getLongitudeAscendingNode(),
                bodies[j]->getArgumentPeriapsis(), bodies[j]->getMeanAnomaly());
        
    }
    fflush(outputfile);
    
}

void System::output2DFluxData(int &snapshotNumber, double &tSnap, string prefixString)
{
    
    /*
     * Written 10/1/14 by dh4gan
     * Requires Worlds to write their LEBM data to files
     *
     */
    
    for (int b=0; b < bodyCount; b++)
    {
        if(bodies[b]->getType()=="PlanetSurface")
        {
            
            if(FullOutput){
                bodies[b]->writeFluxFile(snapshotNumber, nTime, tSnap, prefixString);
                //bodies[b]->writeAverageFile(snapshotNumber, nTime, tSnap, prefixString);
            }
            
            bodies[b]->writeToLocationFiles(tSnap, bodies);
            
        }
        
    }
    
}

void System::outputSummaryFluxData() {
    /*
     * Written 4/12/14 by dh4gan
     * Updated 3/12/24 by LSP7654
     * Writes the integrated and average flux data to files
     *
     */
    
    for (int b = 0; b < bodyCount; b++) {
        if (bodies[b]->getType() == "PlanetSurface") {
            
            bodies[b]->writeIntegratedFile();
            bodies[b]->writeAverageFile();
                    
        }
        
    }
    return;   
}

void System::outputPrateData()
{
    //Written 4/17/14 by LSP7654
    //writes photosynthesis rates to file
    for (int b = 0; b < bodyCount; b++) {
        if (bodies[b]->getType() == "PlanetSurface") {
            bodies[b]->writePrateFile();
        }
    }
}

void System::outputInfoFile(int nSnaps)
{
    
    /*
     * Written 17/12/14 by dh4gan
     * Writes an info file for use in plotting 2D flux data
     *
     */
    
    string fileString = getName()+".info";
    infoFile = fopen(fileString.c_str(), "w");
    if (NULL == infoFile) 
    {
        cout << "Error opening File" << infoFile << endl;
        return;
    }

    
    double globalFluxMax= 0.0;
    double globalFluxMin= 1.0e50;
    double globalAvgFluxMax= 0.0;
    int nStars = countStars();

    double globalPMax = 0.0;
    double globalPMin = 1.0e50;
    
    fprintf(infoFile,"%i \n", nSnaps);
    fprintf(infoFile,"%i \n", nStars);
    
    for (int s = 0; s < bodyCount; s++)
    {
        if (bodies[s]->getType() == "Star")
        {
            fprintf(infoFile, "%s \n", bodies[s]->getName().c_str());
            fprintf(infoFile, "%+.4E %+.4E %+.4E \n", bodies[s]->getRadius(),
                    bodies[s]->getTeff(),
                    bodies[s]->calculatePeakWavelength());
        }
        
        if (bodies[s]->getType() == "PlanetSurface")
        {
            if (bodies[s]->getFluxMax() > globalFluxMax)
            {
                globalFluxMax = bodies[s]->getFluxMax();
            }
            
            if (bodies[s]->getFluxMin() < globalFluxMin)
            {
                globalFluxMin = bodies[s]->getFluxMin();
            }

            if (bodies[s]->getAvgFluxMax() > globalAvgFluxMax)
            {
                globalAvgFluxMax = bodies[s]->getAvgFluxMax();
            }

            if (bodies[s]->getPMax() > globalPMax)
            {
                globalPMax = bodies[s]->getPMax();
            }
            
            if (bodies[s]->getPMin() < globalPMin)
            {
                globalPMin = bodies[s]->getPMin();
            }
            
        }
    }
    
    fprintf(infoFile, "%+.4E \n", globalFluxMax);
    fprintf(infoFile, "%+.4E \n", globalFluxMin);
    fprintf(infoFile, "%+.4E \n", globalAvgFluxMax);
    fprintf(infoFile, "%+.4E \n", globalPMax);
    fprintf(infoFile, "%+.4E \n", globalPMin);
    
    printf("globalFluxMax: %+.4E \n", globalFluxMax);
    printf("globalFluxMin: %+.4E \n", globalFluxMin);
    // fprintf(infoFile, "%+.4E \n", globalPMax);
    // fprintf(infoFile, "%+.4E \n", globalPMin);
    printf("globalPMax: %+.4E \n", globalPMax);
    printf("globalPMin: %+.4E \n", globalPMin);

    fclose(infoFile);
    
    
}

int System::countStars()
{
    /*
     * Written 11/12/14 by dh4gan
     * Returns the number of Body objects with Type Star
     */
    
    int nStars = 0;
    for (int j = 0; j < bodyCount; j++)
    {
        if (bodies[j]->getType() == "Star")
        nStars++;
    }
    return nStars;
}

int System::countPlanets()
{
    /*
     * Written 11/12/14 by dh4gan
     * Returns the number of Body objects with Type Star
     */
    
    int nPlanets = 0;
    for (int j = 0; j < bodyCount; j++)
    {
        if (bodies[j]->getType() == "Planet")
        nPlanets++;
    }
    return nPlanets;
}

int System::countPlanetSurfaces()
{
    /*
     * Written 11/12/14 by dh4gan
     * Returns the number of Body objects with Type Star
     */
    
    int nPlanetSurfaces = 0;
    for (int j = 0; j < bodyCount; j++)
    {
        if (bodies[j]->getType() == "PlanetSurface")
        nPlanetSurfaces++;
    }
    return nPlanetSurfaces;
}

