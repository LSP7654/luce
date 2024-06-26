 /* Body.h
 *
 *  Created on: Nov 8, 2012
 *      Author: dh4gan
 *
 *      This is the header file setting up the Body Class
 */


//using namespace std;

#ifndef BODY_H_
#define BODY_H_

#include <string>
#include "math.h"
#include "Vector3D.h"
#include "Constants.h"
#include "parFile.h"
class Body {
public:
	Body();
	Body(string &namestring, double &m, double &rad, Vector3D &pos, Vector3D &vel);

	Body(string &namestring, double &m, double &rad, double semimaj, double ecc, double inc,
		double longascend, double argper, double meananom, double G, double totalMass);

    Body(parFile &input, int &bodyIndex, double &G);
	virtual ~Body();

    
    Body* nBodyClone();
    
	/* Accessor methods */

	string getName() { return name; }
	double getMass() { return mass; }
	string getType() { return type; }
	double getRadius(){ return radius;}
	double getTimestep(){ return timestep;}
	bool willBounce(){ return collisionBounce;}

	Vector3D getPosition() { return position; }
	Vector3D getVelocity() { return velocity;}
	Vector3D getAcceleration() { return acceleration;}
	Vector3D getJerk() { return jerk;}
	Vector3D getSnap() { return snap;}
	Vector3D getCrackle() { return crackle;}

	double getSemiMajorAxis(){return semiMajorAxis;}
	double getPeriod(){return period;}

	Vector3D getOrbitalAngularMomentum(){ return orbitalAngularMomentum;}
	double getMagOrbitalAngularMomentum(){ return magOrbitalAngularMomentum;}

	Vector3D getEccentricityVector(){ return eccentricityVector;}
	double getEccentricity(){ return eccentricity;}

	double getInclination(){return inclination;}
	double getTrueAnomaly(){return trueAnomaly;}
	double getMeanAnomaly(){return meanAnomaly;}
	double getEccentricAnomaly(){return eccentricAnomaly;}
	double getArgumentPeriapsis(){return argumentPeriapsis;}
	double getLongitudePeriapsis(){return longitudePeriapsis;}
	double getLongitudeAscendingNode(){return longitudeAscendingNode;}

	/* Variable defining methods  */

	void setName(string newName) { name = newName; }
	void setMass(double m) { mass = m; }
	void setRadius(double r){radius=r;}
	void setTimestep(double dt){timestep=dt;}
	void setBounce(bool willBounce){collisionBounce=willBounce;}

	void setPosition(Vector3D pos) { position = pos; }
	void setVelocity(Vector3D vel) { velocity = vel; }
	void setAcceleration(Vector3D acc) { acceleration = acc; }
	void setJerk(Vector3D jk) { jerk = jk; }
	void setSnap(Vector3D snp) { snap = snp; }
	void setCrackle(Vector3D crack) { crackle = crack; }

	void setSemiMajorAxis(double a){semiMajorAxis=a;}

	void setOrbitalAngularMomentum(Vector3D angmom) { orbitalAngularMomentum = angmom;}
	void setMagOrbitalAngularMomentum(double magangmom) { magOrbitalAngularMomentum = magangmom;}

	void setEccentricityVector(Vector3D ecc){ eccentricityVector = ecc;}
	void setEccentricity(double magecc){ eccentricity = magecc;}

	void setInclination(double inc){inclination=inc;}
	void setTrueAnomaly(double anom){trueAnomaly=anom;}
	void setMeanAnomaly(double anom){meanAnomaly = anom;}
	void setArgumentPeriapsis(double argper){argumentPeriapsis=argper;}
	void setLongitudePeriapsis(double longper){longitudePeriapsis=longper;}
	void setLongitudeAscendingNode(double node){longitudeAscendingNode=node;}


	// Standard cloning method
	virtual Body* Clone() { return new Body(*this); }

	// Orbital Calculation Methods

	void calcOrbitalAngularMomentum();
	void calcEccentricity(double G, double totmass);
	void calcTrueAnomaly(double G, double totalMass, double time);

	void calcEccentricAnomaly();
	void calcTrueAnomaly();
	void calcOrbitFromVector(double G, double totmass); // Calculate orbital elements from position and velocity
	void calcVectorFromOrbit(double G, double totmass); // Calculate position and velocity from orbital elements
	double calcPeriod(double G, double totalMass); // This method is virtual as it is used by World
	void changeFrame(Vector3D framepos, Vector3D framevel);

	void calcOrbitFromVector(double G, Body* parentBody); // Calculate orbit around a specific Body

	// N Body Calculation Methods

	void calcTimestep(double greekEta); // Calculates the preferred timestep given Body's state vectors
	void calcAccelJerk(double G, vector<Body*> bodyarray, double softening_length); // Calculates gravitational acceleration and jerk (first derivative)
	void calcSnapCrackle(double G, vector<Body*> bodyarray, double softening_length); // Calculates next two derivatives of acceleration


    void setHostMass(double m){hostMass = m;}
    double getHostMass(){return hostMass;}
    
    void setHostBody(Body* bod){hostBody = bod;}
    Body* getHostBody(){return hostBody;}
    
	double safeAcos(double x) {
		if (x < -1.0)
			x = -1.0;
		else if (x > 1.0)
			x = 1.0;
		return acos(x);
	}

	/*
	 *  Virtual Methods for Derived Classes to call
	 *  NEVER CALL THESE ON BODY OBJECTS!
	 *
	 */

	// Methods for Star Class
	virtual void setLuminosity(double lum){};
	virtual double getLuminosity() {return -1.0;}
	virtual void calcMainSequenceLuminosity(){}
	virtual double calculatePeakWavelength(){return -1.0;}

	virtual double getTeff(){return -1.0;}

	// Methods for Planet Class
	virtual void setEquilibriumTemperature(double temp){};
	virtual double getEquilibriumTemperature(){return -1.0;}

	virtual void setReflectiveLuminosity(double lum){}

	virtual void calcLuminosity(){};

	virtual void setAlbedo(double alb){}
	virtual double getAlbedo(){return -1.0;}

	// Methods for PlanetSurface Class

	virtual int getNLongitude(){return -1;}
	virtual int getNLatitude(){return -1;}
	virtual int getNStars(){return -1;}
	virtual int getLongPick(){return -1;}
	virtual int getLatPick(){return -1;}

	virtual double getPSpin(){return -1.0;}
	virtual double getObliquity(){return -1.0;}
	virtual double getFluxMax(){return -1.0;}
	virtual double getFluxMin(){return -1.0;}
	virtual double getAvgFluxMax(){return -1.0;}
	virtual double getPMax(){return -1.0;}
	virtual double getPMin(){return -1.0;}

	virtual void setNLongitude(int nlong){}
	virtual void setNLatitude(int nlat){}
	virtual void setNStars(int s){}

	virtual void setPSpin(double spin){}
	virtual void setObliquity(double obliq){}

	virtual void initialiseArrays(){};
	virtual void initialiseOutputVariables(string prefixString, vector<Body*> stars){};

	virtual void resetFluxTotals(){};

	virtual void findSurfaceLocation(double &longitude, double &latitude){};
	virtual void calcLongitudeOfNoon(Body* &star, int &istar){};

	virtual void calcFlux(int &istar, Body* &star, double &eclipseFraction, double &time, double &dt){};
	virtual void calcIntegratedQuantities(double &dt){};

//Written 1/25/24 by LSP7654
	//Updated 3/12/24
	virtual void calcAverageFlux(int snapshotNumber, double &dt){};

	virtual void writeFluxFile(int &snapshotNumber, int &nTime, double &time, string prefixString){};
	virtual void writeSkyFile(FILE* outputSky, int &istar, double &time){};
	virtual void writeToLocationFiles(double &time, vector<Body*> bodies){};

	virtual void writeIntegratedFile(){};

	//Written 1/25/24 by LSP7654
		//Updated 3/12/24
	virtual void writeAverageFile(){};

	//Written 4/17/24 by LSP7654
	virtual void calcPhotosynthRate(double &dt){};
	virtual void writePrateFile(){};


	// Variables that are part of the Body Class and its derivations //
protected:

	// Basic variables

	string name;
	string type;
	double mass;
	double radius;
	double timestep;	// Body's own calculated timestep - used to define global timestep
	bool collisionBounce;  // Does the Body bounce or not when it hits another Body?

	// Primitive Spatial Vectors

	Vector3D position;
	Vector3D velocity;
	Vector3D acceleration;
	Vector3D jerk;   			// Jerk = time derivative of acceleration
	Vector3D snap;			// Derivative of Jerk
	Vector3D crackle;			// Second Derivative of Jerk

	// Orbital Elements

	double semiMajorAxis;
	double period;

	Vector3D eccentricityVector;  // Eccentricity Vector
	double eccentricity; 	  // Eccentricity = Absolute Magnitude of Eccentricity Vector

	Vector3D orbitalAngularMomentum;
	double magOrbitalAngularMomentum;

	double inclination;
	double trueAnomaly;
	double meanAnomaly;
	double eccentricAnomaly;
	double argumentPeriapsis;
	double longitudePeriapsis;
	double longitudeAscendingNode;
    
    Body* hostBody;
    double hostMass;

};


#endif /* BODY_H_ */
