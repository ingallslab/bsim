/**
 * BSimParticle.java
 *
 * Abstract class which defines the basic properties of an object in our simulation. 
 * It cannot be used directly and should instead be extended to include additional
 * properties of a real object.
 *
 * Authors: Thomas Gorochowski
 *          Ian Miles
 *          Mattia Fazzini(Update)
 * Created: 12/07/2008
 * Updated: 07/08/2009
 */
package bsim.particle;

import javax.vecmath.Vector3d;

import bsim.BSimParameters;
import bsim.particle.bacterium.BSimBacterium;
import bsim.particle.bead.BSimBead;
import bsim.particle.vesicle.BSimVesicle;
import bsim.scene.BSimScene;

public abstract class BSimParticle {		

	private Vector3d position; // microns		
	private Vector3d force; // piconewtons	
	private double radius; // microns	
	private BSimScene scene; // the environment that the particle exists in
		
	public BSimParticle(Vector3d newPosition, double newRadius, BSimScene newScene) {	
		position = new Vector3d();
		force = new Vector3d();		
		radius = newRadius;
		scene = newScene;
		
		position.set(newPosition);
	}	
	
	/*
	 * Actions independent of other obstacles: flagellar forces, adding chemicals, etc
	 */		
	public abstract void action();
	
	/*
	 * Update the position of the particle according to Stokes' law	
	 * Payable in force, yarr
	 */
	public void updatePosition() {
		Vector3d velocity = new Vector3d();
		velocity.scale(1/stokesCoefficient(), force); // pN/(micrometers*Pa sec) = micrometers/sec 
		position.scaleAdd(BSimParameters.dt, velocity, position);
		force.set(0,0,0);
		
		if(position.x > BSimParameters.xBound) position.x -= BSimParameters.xBound;
		if(position.x < 0) position.x += BSimParameters.xBound;
		if(position.y > BSimParameters.yBound) position.y -= BSimParameters.yBound;
		if(position.y < 0) position.y += BSimParameters.yBound;
		if(position.z > BSimParameters.zBound) position.z -= BSimParameters.zBound;
		if(position.z < 0) position.z += BSimParameters.zBound;
	}
		
	public Vector3d getPosition() { return position; }	
	public Vector3d getForce() { return force; }	
	public double getRadius() { return radius; }
	public double getSurfaceArea() { return 4*Math.PI*Math.pow(radius,2); }
	public BSimScene getScene() { return scene; }
	public void addForce(Vector3d f) { force.add(f); }
	public void setRadius(double r) { radius = r; }
	public double stokesCoefficient() { return 6.0*Math.PI*radius*BSimParameters.visc; } // micrometers*Pa sec
		
	public double distance(BSimParticle p) {
		Vector3d d = new Vector3d();
        d.sub(this.position, p.position);
        return d.length();
	}
	
	public double outerDistance(BSimParticle p) {
		return this.distance(p) - (this.radius + p.radius);
	}	
        	
    /*
     * Applies a force on this of magnitude m towards this,
     * and a force on p of magnitude m towards p.
     */
	public void reaction(BSimParticle p, double m) {
		Vector3d f = new Vector3d();
		f.sub(this.position, p.position);			
		f.normalize();
		f.scale(m);
		this.addForce(f);
		f.negate();
		p.addForce(f);
	}
		
}

