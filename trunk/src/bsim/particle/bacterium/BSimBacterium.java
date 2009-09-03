/**
 * BSimBacterium.java
 *
 * Class that represents a bacterium in our simulation. This class implements attraction to
 * the goal chemoattractant. This class should be sub-classed to implement other GRNs,
 * such as the recruitment, or time delays, etc.
 *
 * Authors: Thomas Gorochowski
 *          Ian Miles
 *          Charlie Harrison
 *          Mattia Fazzini(Update)
 * Created: 12/07/2008
 * Updated: 07/08/2009
 */
package bsim.particle.bacterium;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;
import java.util.Vector;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Vector3d;

import bsim.BSimParameters;
import bsim.BSimUtils;
import bsim.particle.BSimParticle;
import bsim.particle.vesicle.BSimVesicle;
import bsim.scene.BSimScene;


public class BSimBacterium extends BSimParticle {
	
	// Motion states
	private static int RUNNING  = 1;
	private static int TUMBLING = 2;	
	
	// Run probabilities
	private static double runProbUp = 1 - BSimParameters.dt/BSimParameters.runLengthUp;
	private static double runProbDown = 1 - BSimParameters.dt/BSimParameters.runLengthDown;
	private static double runProbIso = 1 - BSimParameters.dt/BSimParameters.runLengthIso;
	
	// Values for gamma tumbling distribution
	private static double[] gammaVals = readGammaVals();	
	
	private Vector3d direction;
	private int motionState; // RUNNING or TUMBLING
	
	private double replicationRadius;
	private double radiusGrowthRate; // microns/sec
	
	private Vector<BSimVesicle> vesicles; // surface vesicles
	private double pNewVesicleDt; // probability of generating a new surface vesicle over dt
		
	private double shortTermMemoryDuration; // seconds
	private double longTermMemoryDuration; // seconds
	private double sensitivity; // to differences in long term vs short term mean concentrations
	private Vector<Double> memory; // memory of the concentration of the goal field
	
	// Set at onset of tumbling phase:
	private int tumbleSteps; 	// Number of time steps remaining in tumble phase
	private double tumbleAngle; // Angle remaining in tumble phase		
					
	/**
	 * General constructor.
	 */
	public BSimBacterium(Vector3d newPosition, double newRadius, Vector3d newDirection, BSimScene newScene) {
		super(newPosition, newRadius, newScene);		
		direction = newDirection;			
		motionState = RUNNING; // Start off running
		
		replicationRadius = 2*newRadius;
		radiusGrowthRate = 0.1;
				
		vesicles = new Vector<BSimVesicle>();
		pNewVesicleDt = 0.0001;
		
		shortTermMemoryDuration = 1.0;
		longTermMemoryDuration = 3.0; 
		sensitivity = 0.000001;			
		double memorySize = (shortTermMemoryDuration + longTermMemoryDuration) / BSimParameters.dt;
		memory = new Vector((int)(memorySize));		
	}
	
	public void action() {			
		if(motionState == RUNNING) {
			memory.remove(0); // forget the oldest concentration
			memory.add(getScene().getGoalField().getConcentration(this.getPosition())); // remember the newest concentration			
			run();			
		}
		else if(motionState == TUMBLING) {
			tumble();
		}
		
		if(Math.random() < pNewVesicleDt)
			vesicles.add(new BSimVesicle(0d));
		
		grow();		
	}
	
	protected void grow() {		
		double surfaceAreaConsumedByVesicles = 0;
		for(BSimVesicle vesicle : vesicles) {
			surfaceAreaConsumedByVesicles += vesicle.grow(this);			
		}			
		
		setRadius(getRadius() + radiusGrowthRate*BSimParameters.dt - surfaceAreaConsumedByVesicles);	
		if(getRadius() > replicationRadius) replicate();
	}
	
	protected void replicate() {
		getScene().addBacterium(new BSimBacterium(getPosition(), getRadius()/2, direction, getScene()));
		setRadius(getRadius()/2);
	}	

	protected void run() {				
		double shortTermMean = BSimUtils.mean(shortTermMemory());
		double longTermMean = BSimUtils.mean(longTermMemory());
		
		if(Math.random() < runProb(shortTermMean, longTermMean)) {
			Vector3d f = new Vector3d();
						
			if(shortTermMean - longTermMean > sensitivity) f.scale(BSimParameters.bactForceUp, direction);			
			else f.scale(BSimParameters.bactForceDown, direction);
			
			this.addForce(f);
			direction.set(getForce());
			direction.normalize();				
		}
		else switchMotionState();
	}
	
	protected double runProb(double shortTermMean, double longTermMean) {				
		if(shortTermMean - longTermMean > sensitivity) return runProbUp;		
		else if(longTermMean - shortTermMean > sensitivity) return runProbDown;		
		else return runProbIso;
	}	
		
	protected Vector<Double> longTermMemory() {
		return (Vector<Double>)memory.subList(0, (int)(longTermMemoryDuration/BSimParameters.dt));
	}
	
	protected Vector<Double> shortTermMemory() {
		return (Vector<Double>)memory.subList((int)(longTermMemoryDuration/BSimParameters.dt), memory.size());
	}	
	
	protected void tumble() {		
		// Obtain a random direction perpendicular to current direction		
		Vector3d randomVector = new Vector3d(Math.random(),Math.random(),Math.random());
		Vector3d crossVector = new Vector3d();
		crossVector.cross(direction, randomVector);		
		
		// Generate the rotation matrix for rotating about this direction by the tumble angle
		Matrix3d r = new Matrix3d();
		r.set(new AxisAngle4d(crossVector, tumbleAngle/tumbleSteps));
		
		// Apply the rotation			
		r.transform(direction);		
		
		// Decrement the tumble steps and tumble angle
		tumbleSteps-= 1;
		tumbleAngle-= tumbleAngle/tumbleSteps;
		
		if (tumbleSteps == 0) switchMotionState();		
	}
	
		
	protected void switchMotionState() {			
		if (motionState == RUNNING) {					
			tumbleSteps = (int)Math.ceil(BSimUtils.expRandVar(0.14)/BSimParameters.dt);			
			tumbleAngle = approxTumbleAngle();
			motionState = TUMBLING;
		} else if (motionState == TUMBLING) {
			motionState = RUNNING;		
		}
	}
			

	public Vector3d getDirection (){ return direction; }	
	public void removeVesicle(BSimVesicle v){ vesicles.remove(v); }
	
	/**
	 * Approximates the new tumble angle based on gamma distributed RV.
	 */
	protected double approxTumbleAngle() {
		double angle;
		
		// Select a random element from the array gammaVals
		int r = (int)Math.floor(Math.random()*1000.0);
		if(r>999) r=999;
		
		angle = Math.toRadians(gammaVals[r]);
		
		// check size and update sign if required
		if(Math.random()>=0.5) return angle;
		else return -angle;
	}
		
	/**
     * Read in file containing gamma distributed values to precision 0.001
     */
    public static double[] readGammaVals() {
    	double[] g = new double[1000];
    	int i = 0;    	
    	try {
    		Scanner scan = new Scanner(new File("src/bsim/resource/gammaVals.txt"));
    		while (scan.hasNextLine()) {
    			g[i] = Double.parseDouble(scan.nextLine());		 
    			i++;
    		}
    	} catch (IOException ex){
    		ex.printStackTrace();
    	}
    	return g;
    }	
}