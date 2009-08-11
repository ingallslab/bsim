/**
 * BSimCoordBacterium.java
 *
 * Class that represents a bacterium that by default will move randomly, until contact 
 * with a bead is made at which time it will follow the goal chemoattractant. A 
 * co-ordination signal will also be released that will cause any bacteria in a high 
 * enough concentration to also follow the chemoattractant.
 *
 * Authors: Thomas Gorochowski
 * Created: 28/08/2008
 * Updated: 28/08/2008
 */
package bsim.drawable.bacteria;

import java.awt.Color;
import java.awt.Graphics;

import bsim.BSimParameters;
import bsim.BSimScene;
import bsim.drawable.BSimDrawable;
import bsim.logic.BSimLogic;


public class BSimCoordBacterium extends BSimSensingBacterium implements BSimLogic, BSimDrawable {

	
	// Threshold for detecting co-ordination signal (AHL)
	protected double coordThreshold = 0.1;


	/**
	 * General constructor.
	 */
	public BSimCoordBacterium(double newSpeed, double newMass, double newSize,
			double[] newDirection, double[] newPosition, double newForceMagnitudeDown,
			double newForceMagnitudeUp,
			int newState, double newTumbleSpeed, int newRemDt, BSimScene newScene, 
		    BSimParameters newParams, double newSwitchSpeed, double newCoordThreshold) {

		// Call the parent constructor with the basic properties	
		super(newSpeed, newMass, newSize, newDirection, newPosition, newForceMagnitudeDown,
		newForceMagnitudeUp, newState,
		      newTumbleSpeed, newRemDt, newScene, newParams, newSwitchSpeed);
		
		coordThreshold = newCoordThreshold;
	}


	/**
	 * Implements the BSimLogic interface. In this case it merely carries out
	 * the standard chemotaxis toward fGoal gradient. The internal force of the bacterium
	 * at a timestep is returned.
	 */
	public double[] runLogic ( boolean contactBac, 
	                           boolean contactBead,
	                           boolean contactBoundary ) {
		
		if(beadContactTimer > 0){
			beadContactTimer--;
			
			// Generate some co-ordination chemical (AHL) at current location
			scene.getCoordinationField().addChemical (1.0, this.getCentrePos());
		}
		
		if(contactBead){
			beadContactTimer = (int)(switchSpeed / params.getDtSecs());
		}
		
		return  super.runLogic(contactBac, contactBead, contactBoundary);
	}

	
	/**
	 * This is an updated version of the BSimBacterium method to only allow for the
	 * sensed goal concentration to be used if in contact with a bead.
	 */
	protected double senseRunContinueProb() {
		double shortTermMean;
		double longTermMean;
		double shortTermCounter = 0.0;
		double longTermCounter = 0.0;
		double shortTermMemoryLength = 1.0; // seconds
		double longTermMemoryLength = 3.0; // seconds
		double sensitivity = 0.000001;
		
		// Check to see if the bacteria has been in contact with a bead
		if(beadContactTimer > 0 || 
		   scene.getCoordinationField().getConcentration(this.getCentrePos()) > coordThreshold){
			// Perform the normal attraction to the goal chemoattractant
			for(int i=0; i<concMemory.size();i++) {
				if(i <= (longTermMemoryLength/params.getDtSecs())) {
					longTermCounter = longTermCounter + (Double)concMemory.elementAt(i);
				} else shortTermCounter = shortTermCounter + (Double)concMemory.elementAt(i);
			}
			shortTermMean = shortTermCounter / (1 + (shortTermMemoryLength/params.getDtSecs()));
			longTermMean = longTermCounter / (longTermMemoryLength/params.getDtSecs());
		
			if(shortTermMean - longTermMean > sensitivity) {
				runUp = true;
				return upRunProb;
			}
			else if(longTermMean - shortTermMean > sensitivity) {
				runUp = false;
				return downRunProb;
			}
			else {
				runUp = false;
				return isoRunProb;
			}
		}
		else{
			// If not perform random walk
			return isoRunProb;
		}
	}
	
	
	/**
	 * Redraws the bacterium. A small red circle is also drawn to represent the direction
	 * of the bacteria.
	 */
	public void redraw(Graphics g) {

		// Draw the main shape of bacterium
		if(beadContactTimer > 0){
			g.setColor(Color.BLUE);
		}
		else if(scene.getCoordinationField().getConcentration(this.getCentrePos()) > coordThreshold){
			g.setColor(Color.YELLOW);
		}
		else{
			g.setColor(Color.GREEN);
		}
		
		g.fillOval((int)position[0],(int)position[1],(int)(size),(int)(size));

		// Draw an indicator of bacterium's direction
		int x1,x2;
		double littleR = size/5.0;
		x1 = (int)(position[0] + (size/2.0)*(1+direction[0]) - (littleR/Math.sqrt(2.0)));
		x2 = (int)(position[1] + (size/2.0)*(1+direction[1]) - (littleR/Math.sqrt(2.0)));
		g.setColor(Color.RED);
		g.fillOval(x1,x2,(int)(littleR*2.0),(int)(littleR*2.0));
	}
}