/**
 * BSimSolidBoxBoundary.java
 *
 * Class that represents a box of the simulation boundaries.
 *
 * Authors: Mattia Fazzini
 * Created: 13/08/2009
 */
package bsim.drawable.boundary;

import java.awt.Graphics;

import bsim.drawable.BSimDrawable;


public class BSimBoxBoundary implements BSimDrawable {

	
	//plane points
	protected double[] p1 = {0.0,0.0,0.0};
	
	protected double length = 0;
	protected double width = 0;
	protected double depth = 0;
	
	/**
	 * General constructor.
	 */
	public BSimBoxBoundary(double[] args) {
		
		// Update the internal variables
		p1[0] = args[0];
		p1[1] = args[1];
		p1[2] = args[2];
		
		length = args[3];
		width = args[4];
		depth = args[5];
	}	
	
	/**
	 * Standard get methods
	 */
	public double[] getP1() { return p1; }
	public double getLength() { return length; }
	public double getWidth() { return width; }
	public double getDepth() { return depth; }
	
	/**
	 * Method to return the centrePos of the box
	 */
	public double[] getCentrePos() {
		double[] centrePos= {0.0,0.0,0.0};
		centrePos[0]=p1[0]+(length/2);
		centrePos[1]=p1[1]+(width/2);
		centrePos[2]=p1[2]+(depth/2);
		return centrePos; 
		}
	
	
	/**
	 * Draw the boundary.
	 */
	public void redraw(Graphics g) {
	}
}