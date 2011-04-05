package bsim.export;

import java.io.File;
import java.util.Vector;

import javax.vecmath.Vector3d;

import bsim.particle.BSimBacterium;
import bsim.particle.BSimParticle;
import bsim.particle.BSimVesicleAcceptor;
import bsim.scene.BSimScene;


public class BSimParticleFullFileExport extends BSimFileExport {

	// Number of frames to skip between measurements
	private int frameSkip = 1;
	
	// The list to output at each timestep
	private Vector<BSimParticle> particles;
	
	
	/**
	 * General constructor that uses a given file as the bases for the object.
	 */
	public BSimParticleFullFileExport(Vector<BSimParticle> v, File f, int newFrameSkip) {
		
		super(f);

		// Update internal parameter
		frameSkip = newFrameSkip;
		particles = v;
	}
	
	
	/**
	 * Returns the final line (if necessary) in the file. Only use this for calculations
	 * where you only need a single output for the whole simulation.
	 */
	public String finalOutputLine(BSimScene scene){

		// No final line required
		return null;
	}

	
	/**
	 * Abstract method to get the header row for the file. Only called once on object
	 * creation.
	 */
	public String getHeaderLine(BSimScene scene) {
		String outStr;
		
		outStr = "Timestep";
		
		for(int i=0; i<particles.size(); i++){
			outStr += ", " + (i+1) + "_x";
			outStr += ", " + (i+1) + "_y";
			outStr += ", " + (i+1) + "_z";
			if(particles.get(i) instanceof BSimVesicleAcceptor) {
				outStr += ", " + (i+1) + "_fusionCount";
			}
		}
		
		outStr += ", Average_x";
		outStr += ", Average_y";
		outStr += ", Average_z";
		
		return outStr;
	}

	
	/**
	 * Returns the next output line based on a given scene and parameters.
	 */
	public String nextOutputLine(BSimScene scene) {
		String outStr;
		Vector3d particleAvg;
		double vecOutput[];
		double actTime;
		
		particleAvg = new Vector3d();
		vecOutput = new double[3];
		
		// Check to see if the frame needs to be output
		if(scene.getTimeStep() % frameSkip == 0 && particles.size() > 0){

			// Add the timestep
			actTime = scene.getTimeStep() * BSimScene.dt;
			outStr = "";
			outStr += "" + actTime;

			// Loop through all the particles and add their details to the average
			for(BSimParticle p : particles){
				
				p.getPosition().get(vecOutput);
				outStr += ", " + vecOutput[0];
				outStr += ", " + vecOutput[1];
				outStr += ", " + vecOutput[2];
				if(p instanceof BSimVesicleAcceptor) {
					outStr += ", " + ((BSimVesicleAcceptor)p).fusionCount;
				}
				
				// Calculate the sum of all bacteria positions
				particleAvg.add(p.getPosition());
			}

			// Calaculte the average positions and output
			particleAvg.scale(1.0/(double)particles.size());
			particleAvg.get(vecOutput);
			
			outStr += ", " + vecOutput[0];
			outStr += ", " + vecOutput[1];
			outStr += ", " + vecOutput[2];

			return outStr;
		}
		else{
			return null;
		}
	}
}