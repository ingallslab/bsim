package bsim.particle.vesicle;
import java.util.Random;

import javax.vecmath.Vector3d;

import bsim.BSimParameters;
import bsim.particle.BSimParticle;
import bsim.particle.bacterium.BSimBacterium;
import bsim.scene.BSimScene;

public class BSimVesicle extends BSimParticle {
	
	private double budRadius; // microns
	protected double radiusGrowthRate; // microns/sec
	
	public BSimVesicle(double newRadius) {
		super(newRadius);
	}
	
	public BSimVesicle(Vector3d newPosition, double newRadius, BSimScene newScene) {
		super(newPosition, newRadius, newScene);	
		budRadius = 0.1;
		radiusGrowthRate = 0.001;
	}

	public void action() {
		
		Random r = new Random();
		
		double resistance = 6.0*Math.PI*(this.getRadius()*Math.pow(10, -6))*BSimParameters.visc;
		double boltzmann = 1.38 * Math.pow(10,-23);
		double temperature = 300;
		double amplitude = Math.sqrt(2*resistance*boltzmann*temperature/BSimParameters.dt)*Math.pow(10,12);
		
		Vector3d f = new Vector3d(r.nextGaussian()*amplitude, r.nextGaussian()*amplitude, r.nextGaussian()*amplitude);		
		this.addForce(f);
	}
	
	public double grow(BSimBacterium bacterium) {		
		double oldSurfaceArea = getSurfaceArea();
		setRadius(getRadius() + radiusGrowthRate*BSimParameters.dt);
		double newSurfaceArea = getSurfaceArea();
				
		if(getRadius() > budRadius)
			bud(bacterium);
		
		return newSurfaceArea - oldSurfaceArea;
	}
	
	public void bud(BSimBacterium bacterium) {		
		bacterium.removeVesicle(this);
		escape(bacterium);		
		getScene().addVesicle(this);
	}
	
}
