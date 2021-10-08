package testPhysicalContactBacterium;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import processing.core.PConstants;
import processing.core.PGraphics3D;

import java.awt.*;
import java.util.ArrayList;

public class PhysicalContactBacteriumDrawer extends BSimP3DDrawer{

    // Bacterium should be a sub-class of BSimCapsuleBacterium
    final ArrayList<TestPhysicalContactBacterium> bacA_bac;
    final ArrayList<TestPhysicalContactBacterium> bacB_bac;
    final double simX;
    final double simY;
    final double window_height;
    final double window_width;
    
    // Two ways to show the toxin fields on a single screen
    int SINGLE_SCREEN;
    final int CHECKER_BOARD = 1;
    final int MIXED_CONC = 2;
    
    public PhysicalContactBacteriumDrawer(BSim sim, double simX, double simY, int window_width, int window_height,
    		ArrayList<TestPhysicalContactBacterium> bac_to_drawA, ArrayList<TestPhysicalContactBacterium> bac_to_drawB,
    		int SINGLE_SCREEN) {
        super(sim, window_width, window_height);
        this.simX = simX;
        this.simY = simY;
        this.window_height = window_height;
        this.window_width = window_width;
        bacA_bac = bac_to_drawA;
        bacB_bac = bac_to_drawB;

        this.SINGLE_SCREEN = SINGLE_SCREEN;
    }

    /**
     * Draw the default cuboid boundary of the simulation as a partially transparent box
     * with a wireframe outline surrounding it.
     */
    @Override
    public void boundaries() {
        p3d.noFill();
        p3d.stroke(128, 128, 255);
        p3d.pushMatrix();
        p3d.translate((float)boundCentre.x,(float)boundCentre.y,(float)boundCentre.z);
        p3d.box((float)bound.x, (float)bound.y, (float)bound.z);
        p3d.popMatrix();
        p3d.noStroke();
    }
    
    /** Draws a colored box to screen. */
    public void legend( Color color, String title, int [] boxParams, int textX, int textY ) {
    	// Box
    	p3d.fill( color.getRed(), color.getGreen(), color.getBlue() );
    	p3d.stroke(128, 128, 255);
    	p3d.rect(boxParams[0], boxParams[1], boxParams[2], boxParams[3]);
    	
    	// Text
    	p3d.fill(50);
    	p3d.text(title, textX, textY);
    }
    
	/** Returns the associated color from the difference of two concentrations relative to concA. */
	public static Color getColor( double conc_A,double max_conc_difference ) {
		double conc = conc_A ;			// Concentration difference is relative to field A
		float value = 0.0f;						// Between 0 and 1
		
		// The floor of the hue is multiplied by 360 to get hue angle from HSB color model
		float max_hue = 0.1f;
		float min_hue = 0.9f;
		float mid_hue = (max_hue + min_hue) / 2;
		
		if ( conc >= max_conc_difference ) {
			value = 1/2f;// Blue 		
		}
		else if ( conc <= -max_conc_difference ) {
			value = 2/3f;// Green		
		}
		else {
			value = (float) ( (mid_hue)*(conc/max_conc_difference) );
		}
		
		float hue = value * max_hue + ( 1 - value ) * min_hue;
		return new Color( Color.HSBtoRGB(hue, 1f, 1f) );
	}
	
	/** Draw **/
    @Override
    public void draw(Graphics2D g) {
        p3d.beginDraw();

        if(!cameraIsInitialised){
            // camera(eyeX, eyeY, eyeZ, centerX, centerY, centerZ, upX, upY, upZ)
            p3d.camera((float)bound.x*0.5f, (float)bound.y*0.5f,
                    // Set the Z offset to the largest of X/Y dimensions for a reasonable zoom-out distance:
                    simX > simY ? (float)simX : (float)simY,
                    (float)bound.x*0.5f, (float)bound.y*0.5f, 0,
                    0,1,0);
            cameraIsInitialised = true;
        }

        p3d.textFont(font);
        p3d.textMode(PConstants.SCREEN);

        p3d.sphereDetail(10);
        p3d.noStroke();
        p3d.background(255, 255,255);
        

        // Draw only one simulation box is selected
            scene(p3d);
            boundaries();
            time();

        p3d.endDraw();
        g.drawImage(p3d.image, 0,0, null);
    }

    /**
     * Draw the formatted simulation time to screen.
     */
    @Override
    public void time() {
        p3d.fill(0);
        //p3d.text(sim.getFormattedTimeHours(), 50, 50);
        p3d.text(sim.getFormattedTime(), 50, 50);
    }
    
    /** Applies light settings to scene. */
    public void lights( PGraphics3D p3d ) {
        p3d.ambientLight(128, 128, 128);
        p3d.directionalLight(128, 128, 128, 1, 1, -1);
        //p3d.lightFalloff(1, 1, 0);
    }
    
    /** Draw Bacteria. */

    public void scene(PGraphics3D p3d) {
        p3d.ambientLight(128, 128, 128);
        p3d.directionalLight(128, 128, 128, 1, 1, -1);
		
        // Draw sub-population A
        for ( int i = 0; i < bacA_bac.size(); i ++ ) {
        	draw(bacA_bac.get(i), Color.GREEN);
        }
        
        // Draw sub-population B 
        for ( int i = 0; i < bacB_bac.size(); i ++ ) {
        	draw(bacB_bac.get(i),Color.CYAN);
        }
    }
}
