package BsimAtiT6SS;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import processing.core.PConstants;
import processing.core.PGraphics3D;

import java.awt.*;
import java.util.ArrayList;

public class AtiT6SSDrawer extends BSimP3DDrawer{

    // Bacterium should be a sub-class of BSimCapsuleBacterium
    final ArrayList<AtiT6SSBacterium> attacker_bac;
    final ArrayList<AtiT6SSBacterium> suscp_bac;
    final double simX;
    final double simY;
    final double window_height;
    final double window_width;
    
    // Two ways to show the toxin fields on a single screen
    int SINGLE_SCREEN;
    final int CHECKER_BOARD = 1;
    final int MIXED_CONC = 2;
    
    public AtiT6SSDrawer(BSim sim, double simX, double simY, int window_width, int window_height,
    		ArrayList<AtiT6SSBacterium> bac_to_attacker, ArrayList<AtiT6SSBacterium> bac_to_suscp,
    		int SINGLE_SCREEN) {
        super(sim, window_width, window_height);
        this.simX = simX;
        this.simY = simY;
        this.window_height = window_height;
        this.window_width = window_width;
        attacker_bac = bac_to_attacker;
        suscp_bac = bac_to_suscp;

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
            if ( SINGLE_SCREEN == CHECKER_BOARD ) {
            	legend( Color.MAGENTA, "Attacker", new int [] {-8, 0, 3, 3}, 57, 130 );
            	legend( Color.BLUE, "Susceptible", new int [] {-8, 6, 3, 3}, 57, 175 );
            }            
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
        for ( int i = 0; i < attacker_bac.size(); i ++ ) {
        	draw(attacker_bac.get(i), Color.GREEN);
        }
        
        // Draw sub-population B 
        for ( int i = 0; i < suscp_bac.size(); i ++ ) {
        	draw(suscp_bac.get(i),suscp_bac.get(i).isContactHappen()? Color.BLACK : Color.CYAN);
        }
    }
}
