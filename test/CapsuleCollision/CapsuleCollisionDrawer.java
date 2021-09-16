package CapsuleCollision;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.draw.BSimP3DDrawer;
import processing.core.PConstants;
import processing.core.PGraphics3D;

import java.awt.*;
import java.util.ArrayList;

public class CapsuleCollisionDrawer extends BSimP3DDrawer{

    // Bacterium should be a sub-class of BSimCapsuleBacterium
    final ArrayList<CapsuleCollisionBacterium> bacA;
    final ArrayList<CapsuleCollisionBacterium> bacB;
    final double simX;
    final double simY;

    public CapsuleCollisionDrawer(BSim sim, int width, int height,
    		ArrayList<CapsuleCollisionBacterium> bac_to_drawA, ArrayList<CapsuleCollisionBacterium> bac_to_drawB) {
        super(sim, width, height);
        simX = width/10.0; //dividing makes the image bigger
        simY = height/10.0;
        bacA = bac_to_drawA;
        bacB = bac_to_drawB;

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

    @Override
    public void scene(PGraphics3D p3d) {
        p3d.ambientLight(128, 128, 128);
        p3d.directionalLight(128, 128, 128, 1, 1, -1);
		
		// Draw sub-population A
		for(CapsuleCollisionBacterium b : bacA) {
			draw(b, Color.GREEN );
		}	
		
		// Draw sub-population B
		for(CapsuleCollisionBacterium b : bacB) {
			draw(b, Color.ORANGE );
		}	

    }
}
