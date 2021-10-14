package BsimAtiT6SS;
import java.util.Collections;
import bsim.BSim;
import bsim.BSimTicker;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;
import bsim.physical.contact.*;

import java.awt.geom.Area;
import java.util.ArrayList;
import java.util.Random;

import javax.vecmath.Vector2d;
import javax.vecmath.Vector3d;

public class AtiT6SSTicker extends BSimTicker {

    BSim sim;
    // logs data about time taken by ticker every LOG_INTERVAL timesteps
    final int LOG_INTERVAL;
    /** Whether to enable growth in the ticker etc. or not... */
    private static boolean WITH_GROWTH = true;
    
    static Random bacRng;
    //growth rate standard deviation
    public final double growth_stdv;
    //growth rate mean
    public final double growth_mean;
    //elongation threshold standard deviation
    public final double length_stdv;
    //elongation threshold mean
    public final double length_mean;

    final ArrayList<AtiT6SSBacterium> attacker_bac;
    final ArrayList<AtiT6SSBacterium> suscp_bac;
    final ArrayList<BSimCapsuleBacterium> bacteriaAll;
    final ArrayList<AtiT6SSBacterium> bac_bornbacAttacker;
    final ArrayList<AtiT6SSBacterium> bac_bornbacSuscp;
    final ArrayList<AtiT6SSBacterium> bac_deadbacAttacker;
    final ArrayList<AtiT6SSBacterium> bac_deadbacSuscp;
    
	// Initial Conditions
    // For a single screen
    final int MIXED_CONC = 1;
    final int CHECKER_BOARD = 2;
    int SINGLE_SCREEN = 2;
    
    /** Defines the progress of the chemical field flowing through the boundary on the x-axis. */
    int endpoint_x = 0;
    int field_box_num = 50;

    // internal machinery - don't worry about this
    final Mover mover;

    public AtiT6SSTicker(BSim sim, ArrayList<AtiT6SSBacterium> attacker_bac, ArrayList<AtiT6SSBacterium> suscp_bac,
    		ArrayList<BSimCapsuleBacterium> bacteriaAll, int LOG_INTERVAL, Random bacRng, 
    		double growth_stdv, double growth_mean, double length_stdv, double length_mean) {
        this.sim = sim;
        this.LOG_INTERVAL = LOG_INTERVAL;
        this.bacRng = bacRng; //random number generator
        this.growth_stdv = growth_stdv;
        this.growth_mean = growth_mean;
        this.length_stdv = length_stdv;
        this.length_mean = length_mean;
        this.attacker_bac = attacker_bac;
        this.suscp_bac = suscp_bac;
        this.bacteriaAll = bacteriaAll;
        bac_bornbacAttacker = new ArrayList();
        bac_bornbacSuscp = new ArrayList();
        bac_deadbacAttacker = new ArrayList();
        bac_deadbacSuscp = new ArrayList();
        mover = new RelaxationMoverGrid(bacteriaAll, sim);
        
    }

	/** Sets the flag for growth. **/
    public void setGrowth(boolean b) { WITH_GROWTH = b; }
    
    /** Function for bacteria growth activities. */
    public void growAttacker( ArrayList<AtiT6SSBacterium> bac, ArrayList<AtiT6SSBacterium> bacBorn ) {
    	Random bacRng = new Random(); 			// Random number generator
    	bacRng.setSeed(50); 					// Initializes random number generator
    	
        for (AtiT6SSBacterium b : bac) { 				// Loop over bac array
            b.grow();

            // Divide if grown past threshold
            if (b.L >= b.L_th) {
            	AtiT6SSBacterium daughter = b.divide();
            	
            	bacBorn.add(daughter);  		// Add daughter to newborn class, 'mother' keeps her status
            }
        }
        
        bac.addAll(bacBorn); 					// Adds all the newborn daughters to sub-population
        bacteriaAll.addAll(bacBorn); 			// Adds all the newborn daughters to total population
        
        // Allow daughter cells to grow 
        for ( AtiT6SSBacterium b : bacBorn ) {
        	
            // Assigns a division length to each bacterium according to a normal distribution
            double lengthThreshold = length_stdv*bacRng.nextGaussian()+length_mean;
            b.setElongationThreshold(lengthThreshold);
        }
        
        bacBorn.clear(); 						// Cleared for next time-step       
    }
    
    
    public void growSuscp( ArrayList<AtiT6SSBacterium> bac, ArrayList<AtiT6SSBacterium> bacBorn ) {
    	Random bacRng = new Random(); 			// Random number generator
    	bacRng.setSeed(50); 					// Initializes random number generator
    	
        for (AtiT6SSBacterium b : bac) { 				// Loop over bac array
            b.grow();

            // Divide if grown past threshold
            if (b.L >= b.L_th) {
            	AtiT6SSBacterium daughter = b.divide();
            	
            	bacBorn.add(daughter);  		// Add daughter to newborn class, 'mother' keeps her status
            }
        }
        
        bac.addAll(bacBorn); 					// Adds all the newborn daughters to sub-population
        bacteriaAll.addAll(bacBorn); 			// Adds all the newborn daughters to total population
        
        // Allow daughter cells to grow 
        for ( AtiT6SSBacterium b : bacBorn ) {
        	
            // Assigns a division length to each bacterium according to a normal distribution
            double lengthThreshold = length_stdv*bacRng.nextGaussian()+length_mean;
            b.setElongationThreshold(lengthThreshold);
        }
        
        bacBorn.clear(); 						// Cleared for next time-step
    }    
    
    /** Function to remove bacteria due to cell death or by boundary. */
    public void removeBacteria( BSim sim, ArrayList<AtiT6SSBacterium> bac, ArrayList<AtiT6SSBacterium> bac_dead ) {

        for (AtiT6SSBacterium b : bac) {
        	
            // Kick out if past any boundary
        	// Bacteria out of bounds = dead
            if(b.position.x < 0 || b.position.x > sim.getBound().x || b.position.y < 0 || b.position.y > sim.getBound().y || b.position.z < 0 || b.position.z > sim.getBound().z){
            	bac_dead.add(b);
            } 
            // Remove cell after it shrinks and becomes too small
            if ( b.L <= 1 ) {
            	bac_dead.add(b);          	
            }
        }                  
        
        
        bac.removeAll(bac_dead);
        
        // Remove from total population
        bacteriaAll.removeAll(bac_dead);
        bac_dead.clear();
    }

    // This one is a bit long too. Let's break it up
    // 1. Begins an "action" -> this represents one timestep
    // 2. Tells each bacterium to perform their action() function
    // 3. Updates each chemical field in the simulation
    // 4. Bacteria are then told to grow()
    // 5. bacteria which are longer than their threshold are told to divide()
    // 6. forces are applied and bacteria move around
    // 7. bacteria which are out of bounds are removed from the simulation
    @Override
    public void tick() {
        // increase lifetimes of cells
        for (AtiT6SSBacterium b : attacker_bac) {
            b.lifetime++;
        }
        for (AtiT6SSBacterium b : suscp_bac) {
            b.lifetime++;
        }         
        
        
        /********************************************** Collision */
        /** If bacterium i from bacterium B type comes in contact with bacterium j from Bacterium A type, its growth rate will be zero.
         */        
        
        
        for (AtiT6SSBacterium b1 : suscp_bac) {
        	for (AtiT6SSBacterium b2 : attacker_bac) { 
                //Determine if there is any overlap between the bounding boxes
                if (b1.isColliding(b2,0) == true) {                      
                	b1.setK_growth(0);
                	b1.isContact(true);
                	System.out.println("Collision detected between susceptible bacteria " + b1.id  + " and attacker bacteria " + b2.id ); 
                }
        	}
    }                 
        
        /**        
    	
        /********************************************** Action */
    	
        long startTimeAction = System.nanoTime(); 	// Wall-clock time, for diagnosing timing

        // Calculates and stores the midpoint of the cell.
        // Bacteria does action at each time step
        for(BSimCapsuleBacterium b : bacteriaAll) {
            b.action(); 
        }

        // Outputs how long each step took, once every log interval.
        long endTimeAction = System.nanoTime();
        if((sim.getTimestep() % LOG_INTERVAL) == 0) {
            System.out.println("Action update for " + bacteriaAll.size() + " bacteria took " + (endTimeAction - startTimeAction)/1e6 + " ms.");
        }  

        /********************************************** Chemical fields */
        
        startTimeAction = System.nanoTime();
        
 

        endTimeAction = System.nanoTime();
        if((sim.getTimestep() % LOG_INTERVAL) == 0) {
            System.out.println("Chemical field update took " + (endTimeAction - startTimeAction)/1e6 + " ms.");
        }

        /********************************************** Growth related activities if enabled. */
        if (WITH_GROWTH) {
            // ********************************************** Growth and division
            startTimeAction = System.nanoTime(); 	// Start action timer

            growAttacker( attacker_bac, bac_bornbacAttacker );		// For sub-population A
            growSuscp( suscp_bac, bac_bornbacSuscp );		// For sub-population B

            // Prints out information about bacteria when u want it to
            endTimeAction = System.nanoTime();
            if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                System.out.println("Growth and division took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
            }
            
            /********************************************** Neighbor interactions */
            
            startTimeAction = System.nanoTime();

            mover.move();

            endTimeAction = System.nanoTime();
            if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                System.out.println("Wall and neighbour interactions took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
            }

            /********************************************** Boundaries/removal */
            startTimeAction = System.nanoTime();
            
            removeBacteria( sim, attacker_bac, bac_deadbacAttacker );		// For sub-population A
            removeBacteria( sim, suscp_bac, bac_deadbacSuscp );		// For sub-population B
          
            
            endTimeAction = System.nanoTime();
            if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                System.out.println("Death and removal took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
            }
        }

        startTimeAction = System.nanoTime();

        endTimeAction = System.nanoTime();
        if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
            System.out.println("Switch took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
        }

    }    
    
}
