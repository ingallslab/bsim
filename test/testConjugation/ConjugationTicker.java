package testConjugation;

import bsim.BSim;
import bsim.BSimTicker;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.conjugation.*;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;

import java.util.ArrayList;
import java.util.Random;

public class ConjugationTicker extends BSimTicker {

    BSim sim;
    // logs data about time taken by ticker every LOG_INTERVAL timesteps
    final int LOG_INTERVAL;
    public static boolean WITH_GROWTH;
    static Random bacRng;
    //growth rate standard deviation
    public final double growth_stdv;
    //growth rate mean
    public final double growth_mean;
    //elongation threshold standard deviation
    public final double length_stdv;
    //elongation threshold mean
    public final double length_mean;

    final ArrayList<ConjugationBacterium> donorBac;
    final ArrayList<ConjugationBacterium> transconjugantBac;
    final ArrayList<ConjugationBacterium> recipientBac;
    final ArrayList<BSimCapsuleBacterium> bacteriaAll;
    final ArrayList<ConjugationBacterium> donorBac_born;
    final ArrayList<ConjugationBacterium> transconjugantBac_born;
    final ArrayList<ConjugationBacterium> recipientBac_born;
    final ArrayList<ConjugationBacterium> donorBac_dead;
    final ArrayList<ConjugationBacterium> transconjugantBac_dead;
    final ArrayList<ConjugationBacterium> recipientBac_dead;

    // internal machinery - don't worry about this
    final Mover mover;

    public ConjugationTicker(BSim sim, ArrayList<ConjugationBacterium> donorBac, ArrayList<ConjugationBacterium> transconjugantBac, ArrayList<ConjugationBacterium> recipientBac,
    		ArrayList<BSimCapsuleBacterium> bacteriaAll, int LOG_INTERVAL, Random bacRng, 
    		double growth_stdv, double growth_mean, double length_stdv, double length_mean) {
        this.sim = sim;
        this.LOG_INTERVAL = LOG_INTERVAL;
        this.bacRng = bacRng; //random number generator
        this.growth_stdv = growth_stdv;
        this.growth_mean = growth_mean;
        this.length_stdv = length_stdv;
        this.length_mean = length_mean;

        this.donorBac = donorBac;
        this.transconjugantBac = transconjugantBac;
        this.recipientBac = recipientBac;
        this.bacteriaAll = bacteriaAll;
        donorBac_born = new ArrayList();
        transconjugantBac_born = new ArrayList();
        recipientBac_born = new ArrayList();
        donorBac_dead = new ArrayList();
        transconjugantBac_dead = new ArrayList();
        recipientBac_dead = new ArrayList();
        mover = new RelaxationMoverGrid(bacteriaAll, sim);
    }
    
    /** Sets the flag for growth. **/
    public void setGrowth(boolean b) { WITH_GROWTH = b; }
    
    /** Function for bacteria growth activities. */
    public void growBacteria( ArrayList<ConjugationBacterium> bac, ArrayList<ConjugationBacterium> bacBorn ) {
    	Random bacRng = new Random(); 			// Random number generator
    	bacRng.setSeed(50); 					// Initializes random number generator
    	
        for (ConjugationBacterium b : bac) { 				// Loop over bac array
            b.grow();

            // Divide if grown past threshold
            if (b.L >= b.L_th) {
            	ConjugationBacterium daughter = b.divide();
            	bacBorn.add(daughter);  		// Add daughter to newborn class, 'mother' keeps her status
            }
        }
        
        bac.addAll(bacBorn); 					// Adds all the newborn daughters to sub-population
        bacteriaAll.addAll(bacBorn); 			// Adds all the newborn daughters to total population
        
        // Allow daughter cells to grow 
        for ( ConjugationBacterium b : bacBorn ) {
        	
            // Assigns a division length to each bacterium according to a normal distribution
            double lengthThreshold = length_stdv*bacRng.nextGaussian()+length_mean;
            b.setElongationThreshold(lengthThreshold);
        }
        
        bacBorn.clear(); 						// Cleared for next time-step
    }
    
    /** Function to remove bacteria due to cell death or by boundary. */
    public void removeBacteria( BSim sim, ArrayList<ConjugationBacterium> bac, ArrayList<ConjugationBacterium> bac_dead ) {

        for (ConjugationBacterium b : bac) {
        	
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
        for (ConjugationBacterium b : donorBac) {
            b.lifetime++;
        }
        for (ConjugationBacterium b : transconjugantBac) {
            b.lifetime++;
        }
        for (ConjugationBacterium b : recipientBac) {
            b.lifetime++;
        }

        //More efficient would be a KdTree search with a range of the bounding box. Eventually, it will be good to do this anyways to get other relationships related to neighbours.
        //TODO: see the BSimTicker in BSimTutorial.java for an example of how they detect collisions! Implement here
/*        for(int i = 1; i < bacteria.size(); i++) {
            for(int j = i+1; j < bacteria.size(); j++) {
                bacteria.get(i).interaction(bacteria.get(j));
            }
        }

        // Each bacterium performs its action and updates its position
        for(BSimTutorialBacterium b : bacteria) {
            b.action();
            b.updatePosition();
        }*/

        /**
         * Another implementation... from IteratorMover.java
         * I believe this one prevents double-counting pairs, so will need to update the conjugation condition in ConjugationBacterium to include donor-recipient
         */

/*        // Interaction
        int j = 1;
        for (BSimCapsuleBacterium b1 : allBacteria) {
            for (BSimCapsuleBacterium b2 : allBacteria.subList(j, bacteriaSize)) {
//                            System.out.println("Force between " + b1 + " and " + b2);
                b1.computeNeighbourForce(b2);
            }
            j += 1;
        }*/

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

        /********************************************** Growth related activities if enabled. */
        if (WITH_GROWTH) {

            // ********************************************** Growth and division
            startTimeAction = System.nanoTime(); 	// Start action timer

            growBacteria( donorBac, donorBac_born);
            growBacteria( transconjugantBac, transconjugantBac_born);
            growBacteria( recipientBac, recipientBac_born );

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
            
            removeBacteria( sim, donorBac, donorBac_dead);
            removeBacteria( sim, donorBac, transconjugantBac_dead);
            removeBacteria( sim, recipientBac, recipientBac_dead );
            
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
