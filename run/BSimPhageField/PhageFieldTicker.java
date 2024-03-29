package BSimPhageField;

import bsim.BSim; 
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;

import java.util.ArrayList;
import java.util.Random;

public class PhageFieldTicker extends BSimTicker {

    BSim sim;
    // logs data about time taken by ticker every LOG_INTERVAL timesteps
    final int LOG_INTERVAL;
    
    /** Whether to enable growth in the ticker etc. or not... */
    public boolean WITH_GROWTH = true;
    /** Whether to flow phage in through the side of the boundary. */
    public static boolean FLOW_IN = true; 
    
    public static BSimChemicalField field;
    
    static Random bacRng;
    //growth rate standard deviation
    public final double el_stdv;
    //growth rate mean
    public final double el_mean;
    //elongation threshold standard deviation
    public final double div_stdv;
    //elongation threshold mean
    public final double div_mean;

    final ArrayList<PhageFieldBacterium> bac;
    final ArrayList<BSimCapsuleBacterium> bacteriaAll;
    final ArrayList<PhageFieldBacterium> bac_born;
    final ArrayList<PhageFieldBacterium> bac_dead;

    // internal machinery - don't worry about this
    final Mover mover;
    
    int field_box_num = 50;
    /** Defines the progress of the chemical field flowing through the boundary on the x-axis. */
    int endpoint_x = 0;

    public PhageFieldTicker(BSim sim, ArrayList<PhageFieldBacterium> bac, ArrayList<BSimCapsuleBacterium> bacteriaAll, int LOG_INTERVAL,
                       Random bacRng, double el_stdv, double el_mean, double div_stdv, double div_mean,
                       BSimChemicalField field, int field_box_num) {
        this.sim = sim;
        this.LOG_INTERVAL = LOG_INTERVAL;
        this.bacRng = bacRng; //random number generator
        this.el_stdv = el_stdv;
        this.el_mean = el_mean;
        this.div_stdv = div_stdv;
        this.div_mean = div_mean;
        this.bac = bac;
        this.bacteriaAll = bacteriaAll;
        bac_born = new ArrayList();
        bac_dead = new ArrayList();
        mover = new RelaxationMoverGrid(bacteriaAll, sim);
        this.field = field;
        this.field_box_num = field_box_num;
    }
    
    /** Sets the flag for flowing phage. **/
    public void setFlow(boolean b) { FLOW_IN = b; }
    /** Sets the flag for growth. **/
    public void setGrowth(boolean b) { WITH_GROWTH = b; }

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
        for (PhageFieldBacterium b : bac) {
            b.lifetime++;
        }
    	
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
        
        /** Allow the flow of phage through a boundary. */
        final double phage_conc = 2e3;   
        
        final int VER1 = 1;
        final int VER2 = 2;
        int version = 1;
        
        // State enabled where phage flows in through a boundary
        if ( FLOW_IN ) {
        	if ( version == VER1 ) {
            	for ( int i = 0; i < field_box_num; i ++ ) {
            		field.addQuantity( 0, i, 0, phage_conc * sim.getDt() );
            	}
        	}
        	else if ( version == VER2 ) {
        		final double conc = 1e2; 
            	if ( endpoint_x < field_box_num ) {
                	for ( int x = 0; x < endpoint_x; x ++ ) {
                		for ( int y = 0; y < field_box_num; y ++ ) {
                			field.addQuantity( x, y, 0, conc * sim.getDt() );
                		}
                	}
                	endpoint_x ++;
            	}
            	else {
            		endpoint_x = 0;
            	}
        	}
        }

        // Update the phage field
        field.update(); 

        endTimeAction = System.nanoTime();
        if((sim.getTimestep() % LOG_INTERVAL) == 0) {
            System.out.println("Chemical field update took " + (endTimeAction - startTimeAction)/1e6 + " ms.");
        }

        /********************************************** Growth related activities if enabled. */
        if (WITH_GROWTH) {

            // ********************************************** Growth and division
            startTimeAction = System.nanoTime(); 	// Start action timer

            for (PhageFieldBacterium b : bac) { 				// Loop over bac array
                b.grow();
                
                // Divide if grown past threshold
                if (b.L >= b.L_th) {
                	PhageFieldBacterium daughter = b.divide();
            		
                	// If the cell is infected, the daughter cell is also infected
                	if ( b.isInfected() ) {
                		daughter.setInfected(true);
                	}
                	
                	bac_born.add(daughter);  		// Add daughter to newborn class, 'mother' keeps her status
                }
            }
            
            bac.addAll(bac_born); 					// Adds all the newborn daughters
            bacteriaAll.addAll(bac_born); 			// Adds all the newborn daughters
            
            // Fixes the bug where daughter cells stop growing
            for ( PhageFieldBacterium b : bac_born ) {
            	
                // Assigns a growth rate and a division length to each bacterium according to a normal distribution
                double growthRate = el_stdv*bacRng.nextGaussian() + el_mean;
                b.setK_growth(growthRate);

                double lengthThreshold = div_stdv*bacRng.nextGaussian() + div_mean;
                b.setElongationThreshold(lengthThreshold);
            }
            
            bac_born.clear(); 						// Cleared for next time-step

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
            
            // Removal
            for (PhageFieldBacterium b : bac) {
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
            bacteriaAll.removeAll(bac_dead);
            bac_dead.clear();

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
