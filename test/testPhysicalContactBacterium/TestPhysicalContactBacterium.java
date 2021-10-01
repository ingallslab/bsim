package testPhysicalContactBacterium;

import bsim.BSim; 

import bsim.BSimChemicalField;
import bsim.physical.contact.PhysicalContactBacterium;
import bsim.winter2021.Bacterium;
import javax.vecmath.Vector3d;

import java.util.*;
/**
 * This class represents a bacteria for a cross-protection simulation. 
 * The growth rate of the bacteria depends on it's initial growth rate and the surrounding
 * concentration of toxins. 
 * If a cell is surrounded by an toxin it isn't resistant to, the growth rate will decrease. 
 * If the concentration of toxin around a bacteria is greater than a certain threshold, the cell 
 * will start to shrink and die (turns red in the simulation). 
 */
public class TestPhysicalContactBacterium extends PhysicalContactBacterium {
	/** Chemical field representing the toxin the bacteria is resistant to. */
	protected BSimChemicalField resistant_field;
	/** Chemical field representing the toxin the bacteria is not resistant to. */
	protected BSimChemicalField toxin_field;
	
	/** Initial growth rate of the bacteria. */
	public static double initial_el_mean = 0.2;
	public static double initial_el_stdv = 0.05;
	final private double initial_el_rate;
	
	/** Rate at which the bacteria shrinks when dying. **/
	final private double shrink_stdv = 0.134;			
	final private double shrink_mean = -0.117;
	final private double shrink_rate;
	
	/** Scales the growth rate of the bacteria. */
	final private double scale_1 = 1;//2;
	/** Scales the concentration for the growth rate of the bacteria. */
	final private double scale_2 = 1;
	
	/** Threshold of toxins for bacterium growth to slow [Molecules/(micron)^3]. */
	final private double growth_threshold;  
	/** Threshold of toxins for bacterium to start dying [Molecules/(micron)^3]. */
	final private double conc_threshold;  
	/** Flag for enzyme production. */
	private boolean production;			
	/** Flag to indicate toxin levels are above threshold. */
	private boolean aboveConcThreshold;
	/** Amount of enzymes produced to decay the toxin field. */
	public static double enzymeNum;

	
    // Function you call when you want to make a new bacterium object
    public TestPhysicalContactBacterium(BSim sim, BSimChemicalField toxin, BSimChemicalField resistant, Vector3d px1, Vector3d px2){
    	
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2);
        
        this.toxin_field = toxin;	
        this.resistant_field = resistant;
        
        growth_threshold = 50;
        conc_threshold = 5e2;	
        production = true;
        aboveConcThreshold = false;
        enzymeNum = 3e3;
        
    	Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator
        
        // Calculate the initial elongation rate 
        initial_el_rate = initial_el_stdv * bacRng.nextGaussian() + initial_el_mean;
        setK_growth(initial_el_rate);
        // Calculate the rate at which the bacteria will shrink
        shrink_rate = shrink_stdv * bacRng.nextGaussian() + shrink_mean;
    }

    // Function you call when you want to make a bacterium object that is a child of another
    public TestPhysicalContactBacterium(BSim sim, BSimChemicalField toxin, BSimChemicalField resistant, Vector3d px1, Vector3d px2, long origin_id, long parent_id){

        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2, origin_id, parent_id);

        this.toxin_field = toxin;
        this.resistant_field = resistant;

        growth_threshold = 1e2;
        conc_threshold = 3e2;
        production = true;
        aboveConcThreshold = false;
        enzymeNum = 3e3;

        Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator

        // Calculate the initial elongation rate
        initial_el_rate = initial_el_stdv * bacRng.nextGaussian() + initial_el_mean;
        setK_growth(initial_el_rate);
        // Calculate the rate at which the bacteria will shrink
        shrink_rate = shrink_stdv * bacRng.nextGaussian() + shrink_mean;
    }

    // In case we want our bacteria to do anything special, we can add things to action() here that an ordinary
    // capsulebacterium wouldn't do.
    @Override
    public void action() { 						// Runs at every time step
        super.action();
        
        // If toxin levels are above the threshold, the bacteria starts to die
        if(toxin_field!=null) {
            // If toxin levels are above the threshold, the bacteria starts to die
            if ( toxin_field.getConc(position) > getConcThreshold() ) {
            	aboveConcThreshold = true;
            	
            	// Bacteria starts shrinking 
            	Random bacRng = new Random();	
    	        setK_growth( shrink_rate );
            	//setK_growth(0.0);		// stops growing
    	        //setProduction(false); 	// stops producing enzymes when dying
            }
            // If toxin levels are below the threshold, 
            // cell growth rate is lowered depending on surrounding toxin concentration 
            // (values and scaling are arbitrary for now)
            else {
            	aboveConcThreshold = false;
                if ( toxin_field.getConc(position) > growth_threshold ) {
                	setK_growth( initial_el_rate * scale_1/(scale_1 + toxin_field.getConc(position) * scale_2 * sim.getDt()));
                } 
            }     	
        	
        }
        
		// Production of enzymes to decay the toxin the bacteris is resistant to
		if ( resistant_field != null ) {						// Steady production of enzymes (enzyme/hr)
			resistant_field.addQuantity(position, -enzymeNum * sim.getDt() );	
		}
		
    }
	
    /** Returns the toxin threshold for a bacterium. */
    public double getConcThreshold() { return this.conc_threshold; }
    /** Returns the flag for the internal enzyme production of a bacterium. */
    public boolean isProducing() { return this.production; }
    /** Sets the flag for the internal enzyme production of a bacterium. */
    public void setProduction( boolean production) { this.production = production; }
    /** Returns the flag that indicates toxin levels above threshold. */
    public boolean isAboveThreshold() { return aboveConcThreshold; }
    /** Sets the elongation mean. **/
    public void set_elMean(double el_mean) { this.initial_el_mean = el_mean; }
    /** Sets the elongation stdv. **/
    public void set_elStdv(double el_stdv) { this.initial_el_stdv = el_stdv; }
    
    /** Sets the elongation stdv. **/
    public void set_enzymNUM(double enzymeNum) { this.enzymeNum = enzymeNum; }    
    
    

    // This function is called when the bacterium has passed its' division threshold and is ready to divide.
    public TestPhysicalContactBacterium divide() {
        Vector3d randomVec = new Vector3d(rng.nextDouble(), rng.nextDouble(), rng.nextDouble());
        randomVec.scale(twist);
        System.out.println("Bacterium " + this.id + " is dividing...");

        Vector3d u = new Vector3d(); 
        u.sub(this.x2, this.x1);

        // Decides where to split bacterium
        double divPert = 0.1*L*(rng.nextDouble() - 0.5); // 0.1 is arbitrary scaling here.

        // Get length of bacterium
        double L_actual = u.length();

        // Computes lengths of child bacteria
        double L1 = L_actual*0.5*(1 + divPert) - radius;
        double L2 = L_actual*0.5*(1 - divPert) - radius;

        // Finds new endpoints of mother and daughter
        Vector3d x2_new = new Vector3d();
        x2_new.scaleAdd(L1/L_actual, u, this.x1);
        Vector3d longVec = new Vector3d();
        longVec.scaleAdd(-1,this.x2,this.x1); 			// Push along bacterium length
        longVec.scale(push*rng.nextDouble()); 			// Push is applied to bacterium
        
        // Impulse, not a force.
        longVec.add(randomVec);
        x2_new.add(longVec);

        Vector3d x1_child = new Vector3d();
        x1_child.scaleAdd(-(L2/L_actual), u, this.x2);
        Vector3d longVec2=new Vector3d();
        longVec2.scaleAdd(-1,longVec);
        x1_child.add(longVec2); 						// Push applied to child cell

        // Set the child cell.
        // Creates new bacterium called child and adds it to the lists, gives posns, infected status and chemical field status
        TestPhysicalContactBacterium child = new TestPhysicalContactBacterium(sim, toxin_field, resistant_field, x1_child, new Vector3d(this.x2), this.origin_id, this.id);
        this.parent_id = this.id;
        this.lifetime = 0;
        // this.initialise(L1, this.x1, x2_new); // for symmetric growth
        // Asymmetrical growth occurs at division node
        // so we need to swap x1 and x2 for the mother after division for asymmetrical elongation
        // This does not affect symmetric growth
        this.initialise(L1, x2_new, this.x1);
        child.L = L2;

        // add child to list of children - Sohaib Nadeem
        addChild(child);

        // Calculate angle between daughter cells at division
        angle_initial = coordinate(child);

        // Prints a line whenever a new bacterium is made
        System.out.println("Child ID is " + child.id);
        return child;
    }
}
