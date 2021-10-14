package BsimAtiT6SS;

import bsim.BSim; 

import bsim.physical.contact.PhysicalContactBacterium;
import bsim.winter2021.Bacterium;
import javax.vecmath.Vector3d;

import java.util.*;
/**
 * This class represents a collision between Attacker and susceptible bacteria.
 * If bacterium i from bacterium Susceptible comes in contact with bacterium j from Bacterium Attacker, its growth rate will be zero.
 */
public class AtiT6SSBacterium extends PhysicalContactBacterium {
	
	/** Initial growth rate of the bacteria. */
	public static double initial_el_mean = 0.2;
	public static double initial_el_stdv = 0.05;
	final private double initial_el_rate;
	
	/** Rate at which the bacteria shrinks when dying. **/
	final private double shrink_stdv = 0.134;			
	final private double shrink_mean = -0.117;
	final private double shrink_rate;
	
	/** Flag to indicate toxin levels are above threshold. */
	private boolean contact=false;    
	
	
    // Function you call when you want to make a new bacterium object
    public AtiT6SSBacterium(BSim sim,Vector3d px1, Vector3d px2){
    	
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2);

    	Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator
        
        // Calculate the initial elongation rate 
        initial_el_rate = initial_el_stdv * bacRng.nextGaussian() + initial_el_mean;
        setK_growth(initial_el_rate);
        // Calculate the rate at which the bacteria will shrink
        shrink_rate = shrink_stdv * bacRng.nextGaussian() + shrink_mean;
    }

    // Function you call when you want to make a bacterium object that is a child of another
    public AtiT6SSBacterium(BSim sim,Vector3d px1, Vector3d px2, long origin_id, long parent_id){

        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2, origin_id, parent_id);;

        Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator

        // Calculate the initial elongation rate
        initial_el_rate = initial_el_stdv * bacRng.nextGaussian() + initial_el_mean;
        setK_growth(initial_el_rate);
        // Calculate the rate at which the bacteria will shrink
        shrink_rate = shrink_stdv * bacRng.nextGaussian() + shrink_mean;
    }

    /** Sets the elongation mean. **/
    public void set_elMean(double el_mean) { this.initial_el_mean = el_mean; }
    /** Sets the elongation stdv. **/
    public void set_elStdv(double el_stdv) { this.initial_el_stdv = el_stdv; }   
    
    

    // This function is called when the bacterium has passed its' division threshold and is ready to divide.
    public AtiT6SSBacterium divide() {
        Vector3d randomVec = new Vector3d(rng.nextDouble(), rng.nextDouble(), rng.nextDouble());
        randomVec.scale(twist);

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
        AtiT6SSBacterium child = new AtiT6SSBacterium(sim,x1_child, new Vector3d(this.x2), this.origin_id, this.id);
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

        return child;
    }     
 
    
    //check contact is happened or not
    public void isContact(boolean contact) {
        this.contact = contact;
    } 
    
    /** Returns the flag that indicates contact is happened. */
    public boolean isContactHappen() { return contact; }    
    
    
}
