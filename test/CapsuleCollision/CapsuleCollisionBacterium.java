package CapsuleCollision;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.winter2021.Bacterium;

import javax.vecmath.Vector3d;
import java.util.Random;

/**
 */
public class CapsuleCollisionBacterium extends Bacterium {
	/** Initial growth rate of the bacteria. */
	final private double initial_growth_mean = 1.23; //um/h
	final private double initial_growth_stdv = 0.277; //um/h
	final private double initial_growth_rate;

    /** Function you call when you want to make a new bacterium object. */
    public CapsuleCollisionBacterium(BSim sim, Vector3d px1, Vector3d px2){
    	
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2);

    	Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator
        initial_growth_rate = initial_growth_stdv * bacRng.nextGaussian() + initial_growth_mean;
        setK_growth(initial_growth_rate);
    }

    /** Function you call when you want to make a bacterium object that is a child of another. */
    public CapsuleCollisionBacterium(BSim sim, Vector3d px1, Vector3d px2, long origin_id, long parent_id){
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2, origin_id, parent_id);

        Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator
        initial_growth_rate = initial_growth_stdv * bacRng.nextGaussian() + initial_growth_mean;
        setK_growth(initial_growth_rate);
    }

    // In case we want our bacteria to do anything special, we can add things to action() here that an ordinary
    // capsulebacterium wouldn't do.
    @Override
    public void action() { 						// Runs at every time step
        super.action();

    }

    // This function is called when the bacterium has passed its' division threshold and is ready to divide.
    public CapsuleCollisionBacterium divide() {
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
        CapsuleCollisionBacterium child = new CapsuleCollisionBacterium(sim, x1_child, new Vector3d(this.x2), this.origin_id, this.id);
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
