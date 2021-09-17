package bsim.conjugation;

import bsim.BSim;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.winter2021.Bacterium;
import bsim.winter2021.RectangleIntersection;

import javax.vecmath.Vector2d;
import javax.vecmath.Vector3d;
import java.util.Random;

/**
 */
public class ConjugationBacterium extends Bacterium {

    /**
     * Conjugation related parameters.
     */
    final private double conjugation_contact_range = this.radius;
    final private static double donor_conjugation_frequency = 1.0;
    final private static double transconjugant_conjugation_frequency = 1.0;
    final private static double fitness_cost = 0.0; //fitness cost for carrying conjugating plasmid

    /**
     * These are the properties of the ConjugatingBacterium that we set
     */
    private String HGTstatus = ""; //determines whether the cell is donor/transconjugant/recipient; default is recipient
    private double contact_range = 0.0; //the default value is the contact range for recipients
    private double conjugation_frequency = 1.0; //conjugation attempts per some unit time

    /** Initial growth rate of the bacteria. */
	final private double initial_growth_mean = 1.23; //um/h
	final private double initial_growth_stdv = 0.277; //um/h
	final private double initial_growth_rate;

    /** Function you call when you want to make a new bacterium object.
     *  status is used such that the bacterium can be listed as a donor, transconjugant, or recipient
     *  */
    public ConjugationBacterium(BSim sim, Vector3d px1, Vector3d px2, String initialStatus){
    	
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2);

    	Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator
        initial_growth_rate = initial_growth_stdv * bacRng.nextGaussian() + initial_growth_mean;
        setK_growth(initial_growth_rate);

        if (initialStatus.equals("donor")) {
            this.HGTstatus = "donor";
            this.contact_range = conjugation_contact_range;
            this.conjugation_frequency = donor_conjugation_frequency;
            setK_growth(initial_growth_rate * (1.0 - fitness_cost)); // may be more effective to define donor, recipient, and transconjugant growth rates to avoid unnecessary calculations
        } else if (initialStatus.equals("transconjugant")) {
            this.HGTstatus = "transconjugant";
            this.contact_range = conjugation_contact_range;
            this.conjugation_frequency = transconjugant_conjugation_frequency;
            setK_growth(initial_growth_rate * (1.0 - fitness_cost));
        }
    }

    /** Function you call when you want to make a bacterium object that is a child of another.
     *  When a new bacterium is born in the simulation, the initialStatus should be set to this.HGTstatus
     */
    public ConjugationBacterium(BSim sim, Vector3d px1, Vector3d px2, long origin_id, long parent_id, String initialStatus){
        // "Bacterium" is a type of capsule bacterium, so we need to do all the things that a CapsuleBacterium does first.
        // This is the purpose of super(). The function super() initializes this bacterium object by first
        // referring to it as a capsulebacterium.
        super(sim, px1, px2, origin_id, parent_id);

        Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator
        initial_growth_rate = initial_growth_stdv * bacRng.nextGaussian() + initial_growth_mean;
        setK_growth(initial_growth_rate);

        if (initialStatus.equals("donor")) {
            this.HGTstatus = "donor";
            this.contact_range = conjugation_contact_range;
            this.conjugation_frequency = donor_conjugation_frequency;
            setK_growth(initial_growth_rate * (1.0 - fitness_cost)); // may be more effective to define donor, recipient, and transconjugant growth rates to avoid unnecessary calculations
        } else if (initialStatus.equals("transconjugant")) {
            this.HGTstatus = "transconjugant";
            this.contact_range = conjugation_contact_range;
            this.conjugation_frequency = transconjugant_conjugation_frequency;
            setK_growth(initial_growth_rate * (1.0 - fitness_cost));
        }
    }

    // In case we want our bacteria to do anything special, we can add things to action() here that an ordinary
    // capsulebacterium wouldn't do.
    @Override
    public void action() { 						// Runs at every time step
        super.action();

    }

    /**
     * This method is equivalent to getOrientedBoundingBoxIntersectionPerimeter() from winter2021, except that it can use unique contact ranges for both bacteria
     * @param neighbour_bac is the neighbouring bacterium
     * @return true if the boxes intersect, else false
     */
    protected boolean isColliding(ConjugationBacterium neighbour_bac) {
        boolean flag; //declare output

        /**
         * Contact ranges for this bacterium and neighbour
         */
        double r = this.contact_range;
        double r_neighbour = neighbour_bac.contact_range;

        //Algorithm for creating bounding boxes
        double bac_length = Math.sqrt(Math.pow(x2.x - x1.x, 2) + Math.pow(x2.y - x1.y, 2));
        double neighbour_bac_length = Math.sqrt(Math.pow(neighbour_bac.x2.x - neighbour_bac.x1.x, 2)
                + Math.pow(neighbour_bac.x2.y - neighbour_bac.x1.y, 2));
        double bac_angle = Math.atan((x2.y - x1.y) / (x2.x - x1.x));
        double neighbour_bac_angle = Math.atan((neighbour_bac.x2.y - neighbour_bac.x1.y) / (neighbour_bac.x2.x - neighbour_bac.x1.x));
        Vector2d[] rectangle_bac = RectangleIntersection.rectangle_vertices(position.x, position.y,
                bac_length + 2 * (radius + r), 2 * (radius + r), bac_angle);
        Vector2d[] rectangle_neighbour = RectangleIntersection.rectangle_vertices(neighbour_bac.position.x,
                neighbour_bac.position.y, neighbour_bac_length + 2 * (neighbour_bac.radius + r_neighbour),
                2 * (neighbour_bac.radius + r_neighbour), neighbour_bac_angle);

        //Determine if there is any overlap between the bounding boxes
        if (RectangleIntersection.rectangle_intersection_perimeter(rectangle_bac, rectangle_neighbour) > 1e-6) {
            flag = true;
        } else {
            flag = false;
        }
        return flag;
    }

    /**
     *
     * @param bacterium the bacterium whose HGT status and parameters are being set to transconjugant
     */
    public void setStatusToTransconjugant(ConjugationBacterium bacterium) {
            bacterium.HGTstatus = "transconjugant";
            bacterium.contact_range = conjugation_contact_range;
            bacterium.conjugation_frequency = transconjugant_conjugation_frequency;
            bacterium.setK_growth(initial_growth_rate * (1.0 - fitness_cost));
    }

    /**
     * This method determines if a successful conjugation event occurs, and what happens if conjugation occurs
     * The logic goes like this:
     * If conjugation in this bacterium is "activated", and...
     *      If this bacterium is in contact with a neighbour, and...
     *          If this bacterium is a donor/transconjugant and the neighbour is a recipient
     *             Turn the neighbour into a transconjugant!
     *
     * Sept 17, 2021: eventually, this will need to go through a list of neighbours and randomly pick only one if a conjugation event is successful
     */
    public void conjugate(ConjugationBacterium neighbour_bac) {
        //double conjugation_probability = this.conjugation_frequency*sim.getDt(); //uniform distribution
        double conjugation_probability = 1.0; //for now, assume 100% conjugation if donor contacts recipient
        if (Math.random() > (1.0 - conjugation_probability) ) {
            if (isColliding(neighbour_bac)) {
                if ((this.HGTstatus.equals("donor") || this.HGTstatus.equals("transconjugant")) && (neighbour_bac.HGTstatus.equals(""))) {
                    setStatusToTransconjugant(neighbour_bac);
                }
            }
        }
    }

    // This function is called when the bacterium has passed its' division threshold and is ready to divide.
    public ConjugationBacterium divide() {
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
        ConjugationBacterium child = new ConjugationBacterium(sim, x1_child, new Vector3d(this.x2), this.origin_id, this.id, this.HGTstatus);
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
