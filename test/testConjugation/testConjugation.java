package testConjugation;

import bsim.BSim;
import bsim.BSimUtils;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;
import bsim.conjugation.*;
import bsim.export.BSimLogger;
import bsim.export.BSimPngExporter;
import bsim.winter2021.Bacterium;
import bsim.winter2021.RawReader;
import bsim.winter2021.SimReader;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import javax.vecmath.Vector3d;
import java.io.File;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;

/**
 *
 * This class simulates conjugation between two bacteria
 * Currently, conjugation will occur if physical contact is made between a donor and recipient
 *
 */
public class testConjugation {

    private static String executionPath = testConjugation.class.getProtectionDomain().getCodeSource().getLocation().getPath();
    private static String projectPath = System.getProperty("user.dir");
    private static String pathToSimulation = projectPath+"/test/"+"testConjugation/";

    @Parameter(names = "-input_data", description = "path to input data")
    private static String input_data = pathToSimulation+"positions.csv";

    // Boundaries
    // Boolean flag: specifies whether any walls are needed
    @Parameter(names = "-fixedbounds", description = "Enable fixed boundaries. (If not, one boundary will be leaky as real uf chamber).")
    private boolean fixedBounds = true;

    // @parameter means an optional user-specified value in the command line
    // export mode means output appears
    @Parameter(names = "-export", description = "Enable export mode.")
    private boolean export = true;

    @Parameter(names = "-export_path", description = "export location")
    private String export_path = pathToSimulation+"/results";

    @Parameter(names = "-grow", description = "Enable bacteria growth.")
    private boolean WITH_GROWTH = true;

    // Simulation setup parameters. Set dimensions in um
    final double width = 50.0;
    final double height = 50.0;

    @Parameter(names = "-dim", arity = 3, description = "The dimensions (x, y, z) of simulation environment (um).")
    //public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[] {75.0, 50.0, 1.0} ));
    public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[]{50.0, 50.0, 1.}));

    // Grid ->
    // 52x42 -> 546
    // 100x86 -> 2150
    // Random:
    // 50x42 -> 250
    // 100x85 -> 1000
    // Density (cell number)
    // optional call to a default initial set of cells
    @Parameter(names = "-pop", arity = 1, description = "Initial seed population (n_total).")
    public int initialPopulation = 2;

    // A:R ratio
    // for default set of cells, set ratio of two subpopulations
    @Parameter(names = "-ratio", arity = 1, description = "Ratio of initial populations.")
    public double populationRatio = 0.5;

    //growth rate standard deviation
    @Parameter(names="-el_stdv",arity=1,description = "elongation rate standard deviation")
    public static double el_stdv = 0.277;
    @Parameter(names="-el_mean",arity=1,description = "elongation rate mean")
    public static double el_mean = 1.23;

    //elongation threshold standard deviation
    @Parameter(names="-div_stdv",arity=1,description = "elongation threshold standard deviation")
    public static double div_stdv = 0.1;
    //elongation threshold mean
    @Parameter(names="-div_mean",arity=1,description = "elongation threshold mean")
    public static double div_mean = 7.0;

    // Simulation Time
    @Parameter(names="-simt",arity=1,description = "simulation time")
    public static double sim_time = 5.0;
    @Parameter(names="-simdt",arity=1,description = "simulation time step")
    public static double sim_dt = 0.01;
    @Parameter(names="-export_time",arity=1,description = "export time")
    public static double export_time = 0.1;// Previously was 10, and simulation time was 100

    // internal force
    @Parameter(names="-k_int",arity=1,description = "internal force")
    public static double k_int = 50.0;
    // cell-cell collision force
    @Parameter(names="-k_cell",arity=1,description = "cell-cell collision force")
    public static double k_cell = 500.0;
    // sticking force
    @Parameter(names="-k_stick",arity=1,description = "side-to-side attraction")
    public static double k_sticking = 0.01;

    // sticking range
    @Parameter(names="-rng_stick",arity=1,description = "max range side-to-side attraction")
    public static double range_sticking = 5.0;

    // twist
    @Parameter(names="-twist",arity=1,description = "twist")
    public static double twist = 0.1;
    // push
    @Parameter(names="-push",arity=1,description = "push")
    public static double push = 0.05;

    // asymmetric growth threshold
    @Parameter(names="-l_asym",arity=1,description = "asymmetric growth threshold")
    public static double L_asym = 3.75;
    // value of asymmetry
    @Parameter(names="-asym",arity=1,description = "asymmetry")
    public static double asymmetry = 0.1;
    // symmetric growth
    @Parameter(names="-sym",arity=1,description = "symmetric growth")
    public static double sym_growth = 0.05;

    // Separate lists of bacteria in case we want to manipulate the species individually
    // If multiple subpopulations, they'd be initialized separately, they'd be kept in different
    // need an array for each subpopulation, members can be repeated.
    final ArrayList<ConjugationBacterium> donorBac = new ArrayList();
    final ArrayList<ConjugationBacterium> transconjugantBac = new ArrayList();
    final ArrayList<ConjugationBacterium> recipientBac = new ArrayList();

    /** Track all of the bacteria in the simulation, for use of common methods etc.
    A general class, no sub-population specifics. */
    final ArrayList<ConjugationBacterium> bac = new ArrayList();
    final ArrayList<BSimCapsuleBacterium> bacteriaAll = new ArrayList();

    // Set up stuff for growth. Placeholders for the recently born and dead
    final ArrayList<ConjugationBacterium> bac_born = new ArrayList();
    final ArrayList<ConjugationBacterium> bac_dead = new ArrayList();

    /** Main Function.
     * This is the very first function that runs in the simulation.
     * This runs the simulation. */
    public static void main(String[] args) {

        // Creates new simulation data object
        // Initializing storage unit for simulation
    	testConjugation bsim_ex = new testConjugation();

        // Starts up JCommander which allows you to read options from the command line more easily
        // Command line stuff
        new JCommander(bsim_ex, args);

        // Begins the simulation
        bsim_ex.run();
    }

    public static ConjugationBacterium createBacterium(BSim sim, Vector3d x1, Vector3d x2, String HGTstatus, double growthRate, double lengthThreshold) {
        // creates a new bacterium object whose endpoints correspond to the above data
        ConjugationBacterium bacterium = new ConjugationBacterium(sim, x1, x2, HGTstatus);

        // determine the vector between the endpoints
        // if the data suggests that a bacterium is larger than L_max, that bacterium will instead
        // be initialized with length 1. Otherwise, we set the length in the following if statement
        // Earlier they were all initialised to length 1.
        Vector3d dispx1x2 = new Vector3d();
        dispx1x2.sub(x2, x1); // sub is subtract
        double length = dispx1x2.length(); // determined.
        if (length < bacterium.L_max) {
            bacterium.initialise(length, x1, x2); // redundant to record length, but ok.
        }

        // assign the growth rate and division length to the bacterium
        bacterium.setK_growth(growthRate);
        bacterium.setElongationThreshold(lengthThreshold);

        return bacterium;
    }

    /** Creates a new Bacterium object with random position*/
    public static ConjugationBacterium createRandomBacterium(BSim sim, String HGTstatus) {
    	Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator

        /**
         * Use this to initialize bacteria from random positions
         */
        ArrayList<Vector3d> positions = BSimUtils.randomPosition(sim, div_mean);
        Vector3d pos1 = positions.get(0);
        Vector3d pos2 = positions.get(1);

        ConjugationBacterium bacterium = new ConjugationBacterium(sim, pos1, pos2, HGTstatus);

        // Determine the vector between the endpoints
        // If the data suggests that a bacterium is larger than L_max, that bacterium will instead
        // be initialized with length 1. Otherwise, we set the length in the following if statement
        // Earlier they were all initialized to length 1.
        Vector3d dispx1x2 = new Vector3d();
        dispx1x2.sub(pos2, pos1); 						// Sub is subtract
        double length = dispx1x2.length(); 				// Determined.

        if (length < bacterium.L_max) {
            bacterium.initialise(length, pos1, pos2); 	// Redundant to record length, but ok.
        }

        // Assigns a division length to bacterium according to a normal distribution
        double lengthThreshold = div_stdv * bacRng.nextGaussian() + div_mean;
        bacterium.setElongationThreshold(lengthThreshold);

        // Assigns the specified forces, range, and impulses
        bacterium.setIntForce(k_int);
        bacterium.setCellForce(k_cell);
        bacterium.setStickForce(k_sticking);
        bacterium.setStickingRange(range_sticking);
        bacterium.setTwist(twist);
        bacterium.setPush(push);
        bacterium.setAsym(asymmetry);

        return bacterium;
    }

    public void run() {

		/*********************************************************
		 * Define simulation domain size and start time
		 */
        final double simX = simDimensions.get(0);
        final double simY = simDimensions.get(1);
        final double simZ = simDimensions.get(2);

        long simulationStartTime = System.nanoTime();

		/*********************************************************
		 * Create a new simulation object and set up simulation settings
		 */
        final BSim sim = new BSim();
        sim.setDt(sim_dt);				    // Set simulation timestep in time units
        									// Let the time units be in hours
        sim.setSimulationTime(sim_time);    // Specified in time units, could also specify a termination condition elsewhere
        sim.setTimeFormat("0.00");		    // Time Format for display on images
        sim.setBound(simX, simY, simZ);		// Simulation domain Boundaries

		// NOTE - solid = false sets a periodic boundary condition.
		// This overrides leakiness!
		sim.setSolid(true, true, true);		// Solid (true) or wrapping (false) boundaries

        // Leaky -> bacteria can escape from sides of six faces of the box.
		// Only usable if fixedbounds allows it
		// Set a closed or open boundary for antibiotic diffusion
		if ( fixedBounds ) {
			sim.setLeaky(false, false, false, false, false, false);
			sim.setLeakyRate(0, 0, 0, 0, 0, 0);
		}
		else {
			sim.setLeaky(true, true, true, true, false, false);
			sim.setLeakyRate(1, 1, 1, 1, 0, 0);
		}

        /*********************************************************
         * Create the bacteria
         */

        Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator

        // Gets the location of the file that is currently running
        // Specify output file path
        String systemPath = executionPath+"/results"; //BS...

        //String systemPath = new File("").getAbsolutePath();

        // TODO: HANDLE DEFAULT INPUT DATA FILE
        String initial_data_path = input_data;
        SimReader reader = new SimReader();

        ArrayList<double[]> cell_endpoints = reader.readcsv(initial_data_path);
        for(int i = 0; i < cell_endpoints.size(); i++) {//double[] cell : cell_endpoints) {
            double[] cell = cell_endpoints.get(i);
            // initializes the endpoints of each bacterium from the array of endpoints; z-dimension is 0.5
            Vector3d x1 = new Vector3d(cell[0], cell[1], 0.5);
            Vector3d x2 = new Vector3d(cell[2], cell[3], 0.5);
            // assigns a growth rate and a division length to bacterium according to a normal distribution
            // assigns a growth rate and a division length to bacterium according to a normal distribution
            double growthRate = el_stdv * bacRng.nextGaussian() + el_mean;
            double divThreshold = div_stdv * bacRng.nextGaussian() + div_mean;
            ConjugationBacterium bac0 = createBacterium(sim, x1, x2, "", growthRate, divThreshold);
            // adds the newly created bacterium to our lists for tracking purposes
            bac.add(bac0); // for separate subpopulations
            bacteriaAll.add(bac0);  // for all cells
        }

        bac.get(1).HGTstatus = "donor";

        final Mover mover;
        mover = new RelaxationMoverGrid(bacteriaAll, sim);

        /*********************************************************
         * Set up the ticker
         */
        final int LOG_INTERVAL = 10; // logs data every X timesteps
        ConjugationTicker ticker = new ConjugationTicker(sim, bacteriaAll, bac, LOG_INTERVAL, bacRng,
        		el_stdv, el_mean, div_stdv, div_mean);
        ticker.setGrowth(WITH_GROWTH);
        sim.setTicker(ticker);

        /*********************************************************
         * Set up the drawer
         */
        int drawerWidth = 800;
        int drawerHeight = 800;
        ConjugationDrawer drawer = new ConjugationDrawer(sim, drawerWidth, drawerHeight,
        		bac);
        sim.setDrawer(drawer);

        /*********************************************************
         * Export data
         */
        if(export) {
            String simParameters = "" + BSimUtils.timeStamp() + "__dim_" + simX + "_" + simY + "_" + simZ
                    + "__ip_" + initialPopulation
                    + "__pr_" + populationRatio;

            if(fixedBounds){
                simParameters += "__fixedBounds";
            } else {
                simParameters += "__leakyBounds";
            }

            String filePath;
            if(export_path == "default") {
                filePath = BSimUtils.generateDirectoryPath(systemPath +"/" + simParameters + "/");
                //String filePath = BSimUtils.generateDirectoryPath("/home/am6465/tmp-results/" + simParameters + "/");
            } else {
                filePath = BSimUtils.generateDirectoryPath(export_path + "/");
            }

            /*********************************************************
             * Various properties of the simulation, for future reference.
             */
            BSimLogger metaLogger = new BSimLogger(sim, filePath + "simInfo.txt") {
                @Override
                public void before() {
                    super.before();
                    write("Simulation metadata.");
                    write("Aaron Yip"); //change name when new person :)
                    write("Simulation dimensions: (" + simX + ", " + simY + ", " + simZ + ")");
                    write("Initial population: "+ initialPopulation);
                    write("Ratio " + populationRatio);

                    if(fixedBounds){
                        write("Boundaries: fixed");
                    } else {
                        write("Boundaries: leaky");
                    }
                }

                @Override
                public void during() {

                }
            };
            metaLogger.setDt(export_time);			// Set export time step
            sim.addExporter(metaLogger);

            BSimLogger sumLogger = new BSimLogger(sim, filePath + "summary.csv") {

                @Override
                public void before() {
                    super.before();
                    write("time, id, HGTstatus, p1x, p1y, p1z, p2x, p2y, p2z, px, py, pz, directions");
                }

                @Override
                public void during() {
                    String buffer = new String();
                    buffer = "";
                    for(ConjugationBacterium b : bac) {
                        buffer += sim.getFormattedTime()+","+b.id
                                + "," + b.HGTstatus
                                + "," + b.x1.x
                                + "," + b.x1.y
                                + "," + b.x1.z
                                + "," + b.x2.x
                                + "," + b.x2.y
                                + "," + b.x2.z
                                + "," + b.position.x
                                + "," + b.position.y
                                + "," + b.position.z
                                + "," + b.direction()
                                + "\n";
                    }

                    write(buffer);
                }
            };
            sumLogger.setDt(export_time);			// Set export time step
            sim.addExporter(sumLogger);

            /**
             * Export a rendered image file
             */
            BSimPngExporter imageExporter = new BSimPngExporter(sim, drawer, filePath );
            imageExporter.setDt(export_time); //this is how often it will output a frame
            // separate time-resolution for images vs excel file
            sim.addExporter(imageExporter);

            sim.export();

            /**
             * Drawing a java plot once we're done?
             * See TwoCellsSplitGRNTest
             */

        }

        // Run the simulation
        else {
            sim.preview();
        }

        long simulationEndTime = System.nanoTime();

        System.out.println("Total simulation time: " + (simulationEndTime - simulationStartTime)/1e9 + " sec.");

    }


}
