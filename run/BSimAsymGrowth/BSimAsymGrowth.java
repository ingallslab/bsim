package BSimAsymGrowth;

import bsim.BSim;
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.export.BSimPngExporter;
import bsim.winter2021.Bacterium;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import processing.core.PConstants;
import processing.core.PGraphics3D;

import javax.vecmath.*;
import java.awt.*;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.List;
import java.io.File;

/**
 * This class simulates the growth of cells asymmetrically
 */
public class BSimAsymGrowth {
	
    /** Whether to enable growth in the ticker etc. or not... */
    private static final boolean WITH_GROWTH = true;
    
    /** Wall boundaries of the bacteria (+x, +y, +z, -x, -y, -z). */
    private final boolean [] wall_boundaries = {false, false, true, false, false, true};
    
    // Boundaries
    // Boolean flag: specifies whether any walls are needed
    @Parameter(names = "-fixedbounds", description = "Enable fixed boundaries. (If not, one boundary will be leaky as real uf chamber).")
    private boolean fixedBounds = false;
    
    // @parameter means an optional user-specified value in the command line
    // export mode means output appears
    @Parameter(names = "-export", description = "Enable export mode.")
    private boolean export = true;
    
    // Simulation setup parameters. Set dimensions in um
    @Parameter(names = "-dim", arity = 3, description = "The dimensions (x, y, z) of simulation environment (um).")
    public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[] {75.0, 50.0, 1.0} ));
    // 198.0, 159.0, 1.0	// onecell.csv
    
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
    @Parameter(names = "-ratio", arity = 1, description = "Ratio of initial populations (proportion of activators).")
    public double populationRatio = 0.0;

    //growth rate standard deviation
    @Parameter(names="-gr_stdv",arity=1,description = "growth rate standard deviation")
    public static double growth_stdv = 0.02;		
    //growth rate mean
    @Parameter(names="-gr_mean",arity=1,description = "growth rate mean")
    public static double growth_mean = 0.5;

    //elongation threshold standard deviation
    @Parameter(names="-len_stdv",arity=1,description = "elongation threshold standard deviation")
    public static double length_stdv = 0.1;
    //elongation threshold mean
    @Parameter(names="-len_mean",arity=1,description = "elongation threshold mean")
    public static double length_mean = 7.0;

    /** Main Function. 
     * This is the very first function that runs in the simulation.
     * This runs the simulation. */ 
    public static void main(String[] args) {

        // Creates new simulation data object
        // Initializing storage unit for simulation
    	BSimAsymGrowth bsim_ex = new BSimAsymGrowth();

        // Starts up JCommander which allows you to read options from the command line more easily
        // Command line stuff
        new JCommander(bsim_ex, args);

        // Begins the simulation
        bsim_ex.run();
    }
    
    /** Creates a new Bacterium object. */
    public static Bacterium createBacterium( BSim sim ) {
    	Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator
        
        // Changed position to random 
        Vector3d pos1 = new Vector3d(Math.random()*sim.getBound().x, 
				Math.random()*sim.getBound().y, 
				Math.random()*sim.getBound().z);
        
        Vector3d pos2 = new Vector3d(Math.random()*sim.getBound().x, 
				Math.random()*sim.getBound().y, 
				Math.random()*sim.getBound().z);
        
        // Creates a new bacterium object whose endpoints correspond to the above data
        Bacterium bacterium = new Bacterium(sim, pos1, pos2);

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

        // Assigns a growth rate and a division length to bacterium according to a normal distribution
        double growthRate = growth_stdv * bacRng.nextGaussian() + growth_mean;
        bacterium.setK_growth(growthRate);
        
        double lengthThreshold = length_stdv * bacRng.nextGaussian() + length_mean;
        bacterium.setElongationThreshold(lengthThreshold);

        return bacterium;
    }
    
    // This function does a lot. let's break it down
    // 1. Initializes simulation (how long it runs, etc)
    // 2. Reads initial position data for bacteria and creates new bacterium objects accordingly
    // 3. Runs simulation loop until the simulation time is up
    // |----> Updates bacteria positions according to the forces acting on them
    // |----> Logs all data from the simulation into an excel sheet
    // |----> Saves images of simulation
    public void run() {
    	
		/*********************************************************
		 * Define simulation domain size and start time
		 */
        final double simX = simDimensions.get(0);
        final double simY = simDimensions.get(1);
        final double simZ = simDimensions.get(2);

        // Saves the exact time when the simulation started, real wall clock time
        long simulationStartTime = System.nanoTime();

		/*********************************************************
		 * Create a new simulation object and set up simulation settings
		 */
        final BSim sim = new BSim();
        sim.setDt(0.05);				    // Set simulation timestep in time units (0.01)
        									// Let the time units be in hours
        sim.setSimulationTime(100);       	// Specified in time units, could also specify a termination condition elsewhere
        sim.setTimeFormat("0.00");		    // Time Format for display on images
        sim.setBound(simX, simY, simZ);		// Simulation domain Boundaries

		// NOTE - solid = false sets a periodic boundary condition. 
		// This overrides leakiness!
		sim.setSolid(true, true, true);		// Solid (true) or wrapping (false) boundaries
        // Leaky -> bacteria can escape from sides of six faces of the box. 
		// Only usable if fixedbounds allows it
		
		// Set the wall boundaries
		sim.setWallBoundaries(wall_boundaries);
		
		// Set a closed or open boundary for phage diffusion
		if ( fixedBounds ) {
			sim.setLeaky(false, false, false, false, false, false);
			sim.setLeakyRate(0, 0, 0, 0, 0, 0);
		}
		else {
			sim.setLeaky(true, true, true, true, false, false);
			sim.setLeakyRate(5, 5, 5, 5, 0, 0);
		}
		
        /*********************************************************
         * Create the bacteria
         */
		
        // Separate lists of bacteria in case we want to manipulate the species individually
        // If multiple subpopulations, they'd be initialized separately, they'd be kept in different
        // need an array for each subpopulation, members can be repeated.
        final ArrayList<Bacterium> bac = new ArrayList();
        
        /** Track all of the bacteria in the simulation, for use of common methods etc.
        A general class, no sub-population specifics. */
        final ArrayList<BSimCapsuleBacterium> bacteriaAll = new ArrayList();

        Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator

        // Gets the location of the file that is currently running
        // Specify output file path
        String systemPath = new File("").getAbsolutePath()+"\\SingleCellSims";
        
        /** Empty list which will later contain the endpoint of rectangle positions 4 = x1,y1,x2,y2. */
        /*double[][] initEndpoints = new double[4][];

        // Creates a new csvreader object which can extract data from .csv files
        BufferedReader csvReader = null;
        try {
            // try reading the initial position file
            csvReader = new BufferedReader(new FileReader("C:\\Users\\sheng\\Documents\\CO-OP_W2021\\Research_Assistant\\Catie_BSim\\onecell.csv"));
        } catch (FileNotFoundException e) {
            // If that doesn't work, print out an error
            e.printStackTrace();
        }
 
        // If loading the content works, try reading the data
        // Updated: Previous version gave Exception in thread "main" java.lang.NumberFormatException: empty String
        try {
        	
            // Goes through each row of the excel sheet and pulls out the initial positions
            String row = csvReader.readLine();
            int i = 0;
            while (row != null) {
            	
                // Row.split takes a single line of the excel sheet and chops it up into the columns
                // mapToDouble takes the values in those columns and converts them to Java double data format
                // toArray converts the data into an array
            	
            	// In case of empty string, first check that the row is not blank
            	if (!row.isBlank()) {
            		initEndpoints[i] = Arrays.stream(row.split(",")).mapToDouble(Double::parseDouble).toArray();
            		
            		i ++;
            	}
            	row = csvReader.readLine();
            }
            csvReader.close(); 				// Close the file once all data is extracted    
        }
        
        catch(IOException e) {
            e.printStackTrace(); 			// If there is an error, this will just print out the message
        }

        // Now that the data is extracted, we can create the bacterium objects
       for(int j = 0; j < initEndpoints[0].length; j++){
        	
            // Initializes the endpoints of each bacterium from the array of endpoints
            // z-dimension is a small number, randomly generated, not sure why.
            Vector3d x1 = new Vector3d(initEndpoints[0][j]/13.89,initEndpoints[1][j]/13.89,bacRng.nextDouble()*0.1*(simZ - 0.1)/2.0);
            Vector3d x2 = new Vector3d(initEndpoints[2][j]/13.89,initEndpoints[3][j]/13.89,bacRng.nextDouble()*0.1*(simZ - 0.1)/2.0);
            // Note: the endpoint positions are scaled by 13.89, since the images are a bit more than 2000 pixels wide
            // While the simulation is rougly 200 micrometers. the conversion factor ends up being 13.89
            // Pixel to um scaling

            // Creates a new bacterium object whose endpoints correspond to the above data
            Bacterium bac0 = createBacterium(sim, field, x1, x2);

            // Adds the newly created bacterium to our lists for tracking purposes
            bac.add(bac0); 			// For separate subpopulations
            bacteriaAll.add(bac0);  // For all cells
        	
        } */
        
        /** Creating random bacterium objects to demonstrate phage field */

		// Create new phage sensing bacteria objects randomly in space
		while( bacteriaAll.size() < initialPopulation) {
            
            // Creates a new bacterium object whose endpoints correspond to the above data
			Bacterium bacterium = createBacterium(sim);
			
			// Adds the newly created bacterium to our lists for tracking purposes
			bac.add(bacterium); 			// For separate sub-populations
			bacteriaAll.add(bacterium);		// For all cells	
		}

        // Set up stuff for growth. Placeholders for the recently born and dead
        final ArrayList<Bacterium> bac_born = new ArrayList();
        final ArrayList<Bacterium> bac_dead = new ArrayList();

        // Internal machinery - dont worry about this line
        // Some kind of initialize of mover
        final Mover mover;
        mover = new RelaxationMoverGrid(bacteriaAll, sim);

        /*********************************************************
         * Set up the ticker
         */
        final int LOG_INTERVAL = 100; 				// Logs data every 100 timesteps (1)
        
        // This one is a bit long too. Let's break it up
        // 1. Begins an "action" -> this represents one timestep
        // 2. Tells each bacterium to perform their action() function
        // 3. Updates each chemical field in the simulation
        // 4. Bacteria are then told to grow()
        // 5. bacteria which are longer than their threshold are told to divide()
        // 6. forces are applied and bacteria move around
        // 7. bacteria which are out of bounds are removed from the simulation
        BSimTicker ticker = new BSimTicker() {
            @Override
            public void tick() {
            	
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
                    //System.out.println("Action update for " + bacteriaAll.size() + " bacteria took " + (endTimeAction - startTimeAction)/1e6 + " ms.");
                }  

                /********************************************** Chemical fields */
                
                startTimeAction = System.nanoTime();

                endTimeAction = System.nanoTime();
                if((sim.getTimestep() % LOG_INTERVAL) == 0) {
                    //System.out.println("Chemical field update took " + (endTimeAction - startTimeAction)/1e6 + " ms.");
                }

                /********************************************** Growth related activities if enabled. */
                if (WITH_GROWTH) {

                    // ********************************************** Growth and division
                    startTimeAction = System.nanoTime(); 	// Start action timer

                    for (Bacterium b : bac) { 				// Loop over bac array
                        b.grow();

                        // Divide if grown past threshold
                        if (b.L >= b.L_th) {
                        	Bacterium daughter = b.divide();
                        	bac_born.add(daughter);  		// Add daughter to newborn class, 'mother' keeps her status
                        }
                    }
                    
                    bac.addAll(bac_born); 					// Adds all the newborn daughters
                    bacteriaAll.addAll(bac_born); 			// Adds all the newborn daughters
                    
                    // Fixes the bug where daughter cells stop growing
                    for ( Bacterium b : bac_born ) {
                    	
                        // Assigns a growth rate and a division length to each bacterium according to a normal distribution
                        double growthRate = growth_stdv*bacRng.nextGaussian() + growth_mean;
                        b.setK_growth(growthRate);

                        double lengthThreshold = length_stdv*bacRng.nextGaussian()+length_mean;
                        b.setElongationThreshold(lengthThreshold);
                    }
                    
                    bac_born.clear(); 						// Cleared for next time-step

                    // Prints out information about bacteria when u want it to
                    endTimeAction = System.nanoTime();
                    if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                        //System.out.println("Growth and division took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
                    }
                    
                    /********************************************** Neighbor interactions */
                    
                    startTimeAction = System.nanoTime();

                    mover.move();

                    endTimeAction = System.nanoTime();
                    if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                        //System.out.println("Wall and neighbour interactions took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
                    }

                    /********************************************** Boundaries/removal */
                    startTimeAction = System.nanoTime();
                    
                    // Removal
                    for (Bacterium b : bac) {
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
                        //System.out.println("Death and removal took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
                    }
                }

                startTimeAction = System.nanoTime();

                endTimeAction = System.nanoTime();
                if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                    //System.out.println("Switch took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
                }

            }
        };

        sim.setTicker(ticker);

        /*********************************************************
         * Set up the drawer
         */
        BSimDrawer drawer = new BSimP3DDrawer(sim, 800, 600) {	
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
				
				// Draw the infected bacteria in red and the non-infected bacteria in green
				for(Bacterium b : bac) {
					draw(b, Color.GREEN);
				}		

            }
        };
        sim.setDrawer(drawer);
        
        export = false;
        if(export) {
            String simParameters = "" + BSimUtils.timeStamp() + "__dim_" + simX + "_" + simY + "_" + simZ
                    + "__ip_" + initialPopulation
                    + "__pr_" + populationRatio;

            if(fixedBounds){
                simParameters += "__fixedBounds";
            } else {
                simParameters += "__leakyBounds";
            }

            String filePath = BSimUtils.generateDirectoryPath(systemPath +"/" + simParameters + "/");
//            String filePath = BSimUtils.generateDirectoryPath("/home/am6465/tmp-results/" + simParameters + "/");

            /*********************************************************
             * Various properties of the simulation, for future reference.
             */
            BSimLogger metaLogger = new BSimLogger(sim, filePath + "simInfo.txt") {
                @Override
                public void before() {
                    super.before();
                    write("Simulation metadata.");
                    write("Catie Terrey Fall 2020."); //change name when new person :)
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
            metaLogger.setDt(10);//3600);			// Set export time step
            sim.addExporter(metaLogger);

            BSimLogger posLogger = new BSimLogger(sim, filePath + "position.csv") {
                DecimalFormat formatter = new DecimalFormat("###.##", DecimalFormatSymbols.getInstance( Locale.ENGLISH ));

                @Override
                public void before() {
                    super.before();
                    write("per Act; per Rep; id, p1x, p1y, p1z, p2x, p2y, p2z, growth_rate,directions");
                }

                @Override
                public void during() {
                    String buffer = new String();

                    buffer += sim.getFormattedTime() + "\n";
                    write(buffer);

                    write("acts");

                    buffer = "";
                    for(BSimCapsuleBacterium b : bacteriaAll) {
                        buffer += b.id + "," + formatter.format(b.x1.x)
                                + "," + formatter.format(b.x1.y)
                                + "," + formatter.format(b.x1.z)
                                + "," + formatter.format(b.x2.x)
                                + "," + formatter.format(b.x2.y)
                                + "," + formatter.format(b.x2.z)
                                + "," + formatter.format(b.getK_growth())
                                + "\n";
                    }

                    write(buffer);

                }
            };
            posLogger.setDt(10);			// set export time step for csv file
            sim.addExporter(posLogger);


            BSimLogger sumLogger = new BSimLogger(sim, filePath + "summary.csv") {


                @Override
                public void before() {
                    super.before();
                    write("time,id, status, p1x, p1y, p1z, p2x, p2y, p2z, px, py, pz, growth_rate, directions");
                }

                @Override
                public void during() {
                    String buffer = new String();
                    buffer = "";
                    for(BSimCapsuleBacterium b : bacteriaAll) {
                        buffer += sim.getFormattedTime()+","+b.id
                                + "," + b.getInfected()
                                + "," + b.x1.x
                                + "," + b.x1.y
                                + "," + b.x1.z
                                + "," + b.x2.x
                                + "," + b.x2.y
                                + "," + b.x2.z
                                + "," + b.position.x
                                + "," + b.position.y
                                + "," + b.position.z
                                + "," + b.getK_growth()
                                + "," + b.direction()
                                + "\n";
                    }

                    write(buffer);
                }
            };
            sumLogger.setDt(10);			// Set export time step
            sim.addExporter(sumLogger);

            /**
             * Export a rendered image file
             */
            BSimPngExporter imageExporter = new BSimPngExporter(sim, drawer, filePath );
            imageExporter.setDt(10); //this is how often it will output a frame
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

        //System.out.println("Total simulation time: " + (simulationEndTime - simulationStartTime)/1e9 + " sec.");
        
    }
    
    
}


