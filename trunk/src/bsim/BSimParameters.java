/**
 * BSimParameters.java
 *
 * Class containing all simulation parameters of interest and get/set methods for each one
 * 
 * Authors: Ian Miles
 *          Thomas Gorochowski (Updates)
 *          Mattia Fazzini(Update)
 * Created: 05/08/2008
 * Updated: 09/08/2009
 */
package bsim;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.Vector;


public class BSimParameters {	

	// bsim
	public static double dt = 0.001; // seconds	
	public static int      screenWidth     	 	= 1025; // pixels
	public static int      screenHeight    		= 700; // pixels	
		
	// bsim.scene
	public static Vector<double[]> bacteria = new Vector();	
	public static Vector<double[]> coordBacteria = new Vector();
	public static Vector<double[]> recruitBacteria = new Vector();
	public static Vector<double[]> repBacteria = new Vector();
	public static Vector<double[]> sensingBacteria = new Vector();
	public static Vector<double[]> beads = new Vector();
	public static double[] fGoal;
	public static double[] fCoord;
	public static double[] fRecruit;
	public static double[] fQuorum;	
	public static double   screenZoom = 1; 
	public static double[] screenMove = {0.0, 0.0};	
	public static double bactRadius       = 1.4;	// microns
	public static double bactForceUp      = 0.4387; // pico newtons
	public static double bactForceDown    = 0.41; // pico newtons	
	public static double beadRadius       = 10.0;	// microns
	
	// bsim.particle
	public static double reactForce = 0.0;
	public static double visc = Math.pow(10.0,-3.0); // Pascal seconds
	public static double wellWidthBactBead = 0.0;
	public static double wellDepthBactBead = 0.0;
	
	// bsim.particle.bacterium	
	public static double runLengthUp      = 1.07;	// seconds
	public static double runLengthDown    = 0.8;	// seconds
	public static double runLengthIso     = 0.86; // seconds
			
	// bsim.field
	public static int      numOfThreads     = 2;	
		
	// bsim.batch
	public static int      dataFramesSkip   = 1;
	public static int      videoFramesSkip  = 1;
	public static int      simRuns          = 1;
	public static int      simLength        = 1;
	public static boolean  recordVideo      = true;
	public static int      frameRecordForSec	= 25;
	public static String   exportDir;	
				
	public BSimParameters(File f) {
		
		Scanner scanner;
		int lineCounter = 1;
		try { scanner = new Scanner(f);
			try {
				while(scanner.hasNextLine()) {
					processLine(scanner.nextLine().split("\t"), lineCounter);
					lineCounter++;
				}
			} finally { scanner.close(); }
		} catch(FileNotFoundException e) {System.err.println("Parameter file not found"); }
		
	}
		
	private void processLine(String[] line, int lineNo) {
		double[] args = parseLine(line);
				
		// bsim
		if     (line[0].equals("dt:")) dt = args[0];
		
		// bsim.scene
		else if(line[0].equals("bacterium:")) bacteria.add(args);
		else if(line[0].equals("coordBacterium:")) coordBacteria.add(args);
		else if(line[0].equals("recruitBacterium:")) recruitBacteria.add(args);
		else if(line[0].equals("repBacterium:")) repBacteria.add(args);
		else if(line[0].equals("sensingBacterium:")) sensingBacteria.add(args);
		else if(line[0].equals("bead:")) beads.add(args);
		else if(line[0].equals("fGoal:")) fGoal = args;
		else if(line[0].equals("fCoord:")) fCoord = args;
		else if(line[0].equals("fRecruit:")) fRecruit = args;
		else if(line[0].equals("fQuorum:")) fQuorum = args;				
		else if(line[0].equals("screenZoom:")) screenZoom = args[0];		
		else if(line[0].equals("screenMove:")) screenMove = args;
		else if(line[0].equals("bactRadius:")) bactRadius = args[0];
		else if(line[0].equals("bactForceUp:")) bactForceUp = args[0];
		else if(line[0].equals("bactForceDown:")) bactForceDown = args[0];	
		else if(line[0].equals("beadRadius:")) beadRadius = args[0];
				
		// bsim.particle
		else if(line[0].equals("reactForce:")) reactForce = args[0];
		else if(line[0].equals("visc:")) visc = args[0];
		else if(line[0].equals("wellWidthBactBead:")) wellWidthBactBead = args[0];
		else if(line[0].equals("wellDepthBactBead:")) wellDepthBactBead = args[0];
		
		// bsim.particle.bacteria		
		else if(line[0].equals("runLengthUp:")) runLengthUp = args[0];
		else if(line[0].equals("runLengthDown:")) runLengthDown = args[0];
		else if(line[0].equals("runLengthIso:")) runLengthIso = args[0];
				
		// bsim.field
		else if(line[0].equals("numOfThreads:")) numOfThreads = (int)args[0];		
		
		// bsim.app					
		else if(line[0].equals("screenWidth:")) screenWidth = (int)args[0];
		else if(line[0].equals("screenHeight:")) screenHeight = (int)args[0];
		
		// bsim.batch		
		else if(line[0].equals("dataFramesSkip:")) dataFramesSkip = (int)args[0];
		else if(line[0].equals("videoFramesSkip:")) videoFramesSkip = (int)args[0];
		else if(line[0].equals("simRuns:")) simRuns = (int)args[0];
		else if(line[0].equals("simLength:")) simLength = (int)args[0];
		else if(line[0].equals("recordVideo:")) { if((int)args[0] == 1) { recordVideo = true; } else { recordVideo = false; } }
		else if(line[0].equals("frameRecordForSec:")) frameRecordForSec = (int)args[0];
		else if(line[0].equals("exportDir:")) exportDir = line[1];
	
		else if(line[0].equals("***"))  {} // Do nothing
		else System.err.println("Line " + lineNo + " not Read in Parameter File");
	}
	
			
	/**
	 * Convert array of Strings into an array of doubles
	 */
	double[] parseLine(String[] line) {
		// Parse each line converting what can into doubles
		double[] parsedLine = new double[line.length-1];
		for(int i=1; i<line.length; i++) {
			try{
				parsedLine[i-1] = Double.parseDouble(line[i]);
			}
			catch(Exception e){
				// Could not convert type so leave as null;
				parsedLine[i-1] = 0.0;
			}
		}
		return parsedLine;
	}			

}
