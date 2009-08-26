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

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.Vector;

import bsim.field.BSimChemicalField;
import bsim.field.BSimChemicalFieldCreate;
import bsim.particle.bacterium.BSimBacteriaCreate;
import bsim.particle.bead.BSimBeadsCreate;
import bsim.rendering.visualaid.BSimVisualAidCreate;


public class BSimParameters {


	public double   beadRadius       = 10.0;	// microns

	public double   bactRadius       = 1.4;	// microns
	public double   bactForceUp      = 0.4387; // pico newtons
	public double   bactForceDown    = 0.41; // pico newtons
	public double   bactSpeed        = 50.0; // microns per second
	
	public double   runLengthUp      = 1.07;	// seconds
	public double   runLengthDown    = 0.8;	// seconds
	public double   runLengthIso     = 0.86; // seconds
	
	public double   visc             = Math.pow(10.0,-3.0); // Pascal seconds
	public int      screenWidth      = 1025; // pixels
	public int      screenHeight     = 700; // pixels

	public double   dt               = 0.001; // seconds
	
	public Vector bacteriaSingles 	  = new Vector(); 
	public Vector bacteriaSets 		  = new Vector();
	public Vector beadSingles 		  = new Vector();
	public Vector beadSets			  = new Vector(); 	              	               
	public Vector vaBacteriaTraces 	  = new Vector();
	public Vector vaAvgBacteriaTraces = new Vector();
	public Vector vaBeadTraces 		  = new Vector();
	public Vector vaClocks 			  = new Vector();
	public Vector vaScales			  = new Vector();

	public double[] cfGoalDefine = {0, 0, 0, 10, 10, 10, 10, 10, 10, 0.001, 1, 1, 0.8};
	public double[] cfGoalSetup = {0, 0, 0, 0, 0, 0};
	public double[] cfCoordDefine = {0, 0, 0, 10, 10, 10, 10, 10, 10, 0.001, 1, 1, 0.8};
	public double[] cfCoordSetup = {0, 0, 0, 0, 0, 0};
	public double[] cfRecruitDefine = {0, 0, 0, 10, 10, 10, 10, 10, 10, 0.001, 1, 1, 0.8};
	public double[] cfRecruitSetup = {0, 0, 0, 0, 0, 0};
	public double[] cfQuorumDefine = {0, 0, 0, 10, 10, 10, 10, 10, 10, 0.001, 1, 1, 0.8};
	public double[] cfQuorumSetup = {0, 0, 0, 0, 0, 0};
	
	public double[] boundingBoxDefine = {0, 0, 0, 0, 0, 0};	
	
	public double[] magnStrength = {0.0, 0.0, 0.0};
	
	public double   screenZoom       = 1, 
	                screenMove[]     = {0.0, 0.0};
	
	//parameters needed to control the processing camera
	public double   minimumDistance  = 0.001;
	public double   maximumDistance  = 1500;
	public double   defaultDistance  = 1000;
	public int      frameForSec      = 25;
	public int      frameRecordForSec= 25;
	
	public int      dataFramesSkip   = 1;
	public int      videoFramesSkip  = 1;
	public int      simRuns          = 1;
	public int      simLength        = 1;
	public boolean  recordVideo      = true;
	public String   exportDir;
	public int      numOfThreads     = 2;
	
	public double   wellWidthBactBact = 0.0;
	public double   wellDepthBactBact = 0.0;
	public double   wellWidthBactBead = 0.0;
	public double   wellDepthBactBead = 0.0;
	public double   wellWidthBeadBead = 0.0;
	public double   wellDepthBeadBead = 0.0;
	public double   wellWidthBeadBdry = 0.0;
	public double   wellDepthBeadBdry = 0.0;
	public double   wellWidthBactBdry = 0.0;
	public double   wellDepthBactBdry = 0.0;
	public double   wellWidthVesBdry = 0.0;
	public double   wellDepthVesBdry = 0.0;
	public double   wellWidthVesBead = 0.0;
	public double   wellDepthVesBead = 0.0;
	
	public double	reactForce = 0.0;
		
	
	public BSimParameters() {	
	}
		
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
				
		if     (line[0].equals("DT:")) setDtSecs(args[0]);
		
		else if(line[0].equals("BEAD_RADIUS:")) setBeadRadius(args[0]);
		else if(line[0].equals("CREATE_BEAD_SINGLE:")) addSingleBead(args);
		else if(line[0].equals("CREATE_BEAD_SET:"))	addBeadSet(args);

		else if(line[0].equals("BACTERIA_RADIUS:")) setBactRadius(args[0]);
		else if(line[0].equals("CREATE_BACTERIUM_SINGLE:")) addSingleBacterium(args);
		else if(line[0].equals("CREATE_BACTERIA_SET:")) addBacteriaSet(args);
		else if(line[0].equals("BACTERIA_FORCE_UP:")) setBactForceUp(args[0]);
		else if(line[0].equals("BACTERIA_FORCE_DOWN:")) setBactForceDown(args[0]);		
		else if(line[0].equals("UP_RUN_LENGTH:")) setUpRunLength(args[0]);
		else if(line[0].equals("DOWN_RUN_LENGTH:")) setDownRunLength(args[0]);
		else if(line[0].equals("ISO_RUN_LENGTH:")) setIsoRunLength(args[0]);
			
		else if(line[0].equals("VISCOSITY:")) setViscosity(args[0]);
		
		else if(line[0].equals("SCREEN_HEIGHT:")) setScreenHeight((int)args[0]);
		else if(line[0].equals("SCREEN_WIDTH:")) setScreenWidth((int)args[0]);
		else if(line[0].equals("SCREEN_ZOOM:")) setScreenZoom((double)args[0]);
		else if(line[0].equals("SCREEN_MOVE:")) setScreenMove((double)args[0], (double)args[1]);
		
		else if(line[0].equals("SCREEN_MINIMUM_DISTANCE:")) setMinimumDistance((double)args[0]);
		else if(line[0].equals("SCREEN_MAXIMUM_DISTANCE:")) setMaximumDistance((double)args[0]);
		else if(line[0].equals("SCREEN_DEFAULT_DISTANCE:")) setDefaultDistance((double)args[0]);
		else if(line[0].equals("SCREEN_FRAME_FOR_SECOND:")) setFrameForSec((int)args[0]);
		else if(line[0].equals("SCREEN_FRAME_RECORD_FOR_SECOND:")) setFrameRecordForSec((int)args[0]);
		
		else if(line[0].equals("PHYSICS_WELL_WIDTH_BACT_BACT:")) setWellWidthBactBact((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_DEPTH_BACT_BACT:")) setWellDepthBactBact((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_WIDTH_BACT_BEAD:")) setWellWidthBactBead((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_DEPTH_BACT_BEAD:")) setWellDepthBactBead((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_WIDTH_BEAD_BEAD:")) setWellWidthBeadBead((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_DEPTH_BEAD_BEAD:")) setWellDepthBeadBead((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_WIDTH_BEAD_BDRY:")) setWellWidthBeadBdry((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_DEPTH_BEAD_BDRY:")) setWellDepthBeadBdry((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_WIDTH_BACT_BDRY:")) setWellWidthBactBdry((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_DEPTH_BACT_BDRY:")) setWellDepthBactBdry((double)args[0]);
		
		else if(line[0].equals("PHYSICS_WELL_WIDTH_VES_BDRY:")) setWellWidthVesBdry((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_DEPTH_VES_BDRY:")) setWellDepthVesBdry((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_WIDTH_VES_BEAD:")) setWellWidthVesBead((double)args[0]);
		else if(line[0].equals("PHYSICS_WELL_DEPTH_VES_BEAD:")) setWellDepthVesBead((double)args[0]);
		
		else if(line[0].equals("PHYSICS_REACT_FORCE:")) setReactForce((double)args[0]);				
		
		else if(line[0].equals("FIELD_GOAL_DEFINE:")) setCfGoalDefine(args);
		else if(line[0].equals("FIELD_GOAL_SETUP:")) setCfGoalSetup(args);
		else if(line[0].equals("FIELD_COORD_DEFINE:")) setCfCoordDefine(args);
		else if(line[0].equals("FIELD_COORD_SETUP:")) setCfCoordSetup(args);
		else if(line[0].equals("FIELD_RECRUIT_DEFINE:")) setCfRecruitDefine(args);
		else if(line[0].equals("FIELD_RECRUIT_SETUP:")) setCfRecruitSetup(args);
		else if(line[0].equals("FIELD_QUORUM_DEFINE:")) setCfQuorumDefine(args);
		else if(line[0].equals("FIELD_QUORUM_SETUP:")) setCfQuorumSetup(args);
		
		else if(line[0].equals("BOUNDING_BOX_DEFINE:")) setBoundingBoxDefine(args);
		
		else if(line[0].equals("MAGN_FIELD_STRENGTH")) setMagnStrength(args);
		
		else if(line[0].equals("VISUAL_AID_BACTERIA_TRACE:")) addBacteriaTrace(args);
		else if(line[0].equals("VISUAL_AID_AVG_BACTERIA_TRACE:")) addAvgBacteriaTrace(args);
		else if(line[0].equals("VISUAL_AID_BEAD_TRACE:")) addBeadTrace(args);
		else if(line[0].equals("VISUAL_AID_CLOCK:")) addClock(args);
		else if(line[0].equals("VISUAL_AID_SCALE:")) addScale(args);
		
		else if(line[0].equals("VIDEO_FRAMES_SKIP:")) setVideoFramesSkip((int)args[0]);
		else if(line[0].equals("RECORD_VIDEO:")) {
			if((int)args[0] == 1) setRecordVideo(true);
			else setRecordVideo(false);
		}
		
		else if(line[0].equals("DATA_FRAMES_SKIP:")) setDataFramesSkip((int)args[0]);
		
		else if(line[0].equals("SIMULATION_LENGTH:")) setSimLength((int)args[0]);
		else if(line[0].equals("SIMULATION_RUNS:")) setSimRuns((int)args[0]);
		
		else if(line[0].equals("EXPORT_DIR:")) setExportDir(line[1]);
		else if(line[0].equals("NUMBER_OF_THREADS:")) setNumOfThreads((int)args[0]);
		
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
		
	/**
	 * Standard get methods.
	 */
	public double 	getDtSecs() {return dt;}	
	public double	getBactRadius() {return bactRadius;}	
	public double 	getBeadRadius() {return beadRadius;}
	public double 	getBactForceUp() {return bactForceUp;}
	public double 	getBactForceDown() {return bactForceDown;}
	public double 	getBactSpeed() {return bactSpeed;}
	public double 	getViscosity() {return visc;}
	public int 		getScreenHeight() {return screenHeight;}
	public int 		getScreenWidth() {return screenWidth;}
	public double 	getUpRunLength() {return runLengthUp;}
	public double 	getDownRunLength() {return runLengthDown;}
	public double 	getIsoRunLength() {return runLengthIso;}
	public Vector 	getSingleBacteria() { return bacteriaSingles; }
	public Vector 	getSingleBead() { return beadSingles; }
	public Vector 	getBacteriaSet() { return bacteriaSets; }
	public Vector 	getBeadSet() { return beadSets; }	
	
	public Vector   getBacteriaTraces() { return vaBacteriaTraces; }
	public Vector   getAvgBacteriaTraces() { return vaAvgBacteriaTraces; }
	public Vector   getBeadTraces() { return vaBeadTraces; }
	public Vector   getClocks() { return vaClocks; }
	public Vector	getScales() { return vaScales; }
	
	public double	getScreenZoom() {return screenZoom;}
	
	//parameters needed to control the processing camera
	public double   getMinimumDistance() {return minimumDistance;}
	public double   getMaximumDistance() {return maximumDistance;}
	public double   getDefaultDistance() {return defaultDistance;}
	public int      getFrameForSec() {return frameForSec;}
	public int      getFrameRecordForSec() {return frameRecordForSec;}
	
	public double[]	getScreenMove() {return screenMove;}
	public int		getDataFramesSkip() {return dataFramesSkip;}
	public int		getVideoFramesSkip() {return videoFramesSkip;}
	public int		getSimRuns() {return simRuns;}
	public int		getSimLength() {return simLength;}
	public boolean  getRecordVideo() {return recordVideo;}
	public String	getExportDir() {return exportDir;}
	public int		getNumOfThreads() {return numOfThreads;}
	
	public double   getWellWidthBactBact () { return wellWidthBactBact; }
	public double   getWellWidthBactBead () { return wellWidthBactBead; }
	public double   getWellWidthBeadBead () { return wellWidthBeadBead; }
	public double   getWellWidthBeadBdry () { return wellWidthBeadBdry; }
	public double   getWellWidthBactBdry () { return wellWidthBactBdry; }
	
	public double   getWellDepthBactBact () { return wellDepthBactBact; }
	public double   getWellDepthBactBead () { return wellDepthBactBead; }
	public double   getWellDepthBeadBead () { return wellDepthBeadBead; }
	public double   getWellDepthBeadBdry () { return wellDepthBeadBdry; }
	public double   getWellDepthBactBdry () { return wellDepthBactBdry; }
	
	public double   getWellWidthVesBdry () { return wellWidthVesBdry; }
	public double   getWellDepthVesBdry () { return wellDepthVesBdry; }
	public double   getWellWidthVesBead () { return wellWidthVesBead; }
	public double   getWellDepthVesBead () { return wellDepthVesBead; }
	
	
	public double	getReactForce () { return reactForce; }
	
	public double[]	 getMagnStrength () {return magnStrength; }
	
	public double[]  getCfGoalDefine () { return cfGoalDefine; }
	public double[]  getCfGoalSetup () { return cfGoalSetup; }
	public double[]  getCfCoordDefine () { return cfCoordDefine; }
	public double[]  getCfCoordSetup () { return cfCoordSetup; }
	public double[]  getCfRecruitDefine () { return cfRecruitDefine; }
	public double[]  getCfRecruitSetup () { return cfRecruitSetup; }
	
	public double[]  getBoundingBoxDefine( double[] x) { return boundingBoxDefine; }
	
	
	/**
	 * Standard set methods.
	 */
	public void 	setDtSecs(double newDt) {dt = newDt;}	
	public void		setBactRadius(double r) {bactRadius = r;}	
	public void 	setBeadRadius(double r) {beadRadius = r;}
	public void 	setBactForceUp(double f) {bactForceUp = f;}
	public void 	setBactForceDown(double f) {bactForceDown = f;}
	public void 	setBactSpeed(double s) {bactSpeed = s;}
	public void 	setViscosity(double v) {visc = v;}
	public void 	setScreenHeight(int h) {screenHeight =h;}
	public void 	setScreenWidth(int w) {screenWidth = w;}
	public void 	setUpRunLength(double l) {runLengthUp = l;}
	public void 	setDownRunLength(double l) {runLengthDown = l;}
	public void 	setIsoRunLength(double l) {runLengthIso = l;}

	public void 	addSingleBacterium(double[] b) { bacteriaSingles.add(b); }
	public void 	addSingleBead(double[] p) { beadSingles.add(p); }
	public void 	addBacteriaSet(double[] b) { bacteriaSets.add(b); }
	public void 	addBeadSet(double[] p) { beadSets.add(p); }
		
	public void 	addBacteriaTrace(double[] p) { vaBacteriaTraces.add(p); }
	public void 	addAvgBacteriaTrace(double[] p) { vaAvgBacteriaTraces.add(p); }
	public void 	addBeadTrace(double[] p) { vaBeadTraces.add(p); }
	public void 	addClock(double[] p) { vaClocks.add(p); }
	public void		addScale(double[] p) { vaScales.add(p); }
	
	public void		setScreenZoom(double z) { screenZoom = z;}
	public void		setScreenMove(double x, double y) { screenMove[0] = x; screenMove[1] = y; }
	public void		setDataFramesSkip(int x) { dataFramesSkip = x;}
	public void	 	setVideoFramesSkip(int x) { videoFramesSkip = x;}
	public void	 	setSimRuns(int x) { simRuns = x;}
	public void	 	setSimLength(int x) { simLength = x;}
	public void  	setRecordVideo(boolean x) { recordVideo = x;}
	public void	 	setExportDir(String x) { exportDir = x;}
	public void 	setNumOfThreads(int x) { numOfThreads = x;}
	
	//parameters needed to control the processing camera
	public void     setMinimumDistance( double x) {minimumDistance = x;}
	public void     setMaximumDistance( double x) {maximumDistance = x;}
	public void     setDefaultDistance( double x) {defaultDistance = x;}
	public void     setFrameForSec(int x) {frameForSec = x;}
	public void     setFrameRecordForSec(int x) {frameRecordForSec = x;}
	
	public void     setWellWidthBactBact ( double x) { wellWidthBactBact = x; }
	public void     setWellWidthBactBead ( double x) { wellWidthBactBead = x; }
	public void     setWellWidthBeadBead ( double x) { wellWidthBeadBead = x; }
	public void     setWellWidthBeadBdry ( double x) { wellWidthBeadBdry = x; }
	public void     setWellWidthBactBdry ( double x) { wellWidthBactBdry = x; }
	
	public void     setWellDepthBactBact ( double x) { wellDepthBactBact = x; }
	public void     setWellDepthBactBead ( double x) { wellDepthBactBead = x; }
	public void     setWellDepthBeadBead ( double x) { wellDepthBeadBead = x; }
	public void     setWellDepthBeadBdry ( double x) { wellDepthBeadBdry = x; }
	public void     setWellDepthBactBdry ( double x) { wellDepthBactBdry = x; }
	
	public void     setWellWidthVesBdry (double x) { wellWidthVesBdry = x; }
	public void     setWellDepthVesBdry (double x) { wellDepthVesBdry = x; }
	public void     setWellWidthVesBead (double x) { wellWidthVesBead = x; }
	public void     setWellDepthVesBead (double x) { wellDepthVesBead = x; }
	
	public void 	setReactForce ( double x) { reactForce = x; }
	
	public void		setMagnStrength (double[] x) {magnStrength = x; }
	
	public void     setCfGoalDefine ( double[] x) { cfGoalDefine = x; }
	public void     setCfGoalSetup ( double[] x) { cfGoalSetup = x; }
	public void     setCfCoordDefine ( double[] x) { cfCoordDefine = x; }
	public void     setCfCoordSetup ( double[] x) { cfCoordSetup = x; }
	public void     setCfRecruitDefine ( double[] x) { cfRecruitDefine = x; }
	public void     setCfRecruitSetup ( double[] x) { cfRecruitSetup = x; }
	public void     setCfQuorumDefine ( double[] x) { cfQuorumDefine = x; }
	public void     setCfQuorumSetup ( double[] x) { cfQuorumSetup = x; }
	
	public void     setBoundingBoxDefine( double[] x) { boundingBoxDefine = x; }

	
	/**
	 * Create methods to create the sets of items that are required by the scene.
	 */
	
	public Vector createNewBacteriaVec(BSimScene scene) {
		int i;
		
		// Vector to hold the new objects
		Vector newVec = new Vector();
		
		// Create a new single bead for every item in the list
		for(i=0; i<bacteriaSingles.size(); i++){
			newVec.add(BSimBacteriaCreate.createBacterium((double[])bacteriaSingles.elementAt(i), scene, this));
		}
		
		// Create a new bead set for every item in the list
		for(i=0; i<bacteriaSets.size(); i++){
			newVec.addAll(BSimBacteriaCreate.createBacteriaSet((double[])bacteriaSets.elementAt(i), scene, this));
		}
		
		// Return the new vector
		return newVec;
	}
	
	public Vector createNewBeadVec() {
		int i;
		
		// Vector to hold the new objects
		Vector newVec = new Vector();
		
		// Create a new single bead for every item in the list
		for(i=0; i<beadSingles.size(); i++){
			newVec.add(BSimBeadsCreate.createBead((double[])beadSingles.elementAt(i), this));
		}
		
		// Create a new bead set for every item in the list
		for(i=0; i<beadSets.size(); i++){
			newVec.addAll(BSimBeadsCreate.createBeadSet((double[])beadSets.elementAt(i), this));
		}
		
		// Return the new vector
		return newVec;
	}
					
	public Vector createNewVisualAidsVec(BSimScene scene) {
		int i;
		
		// Vector to hold the new objects
		Vector newVec = new Vector();
		
		// Loop through each of the visual aid types, create them and add to the same vector
		for(i=0; i<vaBacteriaTraces.size(); i++){
			newVec.add(BSimVisualAidCreate.createBacteriaTrace(scene, (double[])vaBacteriaTraces.elementAt(i)));
		}
		for(i=0; i<vaAvgBacteriaTraces.size(); i++){
			newVec.add(BSimVisualAidCreate.createAvgBacteriaTrace(scene, (double[])vaAvgBacteriaTraces.elementAt(i)));
		}
		for(i=0; i<vaBeadTraces.size(); i++){
			newVec.add(BSimVisualAidCreate.createBeadTrace(scene, (double[])vaBeadTraces.elementAt(i)));
		}
		for(i=0; i<vaClocks.size(); i++){
			newVec.add(BSimVisualAidCreate.createSceneClock(scene));
		}
		for(i=0; i<vaScales.size(); i++){
			newVec.add(BSimVisualAidCreate.createSceneScale(scene, (double[])vaScales.elementAt(i)));
		}
		
		// Return the new vector
		return newVec;
	}
	
	public BSimChemicalField createNewGoalChemicalField() {
		
		// Create the new chemical field
		return BSimChemicalFieldCreate.createChemicalField (cfGoalDefine, cfGoalSetup,
		                                               new Color(0.8f, 0.1f, 0.1f), this);
	}
	
	public BSimChemicalField createNewCoordChemicalField() {
		
		// Create the new chemical field
		return BSimChemicalFieldCreate.createChemicalField (cfCoordDefine, cfCoordSetup, 
		                                               new Color(0.1f, 0.1f, 0.8f), this);
	}
	
	public BSimChemicalField createNewRecruitChemicalField() {
		
		// Create the new chemical field
		return BSimChemicalFieldCreate.createChemicalField (cfRecruitDefine, cfRecruitSetup, 
		                                               new Color(0.1f, 0.8f, 0.1f), this);
	}
	
	public BSimChemicalField createNewQuorumChemicalField() {
		
		// Create the new chemical field
		return BSimChemicalFieldCreate.createChemicalField (cfQuorumDefine, cfQuorumSetup, 
		                                               new Color(0.1f, 0.8f, 0.1f), this);
	}
	
	public BSimBoundingBox createNewBoundingBox() {
		
		// Create the new bounding box
		return new BSimBoundingBox(boundingBoxDefine);
	}
}
