/**
 * BSimToolbar.java
 *
 * Class to generate toolbar for the BSim user interface. This is the main method to 
 * interact with the simulation.
 *
 * Authors: Thomas Gorochowski
 * Created: 14/07/2008
 * Updated: 24/07/2008
 */
package bsim.app.gui;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JToolBar;

import bsim.app.BSimApp;
import bsim.scene.BSimScene;


public class BSimToolbar extends JToolBar implements ActionListener{

	// Other BSim references
	private BSimApp app;
	private BSimScene scene;
	
	// State variable for the playback controls
	private int playState = 0;
	
	// GUI objects
	private JButton btnPlayPause, btnReset, btnSaveRecord, btnStartRecord, btnEndRecord, btnSaveScreenshot, btnTakeScreenshot;
	private JFileChooser fc;
	
	// Images used on the buttons
	private static final ImageIcon iconPlay = new ImageIcon(BSimToolbar.class.getResource("../../resource/icons/play.png"));
	private static final ImageIcon iconPause = new ImageIcon(BSimToolbar.class.getResource("../../resource/icons/pause.png"));
	private static final ImageIcon iconReset = new ImageIcon(BSimToolbar.class.getResource("../../resource/icons/reset.png"));
	
	private static final ImageIcon iconSaveRecord = new ImageIcon(BSimToolbar.class.getResource("../../resource/icons/saveRecord.png"));
	private static final ImageIcon iconStartRecord = new ImageIcon(BSimToolbar.class.getResource("../../resource/icons/startRecord.png"));
	private static final ImageIcon iconEndRecord = new ImageIcon(BSimToolbar.class.getResource("../../resource/icons/endRecord.png"));
	private static final ImageIcon iconStartRecordDisabled = new ImageIcon(BSimToolbar.class.getResource("../../resource/icons/startRecordDisabled.png"));
	private static final ImageIcon iconEndRecordDisabled = new ImageIcon(BSimToolbar.class.getResource("../../resource/icons/endRecordDisabled.png"));
	
	private static final ImageIcon iconSaveScreenshot = new ImageIcon(BSimToolbar.class.getResource("../../resource/icons/saveScreenshot.png"));
	private static final ImageIcon iconTakeScreenshot = new ImageIcon(BSimToolbar.class.getResource("../../resource/icons/takeScreenshot.png"));
	private static final ImageIcon iconTakeScreenshotDisabled = new ImageIcon(BSimToolbar.class.getResource("../../resource/icons/takeScreenshotDisabled.png"));	
	
	
	/**
	 * Creates a new toolbar for a given BSimApp and BSimScene.
	 */
	public BSimToolbar(BSimApp newApp){
		super();
		// Update internal variables
		app = newApp;
		scene = app.getScene();
		
		// Create the toolbar
		setupToolBar();
		setOrientation(HORIZONTAL);
		setFloatable(false);		
	}
	
	
	/**
	 * Create the toolbar.
	 */
	private void setupToolBar(){
		// Create the default file chooser
		fc = new JFileChooser();
		
		// Create GUI controls with initial properties
		btnPlayPause = new JButton("Play");
		btnPlayPause.setIcon(iconPlay);
		btnPlayPause.addActionListener(this);
		
		btnReset = new JButton("Reset");
		btnReset.setIcon(iconReset);
		btnReset.addActionListener(this);

		// Create a folder with the date & time: export_01_09_09_1244\screenshots, .\movies, .\data for example
		btnSaveScreenshot = new JButton("Screenshot");
		btnSaveScreenshot.setIcon(iconSaveScreenshot);
		btnSaveScreenshot.addActionListener(this);
		btnTakeScreenshot = new JButton("Save Screenshot");
		btnTakeScreenshot.setIcon(iconTakeScreenshotDisabled);
		btnTakeScreenshot.addActionListener(this);
		
		btnSaveRecord = new JButton("Record");
		btnSaveRecord.setIcon(iconSaveRecord);
		btnSaveRecord.addActionListener(this);
		btnStartRecord = new JButton("Start Record");
		btnStartRecord.setIcon(iconStartRecordDisabled);
		btnStartRecord.addActionListener(this);
		btnEndRecord = new JButton("End Record");
		btnEndRecord.setIcon(iconEndRecordDisabled);
		btnEndRecord.addActionListener(this);
				
		// Add objects to the toolbar
		// Playback controls
		this.add(btnPlayPause);
		this.add(btnReset);
		this.addSeparator();
		
		// Screenshot controls
		this.add(btnSaveScreenshot);
		this.add(btnTakeScreenshot);
		btnTakeScreenshot.setEnabled(false);
		this.addSeparator();
		// Movie record controls
		this.add(btnSaveRecord);
		this.add(btnStartRecord);
		this.add(btnEndRecord);	
				
		btnStartRecord.setEnabled(false);
		btnEndRecord.setEnabled(false);
		btnPlayPause.setEnabled(true);
		btnReset.setEnabled(false);		
		btnSaveScreenshot.setEnabled(false);
		btnSaveRecord.setEnabled(false);
	}

	
	/**
	 * Handle events for the toolbar.
	 */
	public void actionPerformed(ActionEvent e){
		// Get the object that was clicked
		Object s = e.getSource();
		
		// Play button - switch between play and pause
		if(s == btnPlayPause){
			if(playState == 0){
				playState = 1;
				btnPlayPause.setIcon(iconPause);
				btnPlayPause.setText("Pause");
				app.play();
			}
			else{
				playState = 0;
				btnPlayPause.setIcon(iconPlay);
				btnPlayPause.setText("Play");
				app.pause();
			}
			
		// Reset button
		}else if(s == btnReset){
			// Reset the current simulation
			playState = 0;
			btnPlayPause.setIcon(iconPlay);
			btnPlayPause.setText("Play");
			app.reset(0);
			
		//Record button
		}else if (e.getSource() == btnSaveRecord) {		
			// Variable to check if cancel is pressed (also displays the file dialog)
			int returnVal = fc.showSaveDialog(this);
				
			// If OK is pressed
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				// Get the file that has been entered
				File file = fc.getSelectedFile();
				String videoFileName = file.getPath();
				videoFileName = videoFileName + ".mov";
				file.delete();
				file=null;
				btnSaveRecord.setEnabled(false);
				btnStartRecord.setEnabled(true);
				btnStartRecord.setIcon(iconStartRecord);
				//scene.getProcessing().createMovie(videoFileName);
			}
		}else if (e.getSource() == btnStartRecord) {		
			btnStartRecord.setEnabled(false);
			btnStartRecord.setIcon(iconStartRecordDisabled);
			btnEndRecord.setEnabled(true);
			btnEndRecord.setIcon(iconEndRecord);
			scene.setStartVideo(true);
		}else if (e.getSource() == btnEndRecord) {		
			btnEndRecord.setEnabled(false);
			btnEndRecord.setIcon(iconEndRecordDisabled);
			scene.setStartVideo(false);
			scene.setEndVideo(true);
			btnSaveRecord.setEnabled(true);
			
		//Screenshot button
		}else if(s == btnSaveScreenshot){
			int returnVal = fc.showSaveDialog(this);
				
			// If OK is pressed
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				// Get the file that has been entered
				File file = fc.getSelectedFile();
				String imageFileName = file.getPath();
				imageFileName = imageFileName + ".png";
				file.delete();
				file=null;
				scene.setImageFileName(imageFileName);
				btnTakeScreenshot.setEnabled(true);
				btnTakeScreenshot.setIcon(iconTakeScreenshot);
				btnSaveScreenshot.setEnabled(false);
			}
		}else if(s == btnTakeScreenshot){
			//scene.getProcessing().takeImage(scene.getImageFileName());
			btnTakeScreenshot.setEnabled(false);
			btnTakeScreenshot.setIcon(iconTakeScreenshotDisabled);
			btnSaveScreenshot.setEnabled(true);	
			
		//Load simulation
		}
	}	
}
