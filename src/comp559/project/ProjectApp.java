// NAME: Zhiguo(Frank) Zhang
// ID: 260550226
package comp559.project;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import javax.swing.JButton;
import javax.swing.JPanel;

import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.EasyViewer;
import mintools.viewer.Interactor;
import mintools.viewer.SceneGraphNode;


/**
 * Provided code for particle system simulator.
 * This class provides the mouse interface for clicking and dragging particles, and the 
 * code to draw the system.  When the simulator is running system.advanceTime is called
 * to numerically integrate the system forward.
 *
 * @author kry
 */
public class ProjectApp implements SceneGraphNode, Interactor {
	
	private EasyViewer ev;
	
	private ParticleSystem system;
	
	/**
	 * Entry point for application
	 * @param args
	 */
	public static void main(String[] args) {
		new ProjectApp();        
	}
	
	/**
	 * Creates the application / scene instance
	 */
	public ProjectApp() {
		system = new ParticleSystem();
		system.integrator = symplecticEuler;
		ev = new EasyViewer( "COMP 559 W2018 - Project 3D Cloth Simulation", this, new Dimension(640,480), new Dimension(640,480) );
		ev.addInteractor(this);
		ev.trackBall.fovy.setValue(20.0);
	}
	
	@Override
	public void init(GLAutoDrawable drawable) {
		GL2 gl = drawable.getGL().getGL2();
		gl.glEnable( GL.GL_BLEND );
		gl.glBlendFunc( GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA );
		gl.glEnable( GL.GL_LINE_SMOOTH );
		gl.glEnable( GL2.GL_POINT_SMOOTH );
		system.init(drawable);
	}
	
	@Override
	public void display(GLAutoDrawable drawable) {
		
		if ( run.getValue() ) {   
//			if (!) run.setValue(false);;        
			system.advanceTime( stepsize.getValue());
		}
		
		system.display( drawable );
		
		EasyViewer.beginOverlay( drawable );
		String text = system.toString() + "\n" + 
				"h = " + stepsize.getValue() + "\n";
		EasyViewer.printTextLines( drawable, text );
		EasyViewer.endOverlay(drawable);    
		
		if ( run.getValue() || stepped ) {
			stepped = false;        
			if ( record.getValue() ) {
				// write the frame
				File file = new File( "stills/" + dumpName + format.format(nextFrameNum) + ".png" );                                             
				nextFrameNum++;
				file = new File(file.getAbsolutePath().trim());
				ev.snapshot(drawable, file);
			}
		}
	}
	
	private BooleanParameter record = new BooleanParameter( "record (press ENTER in canvas to toggle)", false );
	
	/** 
	 * boolean to signal that the system was stepped and that a 
	 * frame should be recorded if recording is enabled
	 */
	private boolean stepped = false;
	
	private String dumpName = "dump";
	
	private int nextFrameNum = 0;
	private NumberFormat format = new DecimalFormat("00000");
	private BooleanParameter run = new BooleanParameter( "simulate", false );
	private DoubleParameter stepsize = new DoubleParameter( "step size", 0.05, 1e-5, 1 );
	
	@Override
	public JPanel getControls() {
		VerticalFlowPanel vfp = new VerticalFlowPanel();
		JButton create = new JButton("create cloth system");
		vfp.add( create );
		create.addActionListener( new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				system.createSystem();
			}
		});
		
		vfp.add( record.getControls() );
		vfp.add( run.getControls() );
		vfp.add( stepsize.getSliderControls(true) );
		vfp.add( system.getControls() );
		
		return vfp.getPanel();
	}
	
	@Override
	public void attach(Component component) {
		
		component.addKeyListener( new KeyAdapter() {
			@Override
			public void keyPressed(KeyEvent e) {
				if ( e.getKeyCode() == KeyEvent.VK_SPACE ) {
					run.setValue( ! run.getValue() ); 
				} else if ( e.getKeyCode() == KeyEvent.VK_S ) {
					system.advanceTime( stepsize.getValue() );
					stepped = true;
				} else if ( e.getKeyCode() == KeyEvent.VK_R ) {
					system.resetParticles();                    
				} else if ( e.getKeyCode() == KeyEvent.VK_C ) {                   
					system.clearParticles();
				} else if ( e.getKeyCode() == KeyEvent.VK_1 ) {
					system.explicit.setValue(false); 
				} else if ( e.getKeyCode() == KeyEvent.VK_ESCAPE ) {
					ev.stop();
				} else if ( e.getKeyCode() == KeyEvent.VK_ENTER ) {
					record.setValue( ! record.getValue() );
				}
				if ( e.getKeyCode() != KeyEvent.VK_ESCAPE ) ev.redisplay();
			}
		} );
	}
	private SymplecticEuler symplecticEuler = new SymplecticEuler();
}
