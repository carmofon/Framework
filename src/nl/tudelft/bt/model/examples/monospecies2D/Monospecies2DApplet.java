package nl.tudelft.bt.model.examples.monospecies2D;

import javax.swing.JApplet;

/**
 * Applet interface for Monospecies2DGUI
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Monospecies2DApplet extends JApplet {
	/**
	 * Constructs the applet.
	 */
	public Monospecies2DApplet() {
		Monospecies2DGUI gui = new Monospecies2DGUI();
		getContentPane().add(gui);
	}
}