package nl.tudelft.bt.model.particlebased.granule.phbgranule;

import javax.swing.JApplet;

/**
 * @author jxavier
 */
public class Monospecies2DGranuleApplet extends JApplet {
	/**
	 * Constructs the demo applet.
	 */
	public Monospecies2DGranuleApplet() {
		Monospecies2DGranuleGUI gui = new Monospecies2DGranuleGUI();
		getContentPane().add(gui);

	}

}