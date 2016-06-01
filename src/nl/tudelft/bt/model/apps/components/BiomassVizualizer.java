package nl.tudelft.bt.model.apps.components;

import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.exceptions.InvalidValueException;
import nl.tudelft.bt.model.gui.Biomass2DVisualizerPane;

/**
 * The simple vizualization interface to an application
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class BiomassVizualizer extends ApplicationDecorator {
	private Biomass2DVisualizerPane _vizualizer;
	
	/**
	 * Creates a new biomass vizualizer window
	 * 
	 * @param c
	 */
	public BiomassVizualizer(ApplicationComponent c) {
		super(c);
		_vizualizer = new Biomass2DVisualizerPane(this);
	}

	/* (non-Javadoc)
	 * @see org.photobiofilms.phlip.apps.components.ApplicationDecorator#output()
	 */
	protected void output() {
		_vizualizer.update();
	}
	
	/**
	 *  This method is overriden to allow for (re)defining system size 
	 * parameters in vizualizer window
	 */
	public void initializeSystemSpace() throws InvalidValueException {
		super.initializeSystemSpace();
		_vizualizer.setSystemSize();
	}	
}
