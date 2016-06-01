package nl.tudelft.bt.model.apps.components;

import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.exceptions.InvalidValueException;
import nl.tudelft.bt.model.gui.Visualizer2DControler;

/**
 * The simple vizualization interface to an application
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class VizualModelControler extends ApplicationDecorator {
	private Visualizer2DControler _vizualizer;
	
	/**
	 * Creates a new biomass vizualizer window
	 * 
	 * @param c
	 */
	public VizualModelControler(ApplicationComponent c) {
		super(c);
		_vizualizer = new Visualizer2DControler();
		_vizualizer.createFrame(this);	
	}

	/* (non-Javadoc)
	 * @see org.photobiofilms.phlip.apps.components.ApplicationDecorator#output()
	 */
	protected void output() {
	}
	
	/**
	 *  This method is overriden to allow for (re)defining system size 
	 * parameters in vizualizer window
	 */
	public void initializeSystemSpace() throws InvalidValueException {
		super.initializeSystemSpace();
	}
}
