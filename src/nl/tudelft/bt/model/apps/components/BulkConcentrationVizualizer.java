package nl.tudelft.bt.model.apps.components;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.gui.BulkConcentrationVisualizerPane;
/**
 * The simple vizualization interface to an application
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class BulkConcentrationVizualizer extends ApplicationDecorator {
	private BulkConcentrationVisualizerPane _vizualizer;
	/**
	 * Creates a new concentration plot window
	 * 
	 * @param c
	 */
	public BulkConcentrationVizualizer(ApplicationComponent c) {
		super(c);
		_vizualizer = new BulkConcentrationVisualizerPane(this);
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.components.ApplicationDecorator#output()
	 */
	protected void output() {
		_vizualizer.update();
	}
}
