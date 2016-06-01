package nl.tudelft.bt.model.apps.components;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.output.VariableSeries;
import nl.tudelft.bt.model.gui.VariableSeriesVisualizerPane;
/**
 * The simple vizualization interface to an application
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SeriesVizualizer extends ApplicationDecorator {
	private VariableSeriesVisualizerPane _vizualizer;
	/**
	 * Creates a new plot window
	 * 
	 * @param c
	 */
	public SeriesVizualizer(ApplicationComponent c, VariableSeries[] s,
			String title) {
		super(c);
		_vizualizer = new VariableSeriesVisualizerPane(this, s, title);
	}
	/**
	 * 
	 * 
	 * @param c
	 * @param s
	 */
	public SeriesVizualizer(ApplicationComponent c, VariableSeries s) {
		super(c);
		VariableSeries[] vs = {s};
		_vizualizer = new VariableSeriesVisualizerPane(this, vs, s.getTitle());
	}
	/* (non-Javadoc)
	 * @see org.photobiofilms.model.apps.components.ApplicationDecorator#output()
	 */
	protected void output() {
		_vizualizer.update();
	}
}
