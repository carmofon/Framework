package nl.tudelft.bt.model.apps.components;
import java.util.Collection;
import java.util.Iterator;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.output.VariableSeries;
import nl.tudelft.bt.model.gui.VariableSeriesVisualizerPane;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
/**
 * Asimple vizualization interface to an application
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SoluteRateSeriesVizualizer extends ApplicationDecorator {
	private VariableSeriesVisualizerPane _vizualizer;
	/**
	 * Creates a new plot window
	 * 
	 * @param c
	 */
	public SoluteRateSeriesVizualizer(ApplicationComponent c) {
		super(c);
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.components.ApplicationDecorator#output()
	 */
	protected void output() {
		// When evocked for the first time, initialize the vizualizer with the
		// time series for all the solute sepcies rates
		if (_vizualizer == null) {
			Collection ss = Model.model().getSoluteSpecies();
			VariableSeries[] vs = new VariableSeries[ss.size()];
			int i = 0;
			//get the rate for each solute species
			for (Iterator iter = ss.iterator(); iter.hasNext();) {
				SoluteSpecies s = (SoluteSpecies) iter.next();
				vs[i++] = s.getRateTimeSeries();
			}
			_vizualizer = new VariableSeriesVisualizerPane(this, vs,
					"Solutes Rate");
		}
		_vizualizer.update();
	}
}
