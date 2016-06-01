package nl.tudelft.bt.model.apps.components;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.output.VariableSeries;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.gui.VariableSeriesVisualizerPane;
/**
 * The simple vizualization interface to an application
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DetachedBiomassVizualizer extends ApplicationDecorator {
	private VariableSeriesVisualizerPane _vizualizer;
	/**
	 * Creates a new detached biomass plot window
	 * 
	 * @param c
	 */
	public DetachedBiomassVizualizer(ApplicationComponent c) {
		super(c);
		VariableSeries[] vs = {
				Model.model().detachedBiomassContainer()
						.getTotalDetachedBiomassSeries(),
				Model.model().detachedBiomassContainer()
						.getErodedBiomassSeries(),
				Model.model().detachedBiomassContainer()
						.getSloughedBiomassSeries()};
		_vizualizer = new VariableSeriesVisualizerPane(this, vs,
				"Detached biomass");
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.components.ApplicationDecorator#output()
	 */
	protected void outputVizualizer() {
		_vizualizer.update();
	}
	protected void output() {
		//do nothig here, as output should be performed only after detachment
		// is evoked
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.ApplicationComponent#detach()
	 */
	public void detach() throws ModelException {
		super.detach();
		outputVizualizer();
	}
}