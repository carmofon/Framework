package nl.tudelft.bt.model.apps.output;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;

public class DetachedParticulateSpecies extends MassPerAreaTimeSeries {
	private ParticulateSpecies _p;

	private DetachedSeries _totalBiomassDetachedSeries;

	public DetachedParticulateSpecies(ParticulateSpecies p) {
		super("Detached " + p.getName());
		_p = p;
		_totalBiomassDetachedSeries = (DetachedSeries) Model.model()
				.detachedBiomassContainer().getTotalDetachedBiomassSeries();
	}

	/*
	 * This method updates the Y array each time that last Y is called and X is
	 * longer than Y
	 */
	public float getLastY() {
		float d = _p.getDetachedBiomass();
		_p.resetDetachedBiomass();
		return d;
	}
}
