/*
 * Created on 26-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.examplesroland;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.bulkconcentrations.BulkConcentration;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;

/**
 * Implements a list bulk concentration manager. The values of the bulk
 * concentration are read from a file
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ListBulkConcentration extends BulkConcentration {
	private float[] _times = { 0, 1, 10.0f, 20.f };

	private float[] _concentrations = { 45f, 50f, 60f, 70f };

	private float _smallTimeInterval;

	/**
	 * @param v
	 *            value for constant bulk concentration
	 */
	public ListBulkConcentration(String fileName) {
		super();
		// check consistency of arrays
		if (_times.length != _concentrations.length)
			throw new ModelRuntimeException(
					"Time and concentration arrays do not have the same length");
		float previousTimeValue = _times[0];
		for (int i = 1; i < _times.length; i++) {
			if (_times[i] <= previousTimeValue)
				throw new ModelRuntimeException(
						"Time series in list concentration is not increasing");
			previousTimeValue = _times[i];
		}
		// set the value at time 0
		setValue(_concentrations[0]);
	}

	/**
	 * Computes and returns the value for the bulk concetration
	 * 
	 * @param tstep
	 *            time step
	 */
	public void computeBulkConcentration(float tstep) {
		// check if it is time to change
		int i;
		for (i = 0; i < _times.length; i++) {
			if (_times[i] > Model.model().getTime()) {
				setValue(_concentrations[i - 1]);
				return;
			}
		}
		setValue(_concentrations[i - 1]);
	}

	/**
	 * Ensures that a time point right before the change and another right after
	 * are computed
	 * 
	 */
	public TimeStepConstraint getMaximumTimeStep() {
		TimeStepConstraint constraint = super.getMaximumTimeStep();
		// Determine the time until next regime change
		float t = Model.model().getTime();
		float timeToChange = 0;
		// determine the time until the next change
		int i;
		float nextChange = 0;
		for (i = 0; i < _times.length; i++) {
			if (_times[i] > t) {
				nextChange = _times[i];
				break;
			}
		}
		if (nextChange == 0)
			return constraint;
		timeToChange = _times[i] - t;
		// Compute a point right before change to e
		// ensure that step change shows
		// in plot
		_smallTimeInterval = 1e-5f * (t == 0 ? timeToChange : t);
		if (timeToChange <= 2 * _smallTimeInterval) {
			constraint.setTimeStep(_smallTimeInterval);
			constraint.setName("List species (" + _species.getName()
					+ ") in transition");
			return constraint;
		}
		// Add the predefined small value
		constraint.setTimeStep(timeToChange - _smallTimeInterval);
		constraint.setName("List species (" + _species.getName()
				+ ") before transition");
		return constraint;
	}
}
