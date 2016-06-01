/*
 * Created on 26-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.bulkconcentrations;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;

/**
 * Implements a constant bulk concentration that is 0 at all times except in the
 * duration of a single pulse
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SinglePulseBulkConcentration extends BulkConcentration {
	private float _initialTimePulse;
	private float _endPulse;
	private float _concentrationPulse;
	private float _smallTimeInterval;

	/**
	 * @param value
	 *            of concentration during pulse
	 * @param begining
	 *            time for begining of pulse
	 * @param duration
	 *            duration of pulse
	 */
	public SinglePulseBulkConcentration(float value, float begining,
			float duration) {
		super();
		_initialTimePulse = begining;
		_endPulse = begining + duration;
		_concentrationPulse = value;
		_smallTimeInterval = 1e-3f * Math.min(_endPulse, _initialTimePulse);
		setValue();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.bulkconcentrations.BulkConcentration#computeBulkConcentration(float)
	 */
	public void computeBulkConcentration(float tstep) {
		setValue();
	}

	/**
	 * Sets the value of bulk concnetration to _concentrationPulse if time is
	 * within the pulse time
	 */
	private void setValue() {
		float time = Model.model().getTime();
		if ((time >= _initialTimePulse) & (time <= _endPulse)) {
			setValue(_concentrationPulse);
		} else
			setValue(0);
	}

	/**
	 * Ensures that a time point right before each change and another right
	 * after are computed
	 */
	public TimeStepConstraint getMaximumTimeStep() {
		TimeStepConstraint constraint = super.getMaximumTimeStep();
		//Determine the time until next regime change
		float time = Model.model().getTime();
		if (time < _initialTimePulse) {
			// before pulse begining
			constraint.setTimeStep(_initialTimePulse - time
					+ _smallTimeInterval);
			constraint.setName("Single pulse (" + _species.getName()
					+ ") before pulse");
			return constraint;
		} else if (time < _endPulse) {
			// during pulse
			float timeToChange = _endPulse - time;
			float timeSinceChange = time - _initialTimePulse;
			if (timeSinceChange < 2 * _smallTimeInterval) {
				// before pulse begining
				constraint.setTimeStep(_smallTimeInterval);
				constraint.setName("Single pulse (" + _species.getName()
						+ ") during pulse");
				return constraint;
			}
			// Add the predefined small value
			// before pulse begining
			constraint.setTimeStep(timeToChange + _smallTimeInterval);
			constraint.setName("Single pulse (" + _species.getName()
					+ ") time until change");
			return constraint;
		}
		//if only small time has passed since end of pulse, make another pass
		// at timestep small
		float timeSinceChange = time - _endPulse;
		if (timeSinceChange < 2 * _smallTimeInterval) {
			constraint.setTimeStep(_smallTimeInterval);
			constraint.setName("Single pulse (" + _species.getName()
					+ ") right after change");
			return constraint;
		}
		// If pulse has already passed lon ago, no restrictions to the time step
		// are applied
		constraint.setName("Single pulse (" + _species.getName()
				+ ") already after");		
		return constraint;
	}
}