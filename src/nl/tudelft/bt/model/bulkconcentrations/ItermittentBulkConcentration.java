/*
 * Created on 26-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.bulkconcentrations;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;

/**
 * Implements a constant bulk concentration With itermitant cyclic values. A
 * default value for concentration is set uppon construction. Time intervals for
 * feast and famine periods and length of an initial feast period are also
 * defined. During the feast periods value of concentration is the default
 * value. After the initial feast period, regime alternates between famine and
 * feat periods. Concentration in famine periods is 0;
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ItermittentBulkConcentration extends BulkConcentration {
	private float _initialFeastTime;

	private float _feastTime;

	private float _famineTime;

	private float _feastConcentration;

	private float _famineConcentration;

	private float _smallTimeInterval;

	/**
	 * Feast/famine cycles with concentration changing between vfeast and 0
	 * 
	 * @param v
	 *            feast concentration
	 * @param initalFeastTime
	 *            duration of an initial period of feast
	 * @param feastTime
	 *            feast period duration
	 * @param famineTime
	 *            famine period duration
	 */
	public ItermittentBulkConcentration(float v, float initalFeastTime,
			float feastTime, float famineTime) {
		super();
		initializeParameters(v, 0, initalFeastTime, feastTime, famineTime);
	}

	/**
	 * Feast/famine cycles with concentration changing between vfeast and
	 * vfamine
	 * 
	 * @param vfeast
	 *            feast concentration
	 * @param vfamine
	 *            famine concentration
	 * @param initalFeastTime
	 *            duration of an initial period of feast
	 * @param feastTime
	 *            feast period duration
	 * @param famineTime
	 *            famine period duration
	 */
	public ItermittentBulkConcentration(float vfeast, float vfamine,
			float initalFeastTime, float feastTime, float famineTime) {
		super();
		initializeParameters(vfeast, vfamine, initalFeastTime, feastTime,
				famineTime);
	}

	private void initializeParameters(float vfeast, float vfamine,
			float initalFeastTime, float feastTime, float famineTime) {
		_initialFeastTime = initalFeastTime;
		_feastTime = feastTime;
		_famineTime = famineTime;
		_feastConcentration = vfeast;
		_famineConcentration = vfamine;
		setValue(vfeast);
		// small time interval for restriction to time step
		// set to 1/1000 of the minimum among famine and feast times
		_smallTimeInterval = 1e-3f * Math.min(_feastTime, _famineTime);
	}

	/**
	 * Computes and returns the value for the bulk concetration
	 * 
	 * @param tstep
	 *            time step
	 */
	public void computeBulkConcentration(float tstep) {
		float t = Model.model().getTime() - _initialFeastTime;
		if (t >= 0) {
			if ((t % (_feastTime + _famineTime)) >= _famineTime)
				setValue(_feastConcentration);
			else
				setValue(_famineConcentration);
		}
		// do nothing if t < 0 since that means initial feast period is
		// undergoing
	}

	/**
	 * Ensures that a time point right before the change and another right after
	 * are computed
	 *  
	 */
	public TimeStepConstraint getMaximumTimeStep() {
		TimeStepConstraint constraint = super.getMaximumTimeStep();
		//Determine the time until next regime change
		float t = Model.model().getTime() - _initialFeastTime;
		float timeToChange = 0;
		float timeSinceChange = 0;
		if (t < 0) {
			constraint.setTimeStep(-t + _smallTimeInterval);
			constraint.setName("Itermittent species (" + _species.getName()
					+ ") in transition");
			return constraint;
		} else {
			// The time already in the present cycle
			float cycleTime = (t % (_feastTime + _famineTime));
			if (cycleTime < _famineTime) {
				// In famine regime, compute time to end of famine period
				timeToChange = _famineTime - cycleTime;
				timeSinceChange = cycleTime;
			} else {
				// In feast regime, compute time to end of feast period
				timeToChange = (_feastTime + _famineTime) - cycleTime;
				timeSinceChange = cycleTime - _famineTime;
			}
		}
		// Compute a point right before change to e
		//ensure that step change shows
		// in plot
		if (timeToChange <= 2 * _smallTimeInterval) {
			constraint.setTimeStep(_smallTimeInterval);
			constraint.setName("Itermittent species (" + _species.getName()
					+ ") in transition");
			return constraint;

		}
		// Add the predefined small value
		constraint.setTimeStep(timeToChange - _smallTimeInterval);
		constraint.setName("Itermittent species (" + _species.getName()
				+ ") before transition");
		return constraint;
	}
}