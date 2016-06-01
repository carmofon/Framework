/*
 * Created on 26-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.bulkconcentrations;

import java.io.Serializable;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;

/**
 * Implements bulk concentration of a solute in a sequencing batch ractor. A
 * cycle time is defined (cycleDuration). At the begining of each cycle the bulk
 * concentration of the solute is set to the valus of feedConcentration. During
 * The cycle, concentration of the solute decreases (if solute is consumed) or
 * increases (if solute is produced).
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SbrBulkConcentration extends BulkConcentration {
	private SbrCycle _cycle;

	private float _feedConcentration;

	private float _feedFraction;

	private boolean _timeToFeed;

	private boolean _alreadyFed;

	private float _smallTimeInterval;

	private static float _maxFractionToDecrease;

	private static float DEFAULTMAXFRACTION = 0.95f;

	private float _precisionInConcentration;

	/**
	 * Holds the properties of the SBR cycle
	 * 
	 * @author jxavier
	 */
	public static class SbrCycle implements Serializable {
		private float _cycleDuration;

		private int _elapsedCycles;

		/**
		 * Creates a new cycle with duration cycleDuration
		 * 
		 * @param cycleDuration
		 *            the duration of an SBR cycle
		 */
		public SbrCycle(float cycleDuration) {
			_cycleDuration = cycleDuration;
			_elapsedCycles = 0;
		}

		public int elapsedCycles() {
			return (int) Math.floor(Model.model().getTime() / _cycleDuration);
		}

		public float getCycleDuration() {
			return _cycleDuration;
		}
	}

	/**
	 * SBR cycles with duration cycleDuration
	 * 
	 * @param cycleDuration
	 *            The duration of a SBR cycle
	 * @param feedConcentration
	 *            The solute concentration on the feed
	 * @param feedFraction
	 *            the fraction of the liquid in the reaction that is replaced by
	 *            feed liquid
	 * @param precisionInConcentration
	 *            the precision in the concnetration value (below this value,
	 *            concentration is assumed 0)
	 */
	public SbrBulkConcentration(float feedConcentration, float feedFraction,
			float precisionInConcentration, SbrCycle cycle) {
		this(feedConcentration, feedFraction, precisionInConcentration,
				DEFAULTMAXFRACTION, cycle);
		// assign
		_feedConcentration = feedConcentration;
		_feedFraction = feedFraction;
		_precisionInConcentration = precisionInConcentration;
		_cycle = cycle;
		// initialize
		_timeToFeed = false;
		_alreadyFed = true;
		setValue(feedConcentration);
		// small time interval for restriction to time step
		// set to 1/1000 of the minimum among famine and feast times
		_smallTimeInterval = 0.5e-3f * cycle._cycleDuration;
	}

	/**
	 * SBR cycles with duration cycleDuration
	 * 
	 * @param cycleDuration
	 *            The duration of a SBR cycle
	 * @param feedConcentration
	 *            The solute concentration on the feed
	 * @param feedFraction
	 *            the fraction of the liquid in the reaction that is replaced by
	 *            feed liquid
	 * @param precisionInConcentration
	 *            the precision in the concnetration value (below this value,
	 *            concentration is assumed 0)
	 */
	public SbrBulkConcentration(float feedConcentration, float feedFraction,
			float precisionInConcentration, float maxFractionToDecrease,
			SbrCycle cycle) {
		super();
		// assign
		_feedConcentration = feedConcentration;
		_feedFraction = feedFraction;
		_precisionInConcentration = precisionInConcentration;
		_maxFractionToDecrease = maxFractionToDecrease;
		_cycle = cycle;
		// initialize
		_timeToFeed = false;
		_alreadyFed = true;
		setValue(feedConcentration);
		// small time interval for restriction to time step
		// set to 1/1000 of the minimum among famine and feast times
		_smallTimeInterval = 0.5e-3f * cycle._cycleDuration;
	}

	/**
	 * Computes and returns the value for the bulk concetration
	 * 
	 * @param tstep
	 *            time step
	 */
	public void computeBulkConcentration(float tstep) {
		float t = Model.model().getTime();
		if (_timeToFeed & !_alreadyFed) {
			// if its time to feed, set concentration in the reactor to feed
			setValue(getValue() * (1 - _feedFraction) + _feedConcentration
					* _feedFraction);
			_alreadyFed = true;
		} else {
			// otherwise, compute the concentration dynamics from a mass balance
			float c = getValue() + tstep
					* rate(getCurrentGlobalRateFromMassBalance());
			setValue(c > 0 ? c : 0);
		}
		_timeToFeed = false;
	}

	/**
	 * converts rate into an ammount to compute chang in bulk concentration
	 * 
	 * @param r
	 *            the rate in [MT^-1(L of computational volume)^-3]
	 * @return the global rate for this solute species [M(L of reactor)^-3T^-1]
	 */
	private float rate(float r) {
		float a = Model.model().getComputationalVolumeMultiplier()
				/ Model.model().getReactorVolume();
		return a * r;
	}

	/**
	 * Ensures that a iteration step is computed before the feed and one is
	 * computed right after the feed
	 * 
	 * @return the time step contraint
	 */
	public TimeStepConstraint getMaximumTimeStep() {
		TimeStepConstraint constraint = super.getMaximumTimeStep();
		// Determine the time until next regime change
		float t = Model.model().getTime();
		float timeToChange = 0;
		float timeSinceChange = 0;
		// The time already in the present cycle
		float cycleTime = (t % _cycle._cycleDuration);
		// In famine regime, compute time to end of famine period
		timeToChange = _cycle._cycleDuration - cycleTime;
		timeSinceChange = cycleTime;
		if ((timeToChange < _smallTimeInterval)
				| (timeSinceChange < _smallTimeInterval)) {
			// set system to feed SBR
			_timeToFeed = true;
			constraint.setTimeStep(2 * _smallTimeInterval);
			constraint.setName("Bulk concentration (" + _species.getName()
					+ ") during feed");
			return constraint;
		} else {
			// after feed, reset "already fed" switch
			_alreadyFed = false;
			// This makes sure that feeding only occurs once
		}
		// Add the predefined small value
		float cycleConstraint = timeToChange;
		float mBR = getMaximumTimeStepForMassBalance(); // mass balance
		// restriction
		mBR = (mBR < _smallTimeInterval ? _smallTimeInterval : mBR);
		// final constraint
		if (mBR < cycleConstraint) {
			// throw exception if constraint is NaN
			if (Float.isNaN(mBR))
				throw new ModelRuntimeException(
						"getMaximumTimeStep for bulk concentration with feed value "
								+ _feedConcentration + " returned a NaN value");
			//
			constraint.setTimeStep(mBR);
			constraint.setName("Bulk concentration (" + _species.getName()
					+ ") in mass balance, CS = " + getValue());
			return constraint;
		}
		// in case the liniting step is the cycle constraint
		constraint.setTimeStep(cycleConstraint);
		constraint.setName("Bulk concentration (" + _species.getName()
				+ ") in cycle constraint");
		return constraint;
	}

	/**
	 * Guaranty that bulk concentration does not change in steps greater than
	 * FRACTION of the value in previous iteration.
	 * 
	 * @return the maximum time step to fullfil covergence condition
	 */
	private float getMaximumTimeStepForMassBalance() {
		float cBulk = getValue();
		if (cBulk <= _precisionInConcentration)
			// no restrictions are applied if concentration if below the
			// precision
			// threshold
			return Float.POSITIVE_INFINITY;
		float r = rate(getCurrentGlobalRateFromDiffusionReaction());
		// compute mass balance
		// condition 1 - concentration should not change in steps greater than
		// FRACTION
		float t1 = _maxFractionToDecrease * cBulk / Math.abs(r);
		// if t1 is 0 (concentration 0) no time restrain (return infinite)
		return t1;
	}
}