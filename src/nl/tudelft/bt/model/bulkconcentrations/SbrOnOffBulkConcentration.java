/*
 * Created on 26-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.bulkconcentrations;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;

/**
 * Implements bulk concentration of a solute that is turned off at the begining
 * of an SBR cycle and off somewhere along the cycle
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SbrOnOffBulkConcentration extends BulkConcentration {
	private SbrBulkConcentration.SbrCycle _cycle;

	private float _onConcentration;

	private float _timeToTurnOn;

	private boolean _turnOn;

	private float _smallTimeInterval;

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
	public SbrOnOffBulkConcentration(float onConcentration, float timeToTurnOn,
			SbrBulkConcentration.SbrCycle cycle) {
		super();
		// check that time to turn on is less than cycle time
		if (timeToTurnOn >= cycle.getCycleDuration())
			throw new ModelRuntimeException(
					"tying to initialize on/off concentration"
							+ " with time that is too long");
		// assign
		_onConcentration = onConcentration;
		_timeToTurnOn = timeToTurnOn;
		_cycle = cycle;
		// initialize
		_turnOn = false;
		setValue(0);
		// small time interval for restriction to time step
		// set to 1/1000 of the minimum among famine and feast times
		_smallTimeInterval = 0.5e-3f * cycle.getCycleDuration();
	}

	/**
	 * Computes and returns the value for the bulk concetration
	 * 
	 * @param tstep
	 *            time step
	 */
	public void computeBulkConcentration(float tstep) {
		if (!_turnOn) {
			// if its time to feed, set concentration in the reactor to feed
			setValue(0);
		} else if (_turnOn) {
			setValue(_onConcentration);
		}
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
		// The time already in the present cycle
		float cycleTime = (t % _cycle.getCycleDuration());
		// In famine regime, compute time to end of famine period
		float timeToTurnOff = _cycle.getCycleDuration() - cycleTime;
		float timeToTurnOn = _timeToTurnOn - cycleTime;
		if ((timeToTurnOff < _smallTimeInterval)
				| (cycleTime < _smallTimeInterval)) {
			// set system to feed SBR
			_turnOn = false;
			constraint.setTimeStep(2 * _smallTimeInterval);
			constraint.setName("Bulk concentration (" + _species.getName()
					+ ") during process of turning off");
			return constraint;
		} else if (((timeToTurnOn >= 0) & (timeToTurnOn < _smallTimeInterval))
				| ((timeToTurnOn < 0) & (-timeToTurnOn < _smallTimeInterval))) {
			// set system to feed SBR
			_turnOn = true;
			constraint.setTimeStep(2 * _smallTimeInterval);
			constraint.setName("Bulk concentration (" + _species.getName()
					+ ") during process of turning on");
			return constraint;
		}
		if ((timeToTurnOn > 0) & (timeToTurnOn < timeToTurnOff)) {
			constraint.setTimeStep(timeToTurnOn - _smallTimeInterval * 0.5f);
			constraint.setName("Bulk concentration (" + _species.getName()
					+ ") about to be turned on");
			return constraint;
		}
		// in case the liniting step is the cycle constraint
		constraint.setTimeStep(timeToTurnOff - _smallTimeInterval * 0.5f);
		constraint.setName("Bulk concentration (" + _species.getName()
				+ ") about to be turned off");
		return constraint;
	}
}