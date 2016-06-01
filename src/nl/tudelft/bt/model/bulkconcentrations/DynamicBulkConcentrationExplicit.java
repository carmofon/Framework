/*
 * Created on 26-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.bulkconcentrations;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;
/**
 * Implements a bulk concentration manager class for multidimensional biofilm
 * modelling. In this case, the solute enters the reactor
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DynamicBulkConcentrationExplicit extends BulkConcentration {
	private float _inputConcentration;
	private Model _m = Model.model();
	// Maximum fraction to change at each iteration
	private static final float FRACTION = 0.05f;
	/**
	 * @param cin
	 *            input concentration
	 */
	public DynamicBulkConcentrationExplicit(float cin) {
		super();
		_inputConcentration = cin;
		setValue(cin);
	}
	/**
	 * Computes and returns the value for the bulk concetration
	 * 
	 * @param tstep
	 *            time step
	 */
	public void computeBulkConcentration(float tstep) {
		float c = getValue() + tstep * rate();
		setValue(c > 0 ? c : 0);
	}
	/**
	 * @return the global rate for this solute species [ML^-3T^-1]
	 */
	private float rate() {
		float r = getCurrentGlobalRateFromDiffusionReaction();
		float a = _m.getComputationalVolumeMultiplier() / _m.getReactorVolume();
		float t = _m.getResidenceTime();
		return (a * r + (_inputConcentration - getValue()) / t);
	}
	/**
	 * Guaranty that bulk concentration does not change in steps greater than
	 * 10% of the value in previous iteration.
	 * 
	 * @return the maximum time step to fullfil covergence condition
	 */
	public TimeStepConstraint getMaximumTimeStep() {
		// the super method computes the global conversion rates
		TimeStepConstraint constraint = super.getMaximumTimeStep();
		float r = rate();
		float cBulk = getValue();
		// condition 1 - concentration should not change in steps greater than
		//100%
		float t1 = FRACTION * cBulk / Math.abs(r);
		// if t1 is 0 (concentration 0) no time restrain
		t1 = (t1 == 0 ? Float.POSITIVE_INFINITY : t1);
		// condition 2 - if greter than input concentration and decreasing,
		// should not decrease more than the difference between bulk and
		// input
		if (((cBulk > _inputConcentration) & (r < 0))
				| ((cBulk < _inputConcentration) & (r > 0))) {
			float t2 = Math.abs(cBulk - _inputConcentration) / Math.abs(r);
			constraint.setTimeStep(Math.min(t1, t2));
			constraint.setName("Mass balance of " + _species.getName());		
			return constraint;
		}
		constraint.setTimeStep(t1);
		constraint.setName("Mass balance of " + _species.getName());		
		return constraint;
	}
}
