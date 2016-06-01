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
public class DynamicBulkConcentrationImplicit extends BulkConcentration {
	private float _inputConcentration;
	private Model _m = Model.model();
	// Fraction of residence time to define maximum time step
	private static final float FRACTION = 1f;
	/**
	 * @param cin
	 *            input concentration
	 */
	public DynamicBulkConcentrationImplicit(float cin) {
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
		float tau = _m.getResidenceTime();
		// compute the rate
		float a = _m.getComputationalVolumeMultiplier() / _m.getReactorVolume();
		float r = a * getCurrentGlobalRateFromMassBalance();
		// compute equilibrium constant (after liniearization)
		float ceq = _inputConcentration + tau * r;
		// determine preliminary implicit concentration value
		float cprelim = (ceq + getValue() * tau / tstep)
				/ (1 + tau / tstep);
		// if cprelim has negative values, recompute assuming ceq = 0
		if (cprelim < 0)
			cprelim = (getValue() * tau / tstep) / (1 + tau / tstep);
		// set the value, with further correction for negative values
		setValue(cprelim);
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
	 * Maximum time step is a fraction of the reactor residence time
	 * 
	 * @return the maximum time step to fullfil covergence condition
	 */
	public TimeStepConstraint getMaximumTimeStep() {
		TimeStepConstraint t = super.getMaximumTimeStep();
		t.setTimeStep(_m.getResidenceTime() * FRACTION);
		t.setName(_species.getName() + " mass balance");
		return t;
	}
}
