package nl.tudelft.bt.model.work.relatedness;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Implements detachment that is proportional to square of height and a function
 * of the local EPS fraction
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class EPSonlyFractionDetachment extends DetachmentSpeedFunction {
	private float _detachmentRateConstant;

	private ParticulateSpecies _eps;
	private float _epsDensity;
	private float _p;

	/**
	 * @param rate the rate constant
	 * @param eps
	 * @param p sets the impact of EPS on cohesiveness
	 */
	public EPSonlyFractionDetachment(float rate,
			ParticulateSpecies eps, float p) {
		_detachmentRateConstant = rate;
		_eps = eps;
		_p = p;
		_epsDensity = eps.getDensity();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.detachment.DetachmentFunction#getValue(org.photobiofilms.model.ContinuousCoordinate,
	 *      float)
	 */
	public float getValue(ContinuousCoordinate c) {
		float density = Model.model().biomassContainer.getElementDensity(c);
		float epsCohesiveConcentration = Model.model().biomassContainer
				.getConcentrationInElement(c, _eps);
		float fEps = epsCohesiveConcentration
				/ _epsDensity;
		// in the case that concentrations of both EPS and EPS* are 0
		// compute only de detachment based on the dependencies of x
		// and local biomass density
		if (Float.isNaN(fEps))
			return _detachmentRateConstant * ExtraMath.sq(c.x) ;
		return _detachmentRateConstant * ExtraMath.sq(c.x) / (fEps*_p +1);
	}

	/**
	 * Set the detachment rate constant
	 * 
	 * @param r
	 *            The _detachmentRate to set.
	 */
	public void setDetachmentRateConstant(float r) {
		_detachmentRateConstant = r;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.detachment.DetachmentSpeedFunction#detachmentIsOff()
	 */
	public boolean detachmentIsOff() {
		return _detachmentRateConstant == 0;
	}
}