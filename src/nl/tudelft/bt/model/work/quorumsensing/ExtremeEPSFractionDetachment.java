package nl.tudelft.bt.model.work.quorumsensing;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Infinite detachment when local feps = 0, and no detachment when feps = 1
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ExtremeEPSFractionDetachment extends DetachmentSpeedFunction {
	private float _detachmentRateConstant;

	private ParticulateSpecies _eps;
	private float _epsDensity;

	/**
	 * @param rate the rate constant
	 * @param eps
	 * @param p sets the impact of EPS on cohesiveness
	 */
	public ExtremeEPSFractionDetachment(float rate,
			ParticulateSpecies eps, float p) {
		_detachmentRateConstant = rate;
		_eps = eps;
		_epsDensity = eps.getDensity();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.detachment.DetachmentFunction#getValue(org.photobiofilms.model.ContinuousCoordinate,
	 *      float)
	 */
	public float getValue(ContinuousCoordinate c) {
		float epsCohesiveConcentration = Model.model().biomassContainer
				.getConcentrationInElement(c, _eps);
		float fEps = epsCohesiveConcentration
				/ _epsDensity;
		// in the case that concentrations of both EPS and EPS* are 0
		// compute only de detachment based on the dependencies of x
		// and local biomass density
		if (Float.isNaN(fEps))
			return _detachmentRateConstant * ExtraMath.sq(c.x) ;
		if (fEps > 1) 
			fEps = 1f;
		//_detachmentRateConstant = 0;
		return 0.00001f + _detachmentRateConstant * ExtraMath.sq(c.x)*(1f-fEps) / (fEps);
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