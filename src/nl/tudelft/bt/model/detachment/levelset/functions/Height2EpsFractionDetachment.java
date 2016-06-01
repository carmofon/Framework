package nl.tudelft.bt.model.detachment.levelset.functions;

import nl.tudelft.bt.model.*;
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
public class Height2EpsFractionDetachment extends DetachmentSpeedFunction {
	private float _detachmentRateConstant;

	private ParticulateSpecies _epsCohesive;

	private ParticulateSpecies _epsCompromised;
	
	private float _gamma;

	/**
	 * @param rate
	 * @param gamma the order of biofilm cohesiveness on fEPS
	 */
	public Height2EpsFractionDetachment(float rate,
			ParticulateSpecies epsCohesive, ParticulateSpecies epsCompromised, float gamma) {
		_detachmentRateConstant = rate;
		_epsCohesive = epsCohesive;
		_epsCompromised = epsCompromised;
		_gamma = gamma;
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
				.getConcentrationInElement(c, _epsCohesive);
		float epsCompromisedConcentration = Model.model().biomassContainer
				.getConcentrationInElement(c, _epsCompromised);
		float fEps = epsCohesiveConcentration
				/ (epsCohesiveConcentration + epsCompromisedConcentration);
		// perform gamma correction of fEPS
		try {
			fEps = ExtraMath.gammaCorrection(fEps, _gamma);
		} catch (ModelException e) {
			throw new ModelRuntimeException(e.toString());
		}
		// in the case that concentrations of both EPS and EPS* are 0
		// compute only de detachment based on the dependencies of x
		// and local biomass density
		if (Float.isNaN(fEps))
			return _detachmentRateConstant * ExtraMath.sq(c.x) / density;
		return _detachmentRateConstant * ExtraMath.sq(c.x) / density / fEps;
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