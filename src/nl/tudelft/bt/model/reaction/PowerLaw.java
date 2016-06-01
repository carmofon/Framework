/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.multigrid.*;

/**
 * Implements power law kinetics for reaction factor
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class PowerLaw extends ProcessFactor {
	protected SoluteSpecies _chemical;
	private float _k;

	/**
	 * @param c
	 * @param k
	 */
	public PowerLaw(SoluteSpecies c, float k) {
		_chemical = c;
		_k = k;
	}

	/* (non-Javadoc)
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float conc = _chemical.getValue();
		return (float)Math.pow(conc, _k);
	}

	public float getMaximumValue() {
		float conc = _chemical.getMaximumValue();
		return (float)Math.pow(conc, _k);
	}

	/* (non-Javadoc)
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		if (c == _chemical) {
			float conc = _chemical.getValue();
			return _k * (float)Math.pow( conc , _k-1 );
		}
		return 0f;
	}
}
