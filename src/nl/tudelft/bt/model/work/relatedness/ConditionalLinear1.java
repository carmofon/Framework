/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.work.relatedness;

import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.ProcessFactor;

/**
 * Implements a zeroth order dependency of a species S. The difference between
 * zeroth order dependency and no dependency at all is that if concentration of
 * the species S is zero, the rate value is also 0.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ConditionalLinear1 extends ProcessFactor {
	protected MultigridVariable _substrateSpecies;
	protected MultigridVariable _goodSpecies;
	protected float maintenanceFraction;
	/**
	 * @param c
	 *            species that defines the factor
	 */
	public ConditionalLinear1(MultigridVariable G, MultigridVariable S, float m) {
		_goodSpecies = G;
		_substrateSpecies = S;
		maintenanceFraction = m;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float concS = _substrateSpecies.getValue();
		float concG = _goodSpecies.getValue();
		return (concS <= maintenanceFraction ? 0 : concG);
	}

	public float getMaximumValue() {
		float concS = _substrateSpecies.getMaximumValue();
		float concG = _goodSpecies.getMaximumValue();
		return (concS <= maintenanceFraction ? 0 : concG);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		if (c == _goodSpecies){
			float concG = _goodSpecies.getValue();
			return (concG == 0 ? 0 : 1);
		} 
		return 0f;
	}
}