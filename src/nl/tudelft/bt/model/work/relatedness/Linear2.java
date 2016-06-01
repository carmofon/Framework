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
public class Linear2 extends ProcessFactor {
	protected MultigridVariable _substrateSpecies;
	protected MultigridVariable _goodSpecies;

	/**
	 * @param c
	 *            species that defines the factor
	 */
	public Linear2(MultigridVariable S, MultigridVariable G) {
		_substrateSpecies = S;
		_goodSpecies = G;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float concS = _substrateSpecies.getValue();
		float concG = _goodSpecies.getValue();
		return (concS*concG == 0 ? 0 : concS*concG);
	}

	public float getMaximumValue() {
		float concS = _substrateSpecies.getMaximumValue();
		float concG = _goodSpecies.getMaximumValue();
		return (concS*concG == 0 ? 0 : concS*concG);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		if (c == _substrateSpecies){
			float concS = _substrateSpecies.getValue();
			float concG = _goodSpecies.getValue();
			return concG;
		} else if (c == _goodSpecies){
			float concS = _substrateSpecies.getValue();
			float concG = _goodSpecies.getValue();
			return concS;
		}
		return 0f;
	}
}