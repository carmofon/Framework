/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.multigrid.*;

/**
 * Mass action kinetics with an array of species
 * 
 * @author Joao and Vanni
 */
public class LinearWithSpeciesArray extends ProcessFactor {
	protected MultigridVariable[] _substrateSpecies;

	/**
	 * @param c
	 *            species that defines the factor
	 */
	public LinearWithSpeciesArray(MultigridVariable[] S) {
		_substrateSpecies = S;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float concS = 0;
		for (int i = 0; i < _substrateSpecies.length; i++) {
			concS += _substrateSpecies[i].getValue();
		}
		return (concS < 0 ? 0 : concS);
	}

	public float getMaximumValue() {
		float concS = 0;
		for (int i = 0; i < _substrateSpecies.length; i++) {
			concS += _substrateSpecies[i].getMaximumValue();
		}
		return (concS < 0 ? 0 : concS);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		for (int i = 0; i < _substrateSpecies.length; i++) {
			if (c == _substrateSpecies[i]){
				float concS = _substrateSpecies[i].getValue();
				return (concS == 0 ? 0 : 1);
			} 			
		}
		return 0f;
	}
}