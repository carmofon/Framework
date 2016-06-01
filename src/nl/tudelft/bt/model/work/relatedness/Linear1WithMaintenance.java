/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.work.relatedness;

import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.ProcessFactor;

/**
 * Implements a first order dependency of a species S with maintenance
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Linear1WithMaintenance extends ProcessFactor {
	protected MultigridVariable _substrateSpecies;
	protected float maintenanceFraction;

	/**
	 * @param c
	 *            species that defines the factor
	 */
	public Linear1WithMaintenance(MultigridVariable S, float m) {
		_substrateSpecies = S;
		maintenanceFraction = m;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float v = _substrateSpecies.getValue() - maintenanceFraction;
		return (v <= 0 ? 0 : v);
	}

	public float getMaximumValue() {
		float v = _substrateSpecies.getMaximumValue() - maintenanceFraction;
		return (v <= 0 ? 0 : v);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		if (c == _substrateSpecies){
			float concS = _substrateSpecies.getValue();
			return (concS == 0 ? 0 : 1);
		} 
		return 0f;
	}
}