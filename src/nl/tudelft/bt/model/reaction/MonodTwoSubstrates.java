/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.multigrid.*;

/**
 * Implements monod kinetics for the Reaction Factor interface
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class MonodTwoSubstrates extends ProcessFactor {
	protected SoluteSpecies _solute1;
	protected SoluteSpecies _solute2;
	private float _ks1;
	private float _ks2;

	/**
	 * @param c
	 * @param k
	 */
	public MonodTwoSubstrates(
		SoluteSpecies c1,
		float ks1,
		SoluteSpecies c2,
		float ks2) {
		_solute1 = c1;
		_ks1 = ks1;
		_solute2 = c2;
		_ks2 = ks2;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float conc1 = _solute1.getValue();
		float conc2 = _solute2.getValue();
		return conc1 / (_ks1 + conc1) * conc2 / (_ks2 + conc2);
	}

	public float getMaximumValue() {
		float conc1 = _solute1.getMaximumValue();
		float conc2 = _solute2.getMaximumValue();
		return conc1 / (_ks1 + conc1) * conc2 / (_ks2 + conc2);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		if (c == _solute1) {
			float conc1 = _solute1.getValue();
			float conc2 = _solute2.getValue();
			return _ks1
				/ (_ks1 + conc1)
				/ (_ks1 + conc1)
				* conc2
				/ (_ks2 + conc2);
		} else if (c == _solute2) {
			float conc1 = _solute1.getValue();
			float conc2 = _solute2.getValue();
			return conc1
				/ (_ks1 + conc1)
				* _ks2
				/ (_ks2 + conc2)
				/ (_ks2 + conc2);
		}
		return 0f; // derivative for other chemicals
	}
}
