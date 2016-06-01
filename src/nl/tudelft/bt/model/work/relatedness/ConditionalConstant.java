/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.work.relatedness;

import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;
import nl.tudelft.bt.model.reaction.ProcessFactor;

/**
 * Implements constant kinetics for reaction factor
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ConditionalConstant extends ProcessFactor {
	protected MultigridVariable _substrateSpecies;
	protected float _maintenanceFraction;

	/**
	 * @param k
	 *            constant value
	 */
	public ConditionalConstant(MultigridVariable s, float m) {
		_substrateSpecies = s;
		_maintenanceFraction = m;
	}

	public float getValue() {
		float _concS = _substrateSpecies.getValue();
		return (_concS <= _maintenanceFraction ? 0 : 1);
	}

	public float getMaximumValue() {
		float _concS = _substrateSpecies.getMaximumValue();
		return (_concS <= _maintenanceFraction ? 0 : 1);
	}

	public float getDerivative(SoluteSpecies c) {
		return 0f;
	}

}
