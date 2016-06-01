/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

/**
 * Implements constant kinetics for reaction factor
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Constant extends ProcessFactor {
	private float _k;

	/**
	 * @param k
	 *            constant value
	 */
	public Constant(float k) {
		_k = k;
	}

	public float getValue() {
		return _k;
	}

	public float getMaximumValue() {
		return _k;
	}

	public float getDerivative(SoluteSpecies c) {
		return 0f;
	}

}
