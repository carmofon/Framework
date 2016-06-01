/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import java.io.Serializable;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

/**
 * Abstract Reaction factor
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class ProcessFactor implements Serializable {
	abstract public float getValue();

	protected float preComputedValue;

	abstract public float getDerivative(SoluteSpecies c);

	abstract public float getMaximumValue();

	/**
	 * @param c
	 *            composition of particle for which rates are being computed
	 * @return the rate (getValue() is returned by default)
	 */
	public float getValue(BiomassSpecies.Composition c) {
		return getValue();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return this.getClass().getName();
	}
}