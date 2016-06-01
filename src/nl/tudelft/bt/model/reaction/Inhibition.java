/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.multigrid.*;

/**
 * Implements non-competitive inhibition kinetics for reaction factor
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Inhibition extends ProcessFactor {
	protected MultigridVariable _species;
	private float _k;

	/**
	 * @param c
	 * @param k
	 */
	public Inhibition(MultigridVariable c, float k) {
		_species = c;
		_k = k;
	}

	public float getValue() {
		float conc = _species.getValue();
		conc = (conc < 0 ? 0 : conc);
		return _k / (_k + conc);
	}

	public float getMaximumValue() {
		float conc = _species.getMaximumValue();
		conc = (conc < 0 ? 0 : conc);
		return _k / (_k + conc);
	}

	public float getDerivative(SoluteSpecies c) {
		if (c == _species) {
			float conc = _species.getValue();
			conc = (conc < 0 ? 0 : conc);
			return -_k / ((_k + conc) * (_k + conc));
		}
		return 0f;
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return new String(super.toString() + ": species1 " + ", value"
				+ _species.getValue() + ", _k = " + _k);
	}
}