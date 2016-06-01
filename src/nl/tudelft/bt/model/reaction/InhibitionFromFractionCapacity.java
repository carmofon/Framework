/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.multigrid.*;

/**
 * Implements kinetics based on the equation:
 * 
 * (fmax - f)/( (fmax - f) + k )
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class InhibitionFromFractionCapacity extends ProcessFactor {
	protected ParticulateSpecies _species1;

	protected ParticulateSpecies _species2;

	private float _k;

	private float _fmax;

	/**
	 * @param c1
	 *            species that defines saturation
	 * @param c2
	 *            species that defines total
	 * @param k
	 */
	public InhibitionFromFractionCapacity(ParticulateSpecies c1,
			ParticulateSpecies c2, float k, float fmax) {
		_species1 = c1;
		_species2 = c2;
		_k = k;
		_fmax = fmax;
	}

	public float getValue() {
		float v1 = _species1.getValue();
		float v2 = _species2.getValue();
		float f = v1 / v2;
		if ((v2 == 0) | (f > _fmax))
			return 0;
		return (_fmax - f) / (_k + (_fmax - f));
	}

	public float getMaximumValue() {
		float v1 = _species1.getMaximumValue();
		float v2 = _species2.getMaximumValue();
		float f = v1 / v2;
		if ((v2 == 0) | (f > _fmax))
			return 0;
		return (_fmax - f) / (_k + (_fmax - f));
	}

	public float getDerivative(SoluteSpecies c) {
		//always returns 0 since InhibitionFromFractionCapacity is always
		// defined using
		// particulate species as the involved species, so no solute can be
		// involved
		return 0f;
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return new String(super.toString() + ": species1 "
				+ _species1.getName() + _species1.getValue() + ", species2 "
				+ ", value" + _species2.getValue() + ", _k = " + _k
				+ ", _fmax = " + _fmax);
	}
	/* (non-Javadoc)
	 */
	public float getValue(Composition c) {
		float v1 = c.getSpeciesMass(_species1);
		float v2 = c.getSpeciesMass(_species2);
		float f = v1 / v2;
		if ((v2 == 0) | (f > _fmax))
			return 0;
		return (_fmax - f) / (_k + (_fmax - f));
	}
}