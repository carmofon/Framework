/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;
import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.multigrid.*;
/**
 * Implements inhibition kinetics based on a concentration fraction for biomass
 * components (particulate species) reaction factor
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class InhibitionFromFraction extends ProcessFactor {
	protected ParticulateSpecies _species1;
	protected ParticulateSpecies _species2;
	private float _k;
	/**
	 * @param c1
	 *            species that defines inhibition
	 * @param c2
	 *            species that defines total
	 * @param k
	 */
	public InhibitionFromFraction(ParticulateSpecies c1, ParticulateSpecies c2,
			float k) {
		_species1 = c1;
		_species2 = c2;
		_k = k;
	}
	public float getValue() {
		float v2 = _species2.getValue();
		return (v2 > 0) ? (_k / (_k + _species1.getValue() / v2)) : 0;
	}
	public float getMaximumValue() {
		float v2 = _species2.getMaximumValue();
		return (v2 > 0) ? (_k / (_k + _species1.getMaximumValue() / v2)) : 0;
	}
	public float getDerivative(SoluteSpecies c) {
		//always returns 0 since InhibitionFromFraction is always defined using
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
				+ _species1.getName() + ", value" + _species1.getValue()
				+ ", species2 " + ", value" + _species2.getValue()
				+ _species2.getName() + ", _k = " + _k);
	}
	/* (non-Javadoc)
	 */
	public float getValue(Composition c) {
		float v1 = c.getSpeciesMass(_species1);
		float v2 = c.getSpeciesMass(_species2);
		return (v2 > 0) ? (_k / (_k + v1 / v2)) : 0;
	}
}