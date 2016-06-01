/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.multigrid.*;

/**
 * Implements a zeroth order dependency of a species S. The difference between
 * zeroth order dependency and no dependency at all is that if concentration of
 * the species S is zero, the rate value is also 0.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ZerothOrder extends ProcessFactor {
	protected MultigridVariable _species;

	/**
	 * @param c
	 *            species that defines the factor
	 */
	public ZerothOrder(MultigridVariable c) {
		_species = c;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float conc = _species.getValue();
		return (conc == 0 ? 0 : 1);
	}

	public float getMaximumValue() {
		float conc = _species.getMaximumValue();
		return (conc == 0 ? 0 : 1);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		return 0f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return new String(super.toString() + ": species " + _species.getName()
				+ ", value" + _species.getValue());
	}

	public float fractionForFactor(Composition c) {
		if (_species instanceof ParticulateSpecies) {
			ParticulateSpecies p = (ParticulateSpecies) _species;
			if (c.speciesPartOfComposition(p)) {
				float totalMass = p.getTotalMassInPresentVoxel();
				if (totalMass == 0)
					return 0;
				float m = c.getSpeciesMass(p);
				if (m > 0)
					return m / totalMass;
				return 0;
			}
		}
		return 1;
	}
}