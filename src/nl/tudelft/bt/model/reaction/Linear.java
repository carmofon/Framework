/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.multigrid.*;

/**
 * Implements first order kinetics kinetics for the Reaction Factor interface
 * 
 * @author Joao Xavier (j.xavier@mskcc.org)
 */
public class Linear extends ProcessFactor {
	protected MultigridVariable _species;


	/**
	 * @param c
	 * @param k
	 */
	public Linear(MultigridVariable c) {
		_species = c;

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float conc = _species.getValue();
		conc = (conc < 0 ? 0 : conc);
		return conc;
	}

	public float getMaximumValue() {
		float conc = _species.getMaximumValue();
		conc = (conc < 0 ? 0 : conc);
		return conc;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		if (c == _species) {
			float conc = _species.getValue();
			return (conc < 0 ? 0 : 1);
		}
		return 0f; // derivative for other chemicals
	}


	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return new String(super.toString() + ": species " + _species.getName()
				+ ", value" + _species.getValue() );
	}

	
	/* (non-Javadoc)
	 * @see nl.tudelft.bt.model.reaction.ProcessFactor#getValue(nl.tudelft.bt.model.BiomassSpecies.Composition)
	 */
	public float getValue(Composition c) {
		// if the _species is a particulate species, return a rate that
		// reflects the fraction of the species in the grid element (present
		// voxel
		if (_species instanceof ParticulateSpecies) {
			ParticulateSpecies p = (ParticulateSpecies) _species;
			if (c.speciesPartOfComposition(p)) {
				float totalMass = p.getTotalMassInPresentVoxel();
				float m = c.getSpeciesMass(p);
				// case where totalMass is zero or too small for 
				// precision purposes
				if (Float.isInfinite(m / totalMass))
					return 0;
				if (m > 0) {
					return m / totalMass * getValue();
				}
				return 0;
			}
		}
		return getValue();
	}
}