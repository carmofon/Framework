/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.multigrid.*;

/**
 * Implements monod kinetics for the Reaction Factor interface. This class adds
 * maintenance, with rate defined as a fraction f of the maximum specific growth
 * rate. Rather than being negative, the rate is 0 when S < f/(1-f)*k
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu)
 */
public class SaturationWithMaintenance extends ProcessFactor {
	protected MultigridVariable _species;

	private float _k;

	private float _f;

	/**
	 * @param c
	 * @param k
	 * @param f
	 *            the fraction that goes to maintenance
	 */
	public SaturationWithMaintenance(MultigridVariable c, float k, float f) {
		if ((f >= 1) || (f <= 0))
			throw new RuntimeException("Maintenance fraction f = " + f
					+ " must be < 1 and > 0");
		_species = c;
		_k = k;
		_f = f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		float conc = _species.getValue();
		if (conc < (_f / (1 - _f) * _k))
			return 0;
		return conc / (_k + conc) - _f;
	}

	public float getMaximumValue() {
		float conc = _species.getMaximumValue();
		if (conc < (_f / (1 - _f) * _k))
			return 0;
		return conc / (_k + conc) - _f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
		if (c == _species) {
			float conc = _species.getValue();
			if (conc < (_f / (1 - _f) * _k))
				return 0;
			return _k / ((_k + conc) * (_k + conc));
		}
		return 0f; // derivative for other chemicals
	}

	/**
	 * @param _k
	 *            The _k to set.
	 */
	public void setK(float k) {
		this._k = k;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return new String(super.toString() + ": species " + _species.getName()
				+ ", value" + _species.getValue() + ", _k = " + _k);
	}

	/*
	 * (non-Javadoc)
	 * 
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