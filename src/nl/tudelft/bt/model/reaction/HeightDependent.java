/*
 * Created on May 11, 2003
 */
package nl.tudelft.bt.model.reaction;

import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.multigrid.*;

/**
 * Implements monod kinetics for the Reaction Factor interface
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class HeightDependent extends ProcessFactor {
	private float _k;

	/**
	 * @param c
	 * @param k
	 */
	public HeightDependent(float k) {
		_k = k;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalValue(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getValue() {
		ContinuousCoordinate c = MultigridVariable.getCurrentLocation();
		return (c == null ? 0 : c.x / (_k + c.x));
	}

	public float getMaximumValue() {
		float hmax = Model.model().getMaximumBiofilmHeight();
		return hmax / (_k + hmax);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.reaction.ReactionFactor#getLocalDerivative(org.photobiofilms.phlip.ContinuousCoordinate)
	 */
	public float getDerivative(SoluteSpecies c) {
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
		return new String(super.toString() + ", _k = " + _k);
	}

	
	/* (non-Javadoc)
	 * @see nl.tudelft.bt.model.reaction.ProcessFactor#getValue(nl.tudelft.bt.model.BiomassSpecies.Composition)
	 */
	public float getValue(Composition c) {
		return getValue();
	}
}