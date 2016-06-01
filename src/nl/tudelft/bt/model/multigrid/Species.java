/*
 * File created originally on Jan 18, 2005
 */
package nl.tudelft.bt.model.multigrid;

import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;

/**
 * Common base class for SoluteSpecies and ParticulateSpecies
 * 
 * @author jxavier
 */
public abstract class Species extends MultigridVariable {
	/**
	 * @throws MultigridSystemNotSetException
	 */
	public Species() throws MultigridSystemNotSetException {
		super();
	}

	/**
	 * @return the total mass of solute in the present grid node
	 */
	public float getTotalMassInPresentVoxel() {
		return getValue() * _voxelVolume;
	}

}
