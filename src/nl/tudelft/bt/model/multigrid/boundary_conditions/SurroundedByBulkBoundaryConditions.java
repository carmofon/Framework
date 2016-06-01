/*
 * File created originally on Jul 7, 2005
 */
package nl.tudelft.bt.model.multigrid.boundary_conditions;

import java.io.Serializable;

/**
 * Implements boundary conditions where all borders are constant value
 * 
 * @author jxavier
 */
public class SurroundedByBulkBoundaryConditions implements BoundaryConditions,
		Serializable {

	// // CYCLIC BOUNDARIES EVERYWHERE
	public void refreshBoundaryConditions(float[][][] u) {
		// DO NOTHING - this way the padding keeps the original bulk values

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions#
	 * isCarrier(int, int, int)
	 */
	public boolean isCarrier(int i, int j, int k) {
		// returns false since the carrier for planar geometry is outrside
		// the computational domain
		return false;
	}
}