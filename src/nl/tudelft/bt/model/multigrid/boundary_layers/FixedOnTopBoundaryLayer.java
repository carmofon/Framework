package nl.tudelft.bt.model.multigrid.boundary_layers;

import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;

/**
 * Implements fixed boundary layer at top of system
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class FixedOnTopBoundaryLayer extends BoundaryLayer {

	/**
	 * @throws MultigridSystemNotSetException
	 */
	public FixedOnTopBoundaryLayer() throws MultigridSystemNotSetException {
		super();
	}

	public void setBoundaryLayer(ParticulateSpecies[] b,
			BoundaryConditions bc) {
		float[][][] bl = _mg[_order - 1];
		int n = bl.length;
		int m = bl[0].length;
		int l = bl[0][0].length;

		//  boundary layer on top
		for (int i = n - 1; i >= n - 3; i--)
			for (int j = 0; j < m; j++)
				for (int k = 0; k < l; k++)
					bl[i][j][k] = 1;
	}
}