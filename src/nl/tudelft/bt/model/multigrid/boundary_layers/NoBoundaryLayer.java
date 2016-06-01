package nl.tudelft.bt.model.multigrid.boundary_layers;

import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.MultigridUtils;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;

/**
 * Everything is inside the boundary layer. Solves the system everywhere.
 * Substrate must come from influx
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class NoBoundaryLayer extends BoundaryLayer {

	/**
	 * @throws MultigridSystemNotSetException
	 */
	public NoBoundaryLayer() throws MultigridSystemNotSetException {
		super();
	}

	public void setBoundaryLayer(ParticulateSpecies[] b, BoundaryConditions bc) {
		// boundary layer initialized with 0 everywhere - means that calculation
		// will be everywhere
		float[][][] bl = _mg[_order - 1];
		MultigridUtils.setValues(bl, 0f);
	}
}