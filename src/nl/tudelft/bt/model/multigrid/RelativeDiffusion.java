package nl.tudelft.bt.model.multigrid;

import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.BoundaryLayer;

/**
 * Implements multigrid variable for the relative diffusion
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class RelativeDiffusion extends MultigridVariable { 
	/**
	 * @throws MultigridSystemNotSetException
	 */
	public RelativeDiffusion() throws MultigridSystemNotSetException {
		super();
		_name = "relativeDiffusivity";
	}

	/**
	 * For now, diffusivity is the same throughout the space. Lower diffusivity
	 * values may be set by changing this functions
	 * 
	 * @param bac
	 * @param bl
	 */
	protected void computeValues(ParticulateSpecies[] bac, BoundaryLayer bl,
			BoundaryConditions bc) {
		//set the value of the finner grid
		float[][][] d = _mg[_order - 1];
		MultigridUtils.setValues(d, 1.0f);
		for (int i = 0; i < d.length; i++)
			for (int j = 0; j < d[i].length; j++)
				for (int k = 0; k < d[i][j].length; k++)
					// if inside carrier, relative diffusivity should be
					// very close to 0
					if (bc.isCarrier(i, j, k)) {
						d[i][j][k] = 1.0e-20f;
					}
		//update corser grid copies
		updateMultigridCopies();
	}
}