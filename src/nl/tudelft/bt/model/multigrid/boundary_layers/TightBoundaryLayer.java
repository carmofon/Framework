package nl.tudelft.bt.model.multigrid.boundary_layers;

import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.MultigridUtils;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;

/**
 * Constructs a boundary layer which is tight to the biofilm. It shape follows
 * the surface of the biofilm.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class TightBoundaryLayer extends BoundaryLayer {

	/**
	 * @throws MultigridSystemNotSetException
	 */
	public TightBoundaryLayer() throws MultigridSystemNotSetException {
		super();
	}

	public void setBoundaryLayer(ParticulateSpecies[] b,
			BoundaryConditions bc) {
		float[][][] bl = _mg[_order - 1];
		MultigridUtils.setValues(bl, 1.0f);
		for (int i = 0; i < bl.length; i++)
			for (int j = 0; j < bl[i].length; j++)
				for (int k = 0; k < bl[i][j].length; k++)
					for (int sp = 0; sp < b.length; sp++)
						if (b[sp]._mg[_order - 1][i][j][k] > 0) {
							bl[i][j][k] = 0;
							break;
						}
	}
}