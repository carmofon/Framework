package nl.tudelft.bt.model.multigrid.boundary_layers;

import nl.tudelft.bt.model.detachment.cvf.ConnectedToTopCvf;
import nl.tudelft.bt.model.detachment.cvf.ConnectedVolumeFilter;
import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;

/**
 * Boundary layer that follows the shape of the biomass, set at a distance of
 * the biomass.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DilationFillHolesBoundaryLayer extends SphericalDilationBoundaryLayer {
	private ConnectedVolumeFilter cvf;
	private boolean[][][] _boundaryLayer;

	/**
	 * @param dilationRadius
	 *            the dilation radius (distance to set the boundary layer
	 * @throws MultigridSystemNotSetException
	 */
	public DilationFillHolesBoundaryLayer(float dilationRadius)
			throws MultigridSystemNotSetException {
		super(dilationRadius);
		_boundaryLayer = new boolean[_n][_m][_l];
		cvf = new ConnectedToTopCvf(_n, _m, _l);
	}

	/**
	 * Here we define a dilation boundary condition but then fill the holes
	 * inside the cell cluster (the tumor or the biofilm) using a connected
	 * volume filtration
	 */
	public void setBoundaryLayer(ParticulateSpecies[] b, BoundaryConditions bc) {
		float[][][] bl = _mg[_order - 1];
		// get the information on all the biomass
		for (int i = 0; i < _n; i++)
			for (int j = 0; j < _m; j++)
				for (int k = 0; k < _l; k++) {
					_totalBiomass[i][j][k] = false;
					for (int sp = 0; sp < b.length; sp++)
						if (b[sp]._mg[_order - 1][i + 1][j + 1][k + 1] > 0) {
							_totalBiomass[i][j][k] = true;
							break;
						}
				}
		for (int i = 0; i < _n; i++)
			for (int j = 0; j < _m; j++)
				for (int k = 0; k < _l; k++)
					if (_totalBiomass[i][j][k] | bc.isCarrier(i, j, k)) {
						// if this is biomass,
						// bl[i + 1][j + 1][k + 1] = 0;
						_boundaryLayer[i][j][k] = false;
					} else {
						// if liquid, check dilation sphere for biomass
						//bl[i + 1][j + 1][k + 1] = checkDilationRadiusForBiomass(
						//		i, j, k);
						if (checkDilationRadiusForBiomass(i, j, k) == 1)
							_boundaryLayer[i][j][k] = true;
						else
							_boundaryLayer[i][j][k] = false;
					}
		// Compute the connected volume filtration
		boolean[][][] cvf2 = cvf.computeCvf(_boundaryLayer);
		// Update the boundary layer information (matrix bl)
		// to include the connected volume filtration
		for (int i = 0; i < _n; i++)
			for (int j = 0; j < _m; j++)
				for (int k = 0; k < _l; k++)
					if (cvf2[i][j][k])
						bl[i + 1][j + 1][k + 1] = 1;
					else
						bl[i + 1][j + 1][k + 1] = 0;
	}
}