package nl.tudelft.bt.model.multigrid.boundary_layers;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.MultigridSystemNotSetException;
import nl.tudelft.bt.model.multigrid.MultigridUtils;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Boundary layer that follows the shape of the biomass, set at a distance of
 * the biomass.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class SphericalDilationBoundaryLayer extends BoundaryLayer {
	protected float _dilationRadius;
	protected boolean[][][] _totalBiomass;

	/**
	 * 
	 * 
	 * @param dilationRadius
	 *            the dilation radius (distance to set the boundary layer
	 * @throws MultigridSystemNotSetException
	 */
	public SphericalDilationBoundaryLayer(float dilationRadius)
			throws MultigridSystemNotSetException {
		super();
		_dilationRadius = dilationRadius;
		// allocate space for the total biomass matrix
		_totalBiomass = new boolean[_n][_m][_l];
	}

	public void setBoundaryLayer(ParticulateSpecies[] b,
			BoundaryConditions bc) {
		float[][][] bl = _mg[_order - 1];
		// get the information concerning all the biomass
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
		MultigridUtils.setValues(bl, 1.0f);
		for (int i = 0; i < _n; i++)
			for (int j = 0; j < _m; j++)
				for (int k = 0; k < _l; k++)
					if (_totalBiomass[i][j][k] | bc.isCarrier(i, j, k)) {
						//if this is biomass,
						bl[i + 1][j + 1][k + 1] = 0;
					} else {
						//if liquid, check dilation sphere for biomass
						bl[i + 1][j + 1][k + 1] = checkDilationRadiusForBiomass(
								i, j, k);
					}
	}

	/**
	 * Check the dilation sphere for biomass
	 * 
	 * @return 1 if there is no biomass in the dilation radius, 0 otherwise
	 */
	protected final float checkDilationRadiusForBiomass(int n, int m, int l) {
		int nInterval = (int) Math.floor(_dilationRadius / _voxelSide);
		for (int i = (-nInterval); i <= nInterval; i++) {
			//only procedd if neighbor is within computational volume
			// boundaries
			if ((n + i >= 0) & (n + i < _n)) {
				float deltaN = (float) i * _voxelSide;
				float dilationRadiusM = (float) Math.sqrt(ExtraMath.sq(deltaN)
						+ ExtraMath.sq(_dilationRadius));
				int mInterval = (int) Math.floor(dilationRadiusM / _voxelSide);
				for (int j = (-mInterval); j <= mInterval; j++) {
					if (_l == 1) {
						//2D case
						if (_totalBiomass[n + i][cyclicIndex(m + j, _m)][0]) {
							return 0;
						}
					} else {
						//3D case
						float deltaM = (float) j * _voxelSide;
						float dilationRadiusL = (float) Math.sqrt(ExtraMath
								.sq(deltaN)
								+ ExtraMath.sq(deltaM)
								+ ExtraMath.sq(_dilationRadius));
						int lInterval = (int) Math.floor(dilationRadiusL
								/ _voxelSide);
						for (int k = (-lInterval); k <= lInterval; k++) {
							if ((i != 0) | (j != 0) | (k != 0)) {
								if (_totalBiomass[n + i][cyclicIndex(m + j, _m)][cyclicIndex(
										l + k, _l)]) {
									return 0;
								}
							}
						}
					}
				}
			}
		}
		return 1f;
	}

	/**
	 * inplments cyclic indexes
	 * 
	 * @param val
	 *            value of the index to be cycled
	 * @param limit
	 *            of the grid side
	 * @return
	 */
	private final int cyclicIndex(int val, int limit) {
		return (val < 0 ? limit + val : (val >= limit ? val - limit : val));
	}

	/**
	 * @param radius
	 *            The _dilationRadius to set.
	 */
	public void setThickness(float radius) {
		_dilationRadius = radius;
	}

	/**
	 * determine if a given point is out-of-bounds, i.e. if its distance to the
	 * center plus the dilation radius of the boundary layer are greater or
	 * equal to half the reference system side, thus defining a spherical
	 * out of bounds region
	 * 
	 * @param x
	 * @param y
	 * @param z
	 * @return true if c is out of bounds, false otherwise
	 */
	public boolean sphericalOutOfBounds(float x, float y, float z) {
		//get the center's coordinates:
		float centerX = Math.round(((float) _m - 1) / 2f) * _voxelSide;
		float centerY = Math.round(((float) _n - 1) / 2f) * _voxelSide;
		float centerZ = Math.round(((float) _l - 1) / 2f) * _voxelSide;
		// get point distance to center
		float dist = ExtraMath
				.pointDistance(centerX, centerY, centerZ, x, y, z);
		// if there is a maximum radius set (maximum thikness of biofilm)
		// a poarticle shoould be removed if this thickness is reached
		float maximumRadius = Model.model().getMaximumBiofilmHeight();
		if ((maximumRadius > 0) & (dist > maximumRadius)) {
			int lixo = 0;
			return true;
		} else
			return _referenceSystemSide / 2 < (_dilationRadius + dist);
	}
}