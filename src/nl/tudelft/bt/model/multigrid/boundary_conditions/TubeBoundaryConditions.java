/*
 * File created originally on Jul 7, 2005
 */
package nl.tudelft.bt.model.multigrid.boundary_conditions;

import java.io.Serializable;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.particlebased.tube.TubeBiomassParticleContainer;

/**
 * Implements boundary conditions for a tube morphology for the computation of
 * the solute concentrations (diffusion/reaction) where all around the tube is
 * zero-flux border, only z direction is cyclic (to simulate 2D)
 * 
 * @author jxavier
 */
public class TubeBoundaryConditions implements BoundaryConditions, Serializable {
	private float _radiusOfTubularReactor;

	private TubeBiomassParticleContainer _tube;

	private VelocityBoundaryConditions _velocityBoundaryConditions;

	private boolean[][][] _isCarrier;

	/**
	 * Implements the boundary conditions for the computation of velocity fields
	 * for the tube reactor, where the borders of the 2D computational volume
	 * are all constant value (0)
	 * 
	 * @author xavier
	 */
	public class VelocityBoundaryConditions extends TubeBoundaryConditions {
		/**
		 *  
		 */
		private VelocityBoundaryConditions() {
			super(_radiusOfTubularReactor);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions#refreshBoundaryConditions(float[][][])
		 */
		public void refreshBoundaryConditions(float[][][] u) {
			int l = u[0][0].length - 2;
			int m = u[0].length - 2;
			int n = u.length - 2;

			for (int i = 1; i <= n; i++) {
				for (int j = 1; j <= m; j++) {
					// cyclic borders (3rd dimension sides)
					u[i][j][0] = u[i][j][l];
					u[i][j][l + 1] = u[i][j][1];
				}
				for (int k = 1; k <= l; k++) {
					// constant value (y-sides)
					u[i][0][k] = 0;
					u[i][m + 1][k] = 0;
				}
			}
			for (int j = 1; j <= m; j++) {
				for (int k = 1; k <= l; k++) {
					// constant value (bottom and top)
					u[0][j][k] = 0;
					u[n + 1][j][k] = 0;
				}
			}
		}
	}

	/**
	 * Creates a new instance and new instance of tube velocity boundary
	 * conditions
	 */
	public TubeBoundaryConditions(float r) {
		_radiusOfTubularReactor = r;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions#refreshBoundaryConditions(float[][][])
	 */
	public void refreshBoundaryConditions(float[][][] u) {
		int l = u[0][0].length - 2;
		int m = u[0].length - 2;
		int n = u.length - 2;

		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= m; j++) {
				// cyclic borders (sides)
				u[i][j][0] = u[i][j][l];
				u[i][j][l + 1] = u[i][j][1];
			}
			for (int k = 1; k <= l; k++) {
				// zero-flux border (y-sides)
				u[i][0][k] = u[i][1][k];
				u[i][m + 1][k] = u[i][m][k];
			}
		}
		for (int j = 1; j <= m; j++) {
			for (int k = 1; k <= l; k++) {
				// zero flux borders (bottom)
				u[0][j][k] = u[1][j][k];
				u[n + 1][j][k] = u[n][j][k];
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions#isCarrier(int,
	 *      int, int)
	 */
	public boolean isCarrier(int i, int j, int k) {
		// initialize _tube if it was not initialized yet
		if (_tube == null)
			_tube = (TubeBiomassParticleContainer) Model.model().biomassContainer;
		if (_isCarrier == null) {
			//intialize the matrix for carrier
			_isCarrier = MultigridVariable
					.create3DBooleanMatrixWithFinnerResolution();
			int n = _isCarrier.length;
			int m = _isCarrier[0].length;
			int l = _isCarrier[0][0].length;
			for (int i2 = 0; i2 < n; i2++)
				for (int j2 = 0; j2 < m; j2++)
					for (int k2 = 0; k2 < l; k2++) {
						// get the coordinates of the center of the grid element
						// note: take into account the border padding for
						// indexes in
						// multigrid variables
						float x = ((float) i2 - 0.5f)
								* Model.model().getGridElementSide();
						float y = ((float) j2 - 0.5f)
								* Model.model().getGridElementSide();
						// check if they are located inside carrier from the
						// radius coordinate
						if (_tube.rConvertToPolar(x, y) > _radiusOfTubularReactor)
							_isCarrier[i2][j2][k2] = true;
						else
							_isCarrier[i2][j2][k2] = false;
					}
		}
		return _isCarrier[i][j][k];
	}


	/**
	 * The boundary conditions for the flow velocity
	 * 
	 * @return Returns the _velocityBoundaryConditions.
	 */
	public VelocityBoundaryConditions getVelocityBoundaryConditions() {
		if (_velocityBoundaryConditions == null)
			_velocityBoundaryConditions = new VelocityBoundaryConditions();
		return _velocityBoundaryConditions;
	}
}