/*
 * File created originally on Jul 7, 2005
 */
package nl.tudelft.bt.model.multigrid.boundary_conditions;

import java.io.Serializable;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Implements boundary conditions where all borders have a limited flux
 * 
 * @author jxavier
 */
public class TissueBoundaryConditions implements BoundaryConditions,
		Serializable {

	private float _externalTransfer;

	/**
	 * @param externalTransfer
	 *            scale for external mass transfer
	 */
	public TissueBoundaryConditions(float externalTransfer) {
		super();
		if ((externalTransfer > 1) | (externalTransfer < 0))
			throw new ModelRuntimeException(
					"externalTransfer must be > 0 and < 1");
		_externalTransfer = externalTransfer;
	}

	public void refreshBoundaryConditions(float[][][] u) {
		int l = u[0][0].length - 2;
		int m = u[0].length - 2;
		int n = u.length - 2;

		float bulkConcentration = u[0][0][0];
		float aux1 = bulkConcentration * _externalTransfer;
		float aux2 = 1 - _externalTransfer;

		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= m; j++) {
				u[i][j][0] = aux1 + aux2 * u[i][j][1];
				u[i][j][l + 1] = aux1 + aux2 * u[i][j][l];
			}
			for (int k = 1; k <= l; k++) {
				u[i][0][k] = aux1 + aux2 * u[i][1][k];
				u[i][m + 1][k] = aux1 + aux2 * u[i][m][k];
			}
		}
		for (int j = 1; j <= m; j++) {
			for (int k = 1; k <= l; k++) {
				u[0][j][k] = aux1 + aux2 * u[1][j][k];
				u[n + 1][j][k] = aux1 + aux2 * u[n][j][k];
			}
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions#
	 * isCarrier(int, int, int)
	 */
	public boolean isCarrier(int i, int j, int k) {
		// returns false since the carrier for planar geometry is outside
		// the computational domain
		return false;
	}
}