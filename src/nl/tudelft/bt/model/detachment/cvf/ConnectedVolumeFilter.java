/*
 * Created on 29-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.detachment.cvf;

import java.io.Serializable;

/**
 * Implements connected volume filtration operation. Provides a base class for
 * geometry specific connected volume filtrators.
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class ConnectedVolumeFilter implements Serializable {
	protected int _n;

	protected int _m;

	protected int _l;

	protected boolean[][][] _cvf;

	protected boolean[][][] _matrixToFilter;

	protected int _validateInThisIteration;

	/**
	 * Initialize the _cvf matrix
	 */
	public ConnectedVolumeFilter(int n, int m, int l) {
		_n = n;
		_m = m;
		_l = l;
		_cvf = new boolean[_n][_m][_l];
	}

	/**
	 * Implements actual computation of the cvf algorithm. Each class deriving
	 * from ConnectedVolumeFilter implements this method in a diverse way which
	 * most suits the way that the filtrator performs its computation.
	 * 
	 * @param matrixToFilter
	 *            the data to filter, true values indicate valid material
	 * @return Returns the connected volume filtration matrix.
	 */
	public abstract boolean[][][] computeCvf(boolean[][][] matrixToFilter);

	/**
	 * Initiate the cvf process.
	 */
	protected abstract void initializeCvf();

	/**
	 * @param i
	 * @param j
	 * @param k
	 * @return true is element is validated
	 */
	public boolean validateElement(int i, int j, int k) {
		if (!_cvf[i][j][k] & _matrixToFilter[i][j][k]) {
			// check if any of the neighbors is valid
			// bottom neighbor
			if (((i > 0) ? _cvf[i - 1][j][k] : false)
			// y-side neighbors:
					| (_cvf[i][j < _m - 1 ? j + 1 : 0][k])
					| (_cvf[i][j > 0 ? j - 1 : _m - 1][k]) // z-side
					// neighbors:
					| ((_l > 1) & (_cvf[i][j][k < _l - 1 ? k + 1 : 0]))
					| ((_l > 1) & (_cvf[i][j][k > 0 ? k - 1 : _l - 1]))
					// top neighbor:
					| ((i < _n - 1) ? (_cvf[i + 1][j][k]) : false)) {
				_cvf[i][j][k] = true;
				return true;
			}
		}
		return false;
	}
}