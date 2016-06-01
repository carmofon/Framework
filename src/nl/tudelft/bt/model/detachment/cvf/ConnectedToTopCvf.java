/*
 * Created on Sep 14, 2004
 */
package nl.tudelft.bt.model.detachment.cvf;

/**
 * @author jxavier
 */
public class ConnectedToTopCvf extends ConnectedVolumeFilter {

	/**
	 * @param n
	 * @param m
	 * @param l
	 */
	public ConnectedToTopCvf(int n, int m, int l) {
		super(n, m, l);
	}
	/**
	 * The initiator starts with the top values
	 */
	protected void initializeCvf() {
		for (int j = 0; j < _m; j++) {
			for (int k = 0; k < _l; k++) {
				// all points are initialized to false except the points at the
				//top which are true in the _matrixToFilter
				for (int i = 0; i < _n - 1; i++) {
					_cvf[i][j][k] = false;
				}
				// initialize cvf with points on top of system,
				// i.e. copy values of the top layer
				_cvf[_n - 1][j][k] = _matrixToFilter[_n - 1][j][k];
			}
		}

	}

	/**
	 * Create a new instance of the connected volume filtrator, which actually
	 * performs all the computation.
	 * 
	 * @param matrixToFilter
	 *            the data to filter, true values indicate valid material
	 * @return Returns the connected volume filtration matrix.
	 */
	public boolean[][][] computeCvf(boolean[][][] matrixToFilter) {
		// assign
		_matrixToFilter = matrixToFilter;
		//initiate the level set
		initializeCvf();
		//main loop
		//does not work with validateInThisScan
		//boolean validateInThisScan = false;
		_validateInThisIteration = 1;
		int i, j, k; //pre-allocate for optimization
		while (_validateInThisIteration > 0) {
			_validateInThisIteration = 0;
			for (i = _n - 2; i >= 0; i--) {
				//validateInThisScan = false;
				for (j = 0; j < _m; j++) {
					for (k = 0; k < _l; k++) {
						//if (validateElement(i, j, k))
						//	validateInThisScan = true;
						if (validateElement(i, j, k))
							_validateInThisIteration++;
					}
				}
				//if (!validateInThisScan)
					// break iteration if no element was validated in this scan
					//break;
			}
		}
		return _cvf;
	}
}