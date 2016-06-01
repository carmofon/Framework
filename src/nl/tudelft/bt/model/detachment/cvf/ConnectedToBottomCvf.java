/*
 * Created on Sep 14, 2004
 */
package nl.tudelft.bt.model.detachment.cvf;

/**
 * Implements a connceted volume filtrator to remove everything not connected to
 * the bottom
 * 
 * @author jxavier
 */
public class ConnectedToBottomCvf extends ConnectedVolumeFilter {
	/**
	 * @param n
	 * @param m
	 * @param l
	 */
	public ConnectedToBottomCvf(int n, int m, int l) {
		super(n, m, l);
	}

	/**
	 * The initiator starts with the top values
	 */
	protected void initializeCvf() {
		for (int j = 0; j < _m; j++) {
			for (int k = 0; k < _l; k++) {
				//initialize values to false except the ones which are true in
				//bottom of _matrixToFilter
				for (int i = 1; i < _n; i++) {
					_cvf[i][j][k] = false;
				}
				// initialize cvf with points on top of system,
				// i.e. copy values of the top layer
				_cvf[0][j][k] = _matrixToFilter[0][j][k];
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
		//NOTE: use of validate this scan does not work with the
		//cyclic borders!!
		//boolean validateInThisScan = false;
		_validateInThisIteration = 1;
		int i, j, k; //pre-allocate for optimization
		while (_validateInThisIteration > 0) {
			_validateInThisIteration = 0;
			for (i = 1; i < _n; i++) {
				//validateInThisScan = false;
				for (j = 0; j < _m; j++) {
					for (k = 0; k < _l; k++) {
						if (validateElement(i, j, k))
							_validateInThisIteration++;
						//if (validateElement(i, j, k))
						//	validateInThisScan = true;
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