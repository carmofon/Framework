/*
 * Created on Sep 14, 2004
 */
package nl.tudelft.bt.model.detachment.cvf;

/**
 * @author jxavier
 */
public class ConnectedToSideBordersCvf extends ConnectedToTopCvf {
	/**
	 * @param n
	 * @param m
	 * @param l
	 */
	public ConnectedToSideBordersCvf(int n, int m, int l) {
		super(n, m, l);
	}

	/**
	 * The initiator starts with the top values
	 */
	protected void initializeCvf() {
		//initialize center to false
		for (int i = 1; i < _n-1; i++) {
			for (int j = 1; j < _m-1; j++) {
				for (int k = 1; k < _l-1; k++) {
					_cvf[i][j][k] = false; 
				}
			}
		}
		//bottom and top Borders
		for (int j = 0; j < _m; j++) {
			for (int k = 0; k < _l; k++) {
				// top layer
				_cvf[0][j][k] = _matrixToFilter[0][j][k];
				// bottom layer
				_cvf[_n - 1][j][k] = _matrixToFilter[_n - 1][j][k];
			}
		}
		// side borders
		for (int i = 0; i < _n; i++) {
			for (int k = 0; k < _l; k++) {
				// left layer
				_cvf[i][0][k] = _matrixToFilter[i][0][k];
				// right layer
				_cvf[i][_m - 1][k] = _matrixToFilter[i][_m - 1][k];
			}
		}
		// front and back borders
		for (int i = 0; i < _n; i++) {
			for (int j = 0; j < _m; j++) {
				// front layer
				_cvf[i][j][0] = _matrixToFilter[i][j][0];
				if (_l > 1) {
					// back layer (only compute if 3D
					_cvf[i][j][_l - 1] = _matrixToFilter[i][j][_l - 1];
				}
			}
		}
	}
}