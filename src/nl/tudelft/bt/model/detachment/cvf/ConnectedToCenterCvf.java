/*
 * Created on Sep 14, 2004
 */
package nl.tudelft.bt.model.detachment.cvf;

import nl.tudelft.bt.model.DiscreteCoordinate;
import nl.tudelft.bt.model.exceptions.GrowthProblemException;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * @author jxavier
 */
public class ConnectedToCenterCvf extends ConnectedVolumeFilter {
	private DiscreteCoordinate nearestToCenter;

	/**
	 * @param n
	 * @param m
	 * @param l
	 */
	public ConnectedToCenterCvf(int n, int m, int l) {
		super(n, m, l);
	}
	
	/**
	 * The initiator starts with the top values
	 */
	protected void initializeCvf() {
		//start by reseting all cvf values
		for (int i = 0; i < _cvf.length; i++) {
			for (int j = 0; j < _cvf[i].length; j++) {
				for (int k = 0; k < _cvf[i][j].length; k++) {
					_cvf[i][j][k] = false;
				}
			}
		}
		//get the center's coordinates:
		int centerI = (int) Math.round(((float) _m - 1) / 2f);
		int centerJ = (int) Math.round(((float) _n - 1) / 2f);
		int centerK = (int) Math.round(((float) _l - 1) / 2f);
		// if ther center is valid use it as initiator,
		if (_matrixToFilter[centerI][centerJ][centerK]) {
			_cvf[centerI][centerJ][centerK] = true;
			nearestToCenter = new DiscreteCoordinate(centerI, centerJ, centerK);
			return;
		}
		// otherwise find nearest valid point
		float minDistance = Float.POSITIVE_INFINITY;
		nearestToCenter = new DiscreteCoordinate();
		for (int i = 0; i < _n; i++)
			for (int j = 0; j < _m; j++)
				for (int k = 0; k < _l; k++) {
					if (_matrixToFilter[i][j][k]) {
						//compute distance
						float dist = (float) Math.sqrt(ExtraMath
								.sq(centerI - i)
								+ ExtraMath.sq(centerJ - j)
								+ ExtraMath.sq(centerK - k));
						if (dist < minDistance) {
							// update point
							nearestToCenter.i = i;
							nearestToCenter.j = j;
							nearestToCenter.k = k;
							// update distance
							minDistance = dist;
						}
					}
				}
		if (minDistance == Float.POSITIVE_INFINITY) {
			throw new RuntimeException(
					"No biofilm found in connected volume filtration");
		}
		_cvf[nearestToCenter.i][nearestToCenter.j][nearestToCenter.k] = true;
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
		_validateInThisIteration = 1;
		int i, j, k; //pre-allocate for optimization
		while (_validateInThisIteration > 0) {
			_validateInThisIteration = 0;
			//1
			for (i = nearestToCenter.i; i < _n; i++) {
				for (j = nearestToCenter.j; j < _m; j++) {
					for (k = nearestToCenter.k; k < _l; k++) {
						//TODO delete after testing phase
						try {
							validateElement(i, j, k);
						} catch (ArrayIndexOutOfBoundsException e) {
							int lixo = 0;
						}
						//end
						if (validateElement(i, j, k))
							_validateInThisIteration++;
					}
				}
			}
			//2
			for (i = nearestToCenter.i; i >= 0; i--) {
				for (j = nearestToCenter.j; j < _m; j++) {
					for (k = nearestToCenter.k; k < _l; k++) {
						if (validateElement(i, j, k))
							_validateInThisIteration++;
					}
				}
			}
			//3
			for (i = nearestToCenter.i; i < _n; i++) {
				for (j = nearestToCenter.j; j >= 0; j--) {
					for (k = nearestToCenter.k; k < _l; k++) {
						if (validateElement(i, j, k))
							_validateInThisIteration++;
					}
				}
			}
			//4
			for (i = nearestToCenter.i; i >= 0; i--) {
				for (j = nearestToCenter.j; j >= 0; j--) {
					for (k = nearestToCenter.k; k < _l; k++) {
						if (validateElement(i, j, k))
							_validateInThisIteration++;
					}
				}
			}
			if (_l > 1) {
				//Only for 3D
				//5
				for (i = nearestToCenter.i; i < _n; i++) {
					for (j = nearestToCenter.j; j < _m; j++) {
						for (k = nearestToCenter.k; k >= 0; k--) {
							if (validateElement(i, j, k))
								_validateInThisIteration++;
						}
					}
				}
				//6
				for (i = nearestToCenter.i; i >= 0; i--) {
					for (j = nearestToCenter.j; j < _m; j++) {
						for (k = nearestToCenter.k; k >= 0; k--) {
							if (validateElement(i, j, k))
								_validateInThisIteration++;
						}
					}
				}
				//7
				for (i = nearestToCenter.i; i < _n; i++) {
					for (j = nearestToCenter.j; j >= 0; j--) {
						for (k = nearestToCenter.k; k >= 0; k--) {
							if (validateElement(i, j, k))
								_validateInThisIteration++;
						}
					}
				}
				//8
				for (i = nearestToCenter.i; i >= 0; i--) {
					for (j = nearestToCenter.j; j >= 0; j--) {
						for (k = nearestToCenter.k; k >= 0; k--) {
							if (validateElement(i, j, k))
								_validateInThisIteration++;
						}
					}
				}
			}
		}
		return _cvf;
	}
}