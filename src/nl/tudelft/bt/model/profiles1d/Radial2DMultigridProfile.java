/*
 * File created originally on Jul 21, 2005
 */
package nl.tudelft.bt.model.profiles1d;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.MultigridVariable;

/**
 * Implements the a radial profile for a multigrid variable, with center placed
 * at the center of the computational volume. Only works in 2D systems
 * 
 * @author jxavier
 */
public class Radial2DMultigridProfile implements Profile {
	private static final int NANGLES = 20;

	private static final int NSTEPS = 50;

	private MultigridVariable _variable;

	private float[] _radius = new float[NSTEPS];

	private float[] _averageValue = new float[NSTEPS];

	private int[][] _iValues = new int[NANGLES][NSTEPS];

	private int[][] _jValues = new int[NANGLES][NSTEPS];

	/**
	 * Initilize the radial profile
	 * 
	 * @param v
	 *            variable to make profiles from
	 */
	public Radial2DMultigridProfile(MultigridVariable v) {
		if (Model.model().getDimensionality() != 2) {
			throw new ModelRuntimeException(
					"Trying to create Radial2DMultigridProfile"
							+ " but system is not 2D");
		}
		_variable = v;
		// initialize _iValues, _jValues
		//get the gridSize
		int[] c = MultigridVariable.getFinnerGridSize();
		//check if square system
		if (c[0] != c[1])
			throw new ModelRuntimeException(
					"Trying to create Radial2DMultigridProfile"
							+ " but system is not square (" + c[0] + "," + c[1]
							+ ")");
		int n = c[0];
		//initilize an array with radius
		for (int i = 0; i < _radius.length; i++) {
			// radius does not go al the way to end of padding (n-1)
			_radius[i] = (float) i * (float) (n - 1) * 0.5f
					/ (float) (NSTEPS - 1);
		}
		//initilize an array with theta
		float[] theta = new float[NANGLES];
		for (int i = 0; i < theta.length; i++) {
			theta[i] = (float) i * (float) (2f * 3.14f) / (float) (NANGLES - 1);
		}
		//convert radius and theta into cartesian coordinates
		for (int i = 0; i < NANGLES; i++) {
			for (int j = 0; j < NSTEPS; j++) {
				// compute x
				float x = _radius[j] * (float) Math.cos(theta[i]) + n * 0.5f;
				// snap to grid
				_iValues[i][j] = (int) Math.floor(x);
				// compute y
				float y = _radius[j] * (float) Math.sin(theta[i]) + n * 0.5f;
				// snap to grid
				_jValues[i][j] = (int) Math.floor(y);
				//TODO test code
				if ((_iValues[i][j] >= n) | (_jValues[i][j] >= n)) {
					int lixo = 0;
				}
				//END

			}
		}
		//convert the radius into metric units used in the model
		float voxelSide = Model.model().referenceSystemSide / (float) (n - 2);
		for (int i = 0; i < NSTEPS; i++) {
			_radius[i] *= voxelSide;
		}
	}

	/**
	 * Compute the radial profile by performing averages at different radial
	 * increments
	 */
	public void computeProfile() {
		for (int j = 0; j < NSTEPS; j++) {
			_averageValue[j] = 0;
			for (int i = 0; i < NANGLES; i++) {
				_averageValue[j] += _variable.getValueAt(_iValues[i][j],
						_jValues[i][j], 1);
			}
			_averageValue[j] /= NANGLES;
		}
	}

	/**
	 * @return a string with a table with the values in the radial profile
	 */
	public String getFormatedTable() {
		String s = "radius [L]\t" + _variable.getName() + "\n";
		for (int i = 0; i < NSTEPS; i++) {
			s += _radius[i] + "\t" + _averageValue[i] + "\n";
		}
		return s;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.profiles1d.Profile#getName()
	 */
	public String getName() {
		return _variable.getName() + "_radial_profile";
	}

	/**
	 * get the average value at radius r from the profile
	 * 
	 * @param r
	 * @return value at r
	 */
	public float getValueAt(float r) {
		int i;
		for (i = 0; i < NSTEPS - 1; i++) {
			if (_radius[i + 1] > r)
				break;
		}
		if (i == NSTEPS - 1)
			throw new ModelRuntimeException(
					"Error geting radius from profile\n" + "radius " + r
							+ " is out of range");
		return _averageValue[i];
	}
}