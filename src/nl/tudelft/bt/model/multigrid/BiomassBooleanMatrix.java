/*
 * File created originally on Mar 11, 2006
 */
package nl.tudelft.bt.model.multigrid;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import nl.tudelft.bt.model.DiscreteCoordinate;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_conditions.TubeBoundaryConditions;
import nl.tudelft.bt.model.particlebased.BiomassParticleContainer;
import nl.tudelft.bt.model.work.tube.MonoSpeciesTube;

/**
 * Implements a boolean matrix with the resolution of the finnest grid in
 * multigrid. The objects of this class allow queries on entries of the matrix
 * (e.g. if it is occuppied by biomass, by carrier or by either biomass or
 * carrier, or if a point is in the boundary biomass/liquid, etc.)
 * 
 * @author jxavier
 */
public class BiomassBooleanMatrix {
	private boolean[][][] _totalBiomass;

	private boolean[][][] _carrier;

	private boolean[][][] _carrierOrBiomass;

	// sizes of finest grid used in the multigrid
	private int _l;

	private int _m;

	private int _n;

	// other multigrid properties
	private BoundaryConditions _boundaryCondition;

	private int _order;

	// Tube2D
	private ArrayList _borderNodes;

	/**
	 * Part of the testing phase for tube reactors, implements a grid node in 2D
	 * coordinates (Tube2D)
	 * 
	 * @author jxavier
	 */
	private class GridNode2D extends DiscreteCoordinate {

		/**
		 * Initialize with k = 2;
		 * 
		 * @param i
		 * @param j
		 */
		public GridNode2D(int i, int j) {
			this.i = i;
			this.j = j;
			this.k = 1;

		}

		/**
		 * Make copy
		 * 
		 * @param i
		 * @param j
		 */
		public GridNode2D(GridNode2D n) {
			this(n.i, n.j);
		}

		/**
		 * @param movingDirection
		 * @return true if the next neighbor (according to the moving direction)
		 *         is carrier or biomass
		 */
		public boolean nextIsCarrierOrBiomass(Direction movingDirection) {
			return isCarrierOrBiomass2D(i + movingDirection.nextI(), j
					+ movingDirection.nextJ());
		}

		/**
		 * @param movingDirection
		 * @return the next neighbor (according to the moving direction)
		 */
		public GridNode2D next(Direction movingDirection) {
			return new GridNode2D(i + movingDirection.nextI(), j
					+ movingDirection.nextJ());
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Object#toString()
		 */
		public String toString() {
			return "" + i + "\t" + j + "\n";
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		public boolean equals(Object arg0) {
			GridNode2D n = ((GridNode2D) arg0);
			return (n.i == i) & (n.j == j);
		}
	}

	private class Direction {
		private String[] _dirNames = { "E", "NE", "N", "NW", "W", "SW", "S",
				"SE" };

		private int[] _dirI = { 1, 1, 0, -1, -1, -1, 0, 1 };

		private int[] _dirJ = { 0, 1, 1, 1, 0, -1, -1, -1 };

		private int _current;

		/**
		 * Initialize the direction for W (West)
		 */
		public Direction() {
			_current = 5;
		}

		/**
		 * rotate the present direction CCW (1/8 of full rotation)
		 */
		public void rotateCCW() {
			// rotate by incrementing but applying modulus of 8
			_current = (_current + 1) % _dirNames.length;
		}

		/**
		 * rotate the present direction CW (3/8 of full rotation)
		 */
		public void rotateCW() {
			// rotate by decrementing but applying modulus of 8
			_current = (8 + _current - 3) % _dirNames.length;
		}

		/**
		 * The i component of the move in this direction
		 */
		public int nextI() {
			return _dirI[_current];
		}

		/**
		 * The j component of the move in this direction
		 */
		public int nextJ() {
			return _dirJ[_current];
		}
	}

	/**
	 * Initialize the matrices using data from all the particulate species and
	 * using boundary conditions to detect position of carrier
	 */
	public BiomassBooleanMatrix() {
		// allocate space for matrices
		_totalBiomass = MultigridVariable
				.create3DBooleanMatrixWithFinnerResolution();
		_carrier = MultigridVariable
				.create3DBooleanMatrixWithFinnerResolution();
		_carrierOrBiomass = MultigridVariable
				.create3DBooleanMatrixWithFinnerResolution();
		_n = _totalBiomass.length;
		_m = _totalBiomass[0].length;
		_l = _totalBiomass[0][0].length;
		// get values from the multigrid class
		_order = MultigridVariable._order;
		_boundaryCondition = MultigridVariable._boundaryConditions;
	}

	/**
	 * get the information on all the biomass and carrier position
	 * 
	 * @param b
	 *            the array of all particulate species in the simulation
	 */
	public void upateBiomassConcentrations(ParticulateSpecies[] b) {
		for (int i = 0; i < _n; i++)
			for (int j = 0; j < _m; j++)
				for (int k = 0; k < _l; k++) {
					// the biomass
					_totalBiomass[i][j][k] = false;
					// chek if this is carrier
					if (_boundaryCondition.isCarrier(i, j, k)) {
						_carrier[i][j][k] = true;
					} else
						// if not
						for (int sp = 0; sp < b.length; sp++) {
							if (b[sp]._mg[_order - 1][i][j][k] > 0) {
								_totalBiomass[i][j][k] = true;
								break;
							}
						}
					// the carrier
					_carrier[i][j][k] = _boundaryCondition.isCarrier(i, j, k);
					// the carrier or biomass matrix
					_carrierOrBiomass[i][j][k] = (_totalBiomass[i][j][k] | _carrier[i][j][k]);
				}
	}

	/**
	 * Detect if a point is biomass of carrier and if it is located in the the
	 * limits (i.e. if it contains a neighbor which is not carrier or liquid)
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return true if (i, j, k) is biomass or carrier and is located on the
	 *         border
	 */
	private boolean isBorderPoint(int i, int j, int k) {
		if (!_carrierOrBiomass[i][j][k])
			return false;
		for (int i2 = -1; i < 1; i++)
			for (int j2 = -1; j < 1; j++)
				for (int k2 = -1; k < 1; k++) {
					try {
						if (!_carrierOrBiomass[i + i2][j + j2][k + k2])
							return true;
					} catch (IndexOutOfBoundsException e) {
						// do nothing
					}

				}
		return false;
	}

	/**
	 * Auxiliary method for determining if a poisition is occuppied
	 * 
	 * @param i
	 * @param j
	 * @return true if positon is occuppied, false otherwise or if position is
	 *         outside range
	 */
	public boolean isCarrierOrBiomass2D(int i, int j) {
		try {
			return _carrierOrBiomass[i][j][1];
		} catch (ArrayIndexOutOfBoundsException e) {
			// do nothing
		}
		// false is returned if location is outside the range
		return false;
	}

	/**
	 * Auxiliary method for determining if a poisition is occuppied
	 * 
	 * @param n
	 *            the poistion
	 * @return true if positon is occuppied
	 */
	private boolean isCarrierOrBiomass2D(GridNode2D n) {
		return isCarrierOrBiomass2D(n.i, n.j);
	}

	/**
	 * Auxiliary method for determining if a position is in border for 2D case
	 * 
	 * @param i
	 * @param j
	 * @return if position is ocuppied and on the border
	 */
	private boolean isBorderPoint2D(int i, int j) {
		return isBorderPoint(i, j, 1);
	}

	/**
	 * Find the points located in the liquib/biofilm border sorted
	 */
	// TODO change access control back to private
	public ArrayList findTheBorderPointsInTube() {
		// TODO delete after testing:
		BiomassParticleContainer.writeBooleanMatrixToFile(_carrierOrBiomass,
				MonoSpeciesTube.outputDirectory + "/carrierOrBiomass.txt");
		float[][] detectedPath = new float[_n][_m];
		for (int i = 0; i < _n; i++)
			for (int j = 0; j < _m; j++)
				if (_carrierOrBiomass[i][j][1])
					detectedPath[i][j] = 1;
				else
					detectedPath[i][j] = 0;
		int c = 1;
		// END

		// this is code still in test phase, a few tests follow to make sure
		// it is not inadvertelly used wrongly
		if (_l > 3) {
			throw new ModelRuntimeException(
					"Function should only be used in 2D");
		}
		if (!(MultigridVariable._boundaryConditions instanceof TubeBoundaryConditions)) {
			throw new ModelRuntimeException(
					"Function should only be used with tubes");
		}
		// start in center and move to right until first border point is
		// found
		int i = (int) (((float) _n) * 0.5f);
		int j = (int) (((float) _m) * 0.5f);
		//
		if (isCarrierOrBiomass2D(i, j)) {
			throw new ModelRuntimeException("Point used for starting the "
					+ "border detection is not empty");
		}
		// move left
		for (; i < _n + 1; i++) {
			if (i == _n)
				throw new ModelRuntimeException(
						"reached border of computations"
								+ " volume without finfing an occupied"
								+ " location");
			if (isCarrierOrBiomass2D(i, j))
				break;
		}
		// initialize the array and keep the first point to detect
		// the termination of the algorithm
		_borderNodes = new ArrayList();
		GridNode2D initial = new GridNode2D(i, j);
		_borderNodes.add(initial);
		// start the detection
		Direction movingDirection = new Direction();
		GridNode2D current = new GridNode2D(initial);
		try {
			while (true) {
				// TODO delete
				detectedPath[current.i][current.j] = ++c;
				// END
				// rotate until the next border point is found
				while (!current.nextIsCarrierOrBiomass(movingDirection)) {
					movingDirection.rotateCCW();
				}
				current = current.next(movingDirection);
				// check if cycle is finished
				if (current.equals(initial))
					break;
				// store the point
				_borderNodes.add(current);
				// rotate on the opposit direction
				movingDirection.rotateCW();
			}
		} catch (RuntimeException e) {
			// TODO delete this try/catch block
			e.printStackTrace();
		}
		// TODO delete
		writeMatrixToFile(detectedPath, "pathMatrix.txt");
		writePathToFile(_borderNodes, "path.txt");
		// END
		return _borderNodes;
	}

	/**
	 * Auxiliary method for debugging TODO delete after testing
	 * 
	 * @param path
	 * @param fn
	 */
	private void writePathToFile(ArrayList path, String fn) {
		// contruct string with the list
		String toWrite = "";
		for (Iterator iter = path.iterator(); iter.hasNext();) {
			GridNode2D element = (GridNode2D) iter.next();
			toWrite += element.toString();
		}
		// create the file and write matrix to it
		File f = new File(MonoSpeciesTube.outputDirectory + fn);
		try {
			java.io.FileWriter fr = new java.io.FileWriter(f);
			fr.write(toWrite);
			fr.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Auxiliary method for debugging TODO delete after testing
	 * 
	 * @param path
	 * @param fn
	 */
	private void writeMatrixToFile(float[][] m, String fn) {
		// contruct string with the list
		String toWrite = "";
		for (int j = 0; j < _m; j++) {
			for (int i = 0; i < _n; i++) {
				toWrite += "" + m[i][j] + "\t";
			}
			toWrite += "\n";
		}
		// create the file and write matrix to it
		File f = new File(MonoSpeciesTube.outputDirectory + fn);
		try {
			java.io.FileWriter fr = new java.io.FileWriter(f);
			fr.write(toWrite);
			fr.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public String toString() {
		return MultigridUtils.coreMatrixToString(_carrierOrBiomass);
	}
}