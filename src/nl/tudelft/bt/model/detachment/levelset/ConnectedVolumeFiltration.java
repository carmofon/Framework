/*
 * Created on 30-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.detachment.levelset;

import java.util.ArrayList;
import java.util.Iterator;

import nl.tudelft.bt.model.DiscreteCoordinate;
import nl.tudelft.bt.model.exceptions.ModelException;

/**
 * Implements connected volume filtration operation
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ConnectedVolumeFiltration {
	private int _n;
	private int _m;
	private int _l;
	private boolean[][][] _totalVolume;
	private boolean[][][] _initializer;
	private ArrayList _totalVolumeElements;
	private ArrayList _close;

	/**
	 * Implements a reference to a grid element
	 * 
	 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
	 */
	public class GridElementBoolean extends DiscreteCoordinate {
		/**
		 * @param n
		 * @param m
		 * @param l
		 */
		public GridElementBoolean(int n, int m, int l) {
			super(n, m, l);
		}

		/**
		 * @return a list of all neighbors valid accoring to _totalVolume
		 */
		protected ArrayList getNeighbors() {
			ArrayList n = new ArrayList();
			// bottom neighbor
			if (i > 0)
				if (_totalVolume[i - 1][j][k])
					n.add(new GridElementBoolean(i - 1, j, k));
			// y-side neighbors:
			if (_totalVolume[i][j < _m - 1 ? j + 1 : 0][k])
				n.add(new GridElementBoolean(i, j < _m - 1 ? j + 1 : 0, k));
			if (_totalVolume[i][j > 0 ? j - 1 : _m - 1][k])
				n.add(new GridElementBoolean(i, j > 0 ? j - 1 : _m - 1, k));
			// z-side neighbors:
			if (_l > 1) {
				if (_totalVolume[i][j][k < _n - 1 ? k + 1 : 0])
					n.add(new GridElementBoolean(i, j, k < _n - 1 ? k + 1 : 0));
				if (_totalVolume[i][j][k > 0 ? k - 1 : _n - 1])
					n.add(new GridElementBoolean(i, j, k > 0 ? k - 1 : _n - 1));;
			}
			// top neighbor:
			if (i < _n - 1)
				if (_totalVolume[i + 1][j][k])
					n.add(new GridElementBoolean(i + 1, j, k));
			return n;
		}
	}

	/**
	 * Creates a list with all true grid elements in u
	 * 
	 * @param u
	 * @return list of all true elements in u
	 */
	private ArrayList getElements(boolean[][][] u) {
		int n = u.length;
		int m = u[0].length;
		int l = u[0][0].length;
		ArrayList e = new ArrayList();
		for (int i = 0; i < u.length; i++) {
			for (int j = 0; j < u[0].length; j++) {
				for (int k = 0; k < u[0][0].length; k++) {
					if (u[i][j][k])
						e.add(new GridElementBoolean(i, j, k));
				}
			}
		}
		return e;
	}

	/**
	 * Initialize an instance of connected volume filtration
	 * 
	 * @param totalVolume
	 *            the volume to filter
	 * @param initializer
	 *            the initializing volume
	 */
	public ConnectedVolumeFiltration(boolean[][][] totalVolume,
			boolean[][][] initializer) throws ModelException {
		_n = totalVolume.length;
		_m = totalVolume[0].length;
		_l = totalVolume[0][0].length;
		if ((initializer.length != _n) | (initializer[0].length != _m)
				| (initializer[0][0].length != _l))
			throw new ModelException("ConnectedVolumeFiltration "
					+ "initialized with incompatible matrices");
		_totalVolume = totalVolume;
		_initializer = initializer;
		// get a colection of points represent the total volume and initializer
		_totalVolumeElements = getElements(_totalVolume);
		_close = getElements(_initializer);
		if (!_totalVolumeElements.containsAll(_close))
			throw new ModelException("ConnectedVolumeFiltration: "
					+ "not all elements in initializer are in total volumr");
	}

	/**
	 * @return the connected volume elements as a list of DiscreteCoordinate
	 *         elements
	 */
	public boolean[][][] getConnectVolumeElements() {
		ArrayList connectedVolumeElements = new ArrayList();
		while (_close.size() > 0) {
			// get all neighbors of close that belong to total volume but not
			// to _close nor to connectedvolume
			GridElementBoolean trial = (GridElementBoolean) _close.get(0);
			_close.remove(trial);
			connectedVolumeElements.add(trial);
			ArrayList nb = trial.getNeighbors();
			nb.removeAll(connectedVolumeElements);
			nb.removeAll(_close);
			_close.addAll(nb);
		}
		boolean[][][] cvf = new boolean[_n][_m][_l];
		for (Iterator iter = connectedVolumeElements.iterator(); iter.hasNext();) {
			DiscreteCoordinate c = (DiscreteCoordinate) iter.next();
			cvf[c.i][c.j][c.k] = true;
		}
		return cvf;
	}

	public static void computeCVF(boolean[][][] totalVolume, boolean[][][] cvf) {
		int n = totalVolume.length;
		int m = totalVolume[0].length;
		int l = totalVolume[0][0].length;
		// the loop
		int validateInThisIteration = 1;
		boolean validateInThisIScan = false;
		int i, j, k; //pre-allocate for optimization
		while (validateInThisIteration > 0) {
			validateInThisIteration = 0;
			for (i = n-2; i >= 0; i--) {
				validateInThisIScan = false;
				for (j = 0; j < m; j++) {
					for (k = 0; k < l; k++) {
						if (!cvf[i][j][k] & totalVolume[i][j][k]) {
							// check if any of the neighbors is valid
							// bottom neighbor
							if (((i > 0) ? cvf[i - 1][j][k] : false)
									// y-side neighbors:
									| (cvf[i][j < m - 1 ? j + 1 : 0][k])
									| (cvf[i][j > 0 ? j - 1 : m - 1][k]) // z-side
									// neighbors:
									| ((l > 1) & (cvf[i][j][k < l - 1
											? k + 1
											: 0]))
									| ((l > 1) & (cvf[i][j][k > 0
											? k - 1
											: l - 1]))
									// top neighbor:
									| ((i < n - 1) & (cvf[i + 1][j][k]))) {
								cvf[i][j][k] = true;
								validateInThisIteration++;
								validateInThisIScan = true;
							}
						}
					}
				}
				if (!validateInThisIScan)
					// break iteration if no element was validated in this scan
					break;
			}
		}
	}
}