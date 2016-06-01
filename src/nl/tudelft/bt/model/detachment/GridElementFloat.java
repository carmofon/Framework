/*
 * Created on 30-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.detachment;

import nl.tudelft.bt.model.DiscreteCoordinate;

/**
 * Implements a reference to a grid element
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class GridElementFloat extends DiscreteCoordinate {
	/**
	 * Sorts the grid elements according to the T value
	 */
	public static class TValueComparator implements java.util.Comparator {
		public int compare(Object b1, Object b2) {
			float z1 = ((GridElementFloat) b1).getValue();
			float z2 = ((GridElementFloat) b2).getValue();
			return (z1 > z2 ? 1 : -1);
		}
	}

	private float[][][] _values;

	/**
	 * Initialize a new grid element
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @param v
	 */
	public GridElementFloat(float[][][] values, int i, int j, int k) {
		this.i = i;
		this._values = values;
		this.j = j;
		this.k = k;
	}

	/**
	 * @return the levelset value for this grid element
	 */
	public float getValue() {
		return _values[i][j][k];
	}
}