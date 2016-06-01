package nl.tudelft.bt.model;

import java.io.Serializable;

/**
 * Implements 3D vector of discrete spatial coordinates
 * 
 * @author João Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DiscreteCoordinate implements Serializable {
	public int i;
	public int j;
	public int k;

	/**
	 * Constructor (empty vector)
	 */
	public DiscreteCoordinate() {
		i = 0;
		j = 0;
		k = 0;
	}

	/**
	 * Constructor
	 * 
	 * @param l
	 *            depth coordinate
	 * @param m
	 *            horizontal coordinate
	 * @param n
	 *            vertical coordinate
	 */
	public DiscreteCoordinate(int n, int m, int l) {
		this.k = l;
		this.j = m;
		this.i = n;
	}

	public boolean equals(Object obj) {
		DiscreteCoordinate e = (DiscreteCoordinate) obj;
		return ((e.i == i) & (e.j == j) & (e.k == k));
	}
}