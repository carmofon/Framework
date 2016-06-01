/*
 * Created on Apr 23, 2003
 */
package nl.tudelft.bt.model.particlebased;

import java.io.Serializable;

import nl.tudelft.bt.model.exceptions.BacteriaNotInSetException;

/**
 * Implements a set of biomass particles to which particles can be added or
 * removed. The choice of using this class instead of JAVA's native containers
 * was made to improve the speed of shoving iterations.
 * 
 * @author João Xavier (j.xavier@tnw.tudelft.nl)
 */
public class BiomassParticleSet implements Serializable {
	private BiomassParticle[] _set;
	private int _end;
	private int _index;
	// this was changed from 7 for working with the grid 256
	// clsm data
	static final private int N = 15;

	/**
	 * Creates a bacteriaSet with N slots for bacterium
	 */
	public BiomassParticleSet() {
		_end = 0;
		_index = 0;
		_set = new BiomassParticle[N];
	}

	/**
	 * Creates a bacteriaSet with N slots for bacterium with nneib*N slots for
	 * bacterium
	 * 
	 * @param nneib
	 *            factor to multiply by N
	 */
	public BiomassParticleSet(int nneib) {
		_end = 0;
		_index = 0;
		_set = new BiomassParticle[nneib * N];
	}

	/**
	 * Adds a bacterium to bacteria set, increases the size of the set if
	 * necessary
	 * 
	 * @param b
	 *            bacterium to add
	 */
	public void add(BiomassParticle b) {
		if (_end < _set.length) {
			_set[_end++] = b;
		} else {
			// increase storage size in case limit is reached
			BiomassParticle[] _set2 = new BiomassParticle[_set.length + 1];
			// copy values to new array
			for (int i = 0; i < _set.length; i++) {
				_set2[i] = _set[i];
			}
			// add the new value
			_set2[_end++] = b;
			// replace array with extended array
			_set = _set2;
		}
	}

	/**
	 * Removes a bacterium from the set
	 * 
	 * @param b
	 *            bacterium to remove
	 */
	public void remove(BiomassParticle b) throws BacteriaNotInSetException {
		for (int i = 0; i < _end; i++) {
			if (_set[i] == b) {
				// remove particle from entry
				_set[i] = null;
				// move remaining bacteria in list (if any) one
				// place backward
				for (int j = i; j < (_end - 1); j++) {
					_set[j] = _set[j + 1];
				}
				// remove reference from previous end
				_set[_end - 1] = null;
				// update end
				_end--;
				return;
			}
		}
		throw new BacteriaNotInSetException();
	}

	/**
	 * Adds all bacteria from set bs to this set
	 * 
	 * @param bs
	 *            set to get bacteria from
	 */
	public void addAll(BiomassParticleSet bs) {
		if (bs != null)
			for (int i = 0; i < bs._end; i++) {
				this.add(bs._set[i]);
			}
	}

	/**
	 * Reset current index for iterating through bacteria set
	 */
	public void resetIndex() {
		_index = 0;
	}

	/**
	 * Reset the end, which is the same as reseting the whole array. This does
	 * not effectivelly clean the memory space occupied by previous use of the
	 * set, but these will be overritten
	 */
	public void reset() {
		_end = 0;
		_index = 0;
	}

	/**
	 * Get bacteria in current index and increment index (for iteration
	 * purposes)
	 * 
	 * @return next bacterium in the set
	 */
	public BiomassParticle next() {
		return _set[_index++];
	}

	/**
	 * Boolen test
	 * 
	 * @return true is there is a bacterium next in the set and false otherwise
	 */
	public boolean hasNext() {
		return (_end > 0) && (_index < _end);
	}

	/**
	 * @return the number of particles currently in the set
	 */
	public int getNumberOfParticles() {
		return _end;
	}
}