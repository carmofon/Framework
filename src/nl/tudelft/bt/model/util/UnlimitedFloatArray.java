/*
 * Created on 28-jan-2004 by Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
package nl.tudelft.bt.model.util;

import java.io.Serializable;

import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;

/**
 * Implements an increasing size container for floats, with iterating
 * capabilities
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class UnlimitedFloatArray implements Serializable{
	private float[] _set;
	private int _end;
	private int _index;
	static final private int N = 100;

	/**
	 * Creates a bacteriaSet with N slots for bacterium
	 */
	public UnlimitedFloatArray() {
		_end = 0;
		_index = 0;
		_set = new float[N];
	}

	/**
	 * Creates an array with N slots for
	 * 
	 * @param n
	 *            number of slots
	 */
	public UnlimitedFloatArray(int n) {
		_end = 0;
		_index = 0;
		_set = new float[n];
	}

	/**
	 * Adds a float to array, increases the size of the set if necessary
	 * 
	 * @param f
	 *            float to add
	 */
	public void add(float f) {
		if (_end < _set.length) {
			_set[_end++] = f;
		} else {
			// increase storage size in case limit is reached
			float[] _set2 = new float[_set.length + 100];
			// copy values to new array
			for (int i = 0; i < _set.length; i++) {
				_set2[i] = _set[i];
			}
			// add the new value
			_set2[_end++] = f;
			// replace array with extended array
			_set = _set2;
		}
	}

	/**
	 * Reset current index for iterating through bacteria set
	 */
	public void resetIndex() {
		_index = 0;
	}

	/**
	 * Get number in current index and increment index (for iteration purposes)
	 * 
	 * @return next number in the set
	 */
	public float next() {
		return _set[_index++];
	}

	/**
	 * Boolen test
	 * 
	 * @return true is there is a float next in the set and false otherwise
	 */
	public boolean hasNext() {
		return (_end > 0) && (_index < _end);
	}

	/**
	 * @return the array of data (just relevant slora, not the empty ones)
	 */
	public float[] getArray() {
		float[] array = new float[_end];
		for (int i = 0; i < _end; i++) {
			array[i] = _set[i];
		}
		return array;
	}

	/**
	 * @return the last value in the array
	 */
	public float getLastValue() throws ModelRuntimeException {
		if (_end == 0)
			throw new ModelRuntimeException(
				"trying to get value" + " from empty UnlimitedFloatArray");
		return _set[_end - 1];
	}

	/**
	 * @param v
	 * @return the index of the first occurrence of v in this set; returns -1
	 *         if v is not found.
	 */
	public int indexOf(float v) {
		for (int i = 0; i < _end; i++) {
			if (_set[i] == v)
				return i;
		}
		return -1;
	}

	/**
	 * Adds value v to the entry in position i of the set
	 * 
	 * @param i
	 * @param v
	 * @throws ModelException
	 *             if i exceeds set limit
	 */
	public void addValueToEntry(int i, float v) throws ModelException {
		if (i < _end)
			_set[i] += v;
		else
			throw new ModelException(
				"trying to add above" + " UnlimitedFloatArray limit");
	}

	public float getValue(int i) throws ModelRuntimeException {
		if (i >= _end)
			throw new ModelRuntimeException("out of bounds");
		return _set[i];
	}

	/**
	 * @return Returns the _end.
	 */
	public int getSize() {
		return _end;
	}

}
