package nl.tudelft.bt.model.work.will;

import nl.tudelft.bt.model.exceptions.NonMatchingNumberException;
import nl.tudelft.bt.model.multigrid.MultigridVariable;
import nl.tudelft.bt.model.reaction.Saturation;

/**
 * Implements saturation kinetics with minimum relative rate. i.e. Instead
 * of saturating from 0 -> 1 this saturates from f -> 1 for 0 <= f <= 1
 * @author wkc
 *
 */
public class SaturationWithFloor extends Saturation {
	
	private float _min;
	
	/**
	 * 
	 * @param c
	 * 			the dependent species
	 * @param k
	 * 			the saturation constant
	 * @param f
	 * 			the minimum relative value
	 */
	public SaturationWithFloor(MultigridVariable c, float k, float min) 
		throws NonMatchingNumberException {
		super(c, k);
		if (min <= 0 || min >= 1) {
			throw new NonMatchingNumberException(
					"min = " + min + " must be > 0 and < 1");
		}
		_min = min;
	}
	
	/**
	 * Gets the regular saturation value but returns the minimum fractional
	 * value if it is greater.
	 */
	@Override
	public float getValue() {
		float v = super.getValue();
		float out = (v >= _min ? v : _min);
		return out;
	}

}
