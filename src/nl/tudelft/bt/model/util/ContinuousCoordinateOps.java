package nl.tudelft.bt.model.util;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.Model;

/**
 * Some methods useful for vector math using ContinuousCoordinates 
 * within a model.
 * 
 * @author Will K. Chang
 * 
 * 2010
 */

public final class ContinuousCoordinateOps {
	
	/**
	 * Gets the length of a ContinuousCoordinate
	 * 
	 * @param v ContinuousCoordinate
	 * @return the scalar length (float)
	 */
	public static final float getLength(ContinuousCoordinate v) {
		float length = 0;
		length += Math.pow(v.x, 2);
		length += Math.pow(v.y, 2);
		length += Math.pow(v.z, 2);
		length = (float) Math.sqrt(length);
		return length;
	}
	
	/**
	 * Return a new ContinuousCoordinate of same direction but length 1
	 * 
	 * @param v ContinuousCoordinate to be normalized
	 * @return ContinuousCoordinate of same direction and length 1
	 */
	public static final ContinuousCoordinate getUnitVector(ContinuousCoordinate v) {
		ContinuousCoordinate w = new ContinuousCoordinate();
		float length = getLength(v);
		// divide each component by the length
		w.x = v.x / length;
		w.y = v.y / length;
		w.z = v.z / length;
		return w;
	}
	
	/**
	 * Multiplies a ContinuousCoordinate by a scalar and returns
	 * the result as a NEW ContinuousCoordinate
	 * 
	 * If you want to replace the original ContinuousCoordinate with
	 * the scaled version use ContinuousCoordinate.scale()
	 * 
	 * @param v	the ContinuousCoordinate to be scaled
	 * @param l	the scaling factor
	 * @return	the scaled ContinuousCoordinate
	 */
	public static final ContinuousCoordinate newScale(ContinuousCoordinate v, float l) {
		return new ContinuousCoordinate(v.x * l, v.y * l, v.z * l);
	}
	
	/**
	 * Sums an array of ContinuousCoordinates
	 * 
	 * @param vs an array of ContinuousCoordinates
	 * @return the summed vector (ContinuousCoordinate)
	 */
	public static final ContinuousCoordinate sum(ContinuousCoordinate[] vs) {
		ContinuousCoordinate w = new ContinuousCoordinate();
		for (int i = 0; i < vs.length; i++) {
			ContinuousCoordinate v = vs[i];
			w.add(v);
		}
		return w;
	}
	
	/**
	 * @return random ContinuousCoordinate of length 1
	 */
	public static final ContinuousCoordinate getRandomUnitVector() {
		float randomTheta = Model.model().getRandom() * ExtraMath.PI2;
		float randX = (float) Math.cos(randomTheta);
		float randY = (float) Math.sin(randomTheta);
		float randZ = 0;
		if (Model.model().getDimensionality() == 3) {
			// in 3D use theta as the azimuth, phi as the zenith
			float randomPhi = Model.model().getRandom() * ExtraMath.PI;
			randX = randX * (float) Math.sin(randomPhi);
			randY = randY * (float) Math.sin(randomPhi);
			randZ = (float) Math.cos(randomPhi);
		}
		return new ContinuousCoordinate(randX, randY, randZ);
	}
}
