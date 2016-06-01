/*
 * Created on May 13, 2003
 */
package nl.tudelft.bt.model.util;

import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.exceptions.ModelException;

/**
 * Abstract class with some extra math functions
 * 
 * @author TU Delft
 */
public final class ExtraMath {
	/**
	 * the value of PI in float
	 */
	public static final float PI = 3.1416f;

	/**
	 * the value of 2*PI in float
	 */
	public static final float PI2 = 6.2832f;
	
	
	/**
	 * the value of log(2) in float
	 */
	public static final float LOG2 = 0.6931f;	
	

	/**
	 * Computes the logarothm of base 2
	 * 
	 * @param x
	 *            a number greater than 0.0
	 * @return the logarithm base 2 of x
	 */
	public static final double log2(double x) {
		return (Math.log(x) / LOG2);
	}

	/**
	 * Square
	 * 
	 * @param x
	 *            value to square
	 * @return x*x
	 */
	public static final int sq(int x) {
		return x * x;
	}

	/**
	 * Square
	 * 
	 * @param x
	 *            value to square
	 * @return x*x
	 */
	public static final float sq(float x) {
		return x * x;
	}

	/**
	 * Square root
	 * 
	 * @param x
	 *            value to square
	 * @return uses Math.sqrt and casts result to float
	 */
	public static final float sqrt(float x) {
		return (float) (Math.sqrt(x));
	}

	/**
	 * cube of a number
	 * 
	 * @param x
	 *            value to cube
	 * @return x*x*x
	 */
	public static final float cube(float x) {
		return x * x * x;
	}

	/**
	 * power of 2
	 * 
	 * @param x
	 * @return 2^x
	 */
	public static final int exp2(int x) {
		return (int) Math.pow(2, x);
	}

	/**
	 * The volume of a sphere with radius r
	 * 
	 * @param r
	 *            radius
	 * @return volume of sphere
	 */
	public static final float volumeOfASphere(float r) {
		return 4.1888f * r * r * r;
	}

	/**
	 * The volume of a cylinder with radius r and length l
	 * 
	 * @param r
	 *            radius
	 * @param l
	 *            length
	 * @return volume of cylinder
	 */
	public static final float volumeOfACylinder(float r, float l) {
		return 3.14f * r * r * l;
	}

	/**
	 * The area of circle with radius r
	 * 
	 * @param r
	 *            radius
	 * @return area of circle
	 */
	public static final float areaOfACircle(float r) {
		return 3.1416f * r * r;
	}

	/**
	 * Radius of a sphere with volume v
	 * 
	 * @param v
	 *            volume
	 * @return radius
	 */
	public static final float radiusOfASphere(float v) {
		return (float) Math.pow(0.2387 * v, 0.333);
	}

	/**
	 * Returns the radius of a cilinder with volume v and length l
	 * 
	 * @param v
	 * @return
	 * @throws ModelException
	 */
	public static final float radiusOfACilinder(float v, float l) {
		return (float) Math.sqrt(v / (3.14f * l));
	}

	/**
	 * Distance between 2 points
	 * 
	 * @param x1
	 * @param y1
	 * @param z1
	 * @param x2
	 * @param y2
	 * @param z2
	 * @return distance
	 */
	public static final float pointDistance(float x1, float y1, float z1,
			float x2, float y2, float z2) {
		return (float) Math.sqrt(sq(x1 - x2) + sq(y1 - y2) + sq(z1 - z2));
	}

	/**
	 * Distance between 2 points defined in 2D
	 * 
	 * @param x1
	 * @param y1
	 * @param x2
	 * @param y2
	 * @return distance
	 */
	public static final float pointDistance(float x1, float y1, float x2,
			float y2) {
		return (float) Math.sqrt(sq(x1 - x2) + sq(y1 - y2));
	}

	/**
	 * Distance between 2 points
	 * 
	 * @param p1
	 *            point 1
	 * @param p2
	 *            point 2
	 * @return distance
	 */
	public static final float pointDistance(ContinuousCoordinate p1,
			ContinuousCoordinate p2) {
		return pointDistance(p1.x, p1.y, p1.z, p2.x, p2.y, p2.z);
	}

	/**
	 * Perform gamma correction of value v: v^gamma check that v is in the [0,
	 * 1] range
	 * 
	 * @param v
	 * @param gamma
	 * @return v^gamma
	 * @throws ModelException
	 *             is v is not in the [0, 1] range
	 */
	public static float gammaCorrection(float v, float gamma)
			throws ModelException {
		if ((v < 0) | (v > 1))
			throw new ModelException("invalid v for gamma correction, v = " + v);
		return (float) Math.pow(v, gamma);
	}

	/**
	 * The maximum among 2 floats
	 * 
	 * @param a
	 * @param b
	 * @return the maximum among a and b
	 */
	public static float max(float a, float b) {
		return (a > b ? a : b);
	}

	/**
	 * The minimum among 2 floats
	 * 
	 * @param a
	 * @param b
	 * @return the minimum among a and b
	 */
	public static float min(float a, float b) {
		return (a > b ? b : a);
	}

	/**
	 * The maximum squate among 2 floats
	 * 
	 * @param a
	 * @param b
	 * @return the maximum square among a and b
	 */
	public static float maxSquare(float a, float b) {
		float a2 = sq(a);
		float b2 = sq(b);
		return (a2 > b2 ? a2 : b2);
	}

	/**
	 * the absolute value of a
	 * 
	 * @param a
	 * @return (a >= 0 ? a : -a)
	 */
	public static float abs(float a) {
		return (a >= 0 ? a : -a);
	}
}