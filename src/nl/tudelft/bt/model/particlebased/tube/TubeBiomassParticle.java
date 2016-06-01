/*
 * Created on Sep 15, 2004
 */
package nl.tudelft.bt.model.particlebased.tube;


import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * Implements a biomass particle to be used in tube modelling overrides the out
 * of bounds method.
 * 
 * @author jxavier
 */
public class TubeBiomassParticle extends BiomassParticle {
	/**
	 * @param s
	 */
	public TubeBiomassParticle(BiomassSpecies s) {
		super(s);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.particlebased.BiomassParticle#outOfBounds()
	 */
	public boolean outOfBounds() {
		return false;
	}
	/**
	 * Set position of center.
	 * 
	 * @param x
	 *            height [micron]
	 * @param y
	 *            y [micron]
	 * @param z
	 *            z [micron]
	 */
	public void setCenter(float x, float y, float z) {
		//get the tube radius from the particle container
		TubeBiomassParticleContainer c = ((TubeBiomassParticleContainer) _particleContainer);
		// convert cartesian to polar
		float rCoord = c.rConvertToPolar(x, y);
		float thetaCoor = c.thetaConvertToPolar(x, y);
		// check if it is out of bounds
		rCoord = c.putInsideReactor(rCoord, _r);
		// reconvert to cartesian coordinates
		_centerX = c.xConvertToCartesian(rCoord, thetaCoor);
		_centerY = c.yConvertToCartesian(rCoord, thetaCoor);
		// tubular is always 2D case
		_centerZ = 0;
		setShovingMapPositionAndAddToShovingMap();
	}
}