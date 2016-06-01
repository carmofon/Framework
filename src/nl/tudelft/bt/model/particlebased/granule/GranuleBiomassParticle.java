/*
 * Created on Sep 15, 2004
 */
package nl.tudelft.bt.model.particlebased.granule;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.multigrid.boundary_layers.SphericalDilationBoundaryLayer;
import nl.tudelft.bt.model.particlebased.BiomassParticle;

/**
 * Implements a biomass particle to be used in granule modelling overrides the
 * out of bounds method.
 * 
 * @author jxavier
 */
public class GranuleBiomassParticle extends BiomassParticle {
	private static boolean _sheedingOn = false;

	/**
	 * @param s
	 */
	public GranuleBiomassParticle(BiomassSpecies s) {
		super(s);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.particlebased.BiomassParticle#outOfBounds()
	 */
	public boolean outOfBounds() {
		if (_sheedingOn)
			return ((SphericalDilationBoundaryLayer) Model.model()
					.getBoundaryLayer()).sphericalOutOfBounds(_centerX,
					_centerY, _centerZ);
		return false;
	}

	/**
	 * Set position of center. All borders are cyclic
	 * 
	 * @param x
	 *            x [micron]
	 * @param y
	 *            y [micron]
	 * @param z
	 *            z [micron]
	 */
	public void setCenter(float x, float y, float z) {
		// Implement periodic boundaries
		float Lxx = _model.systemSize.x;
		if (x >= Lxx)
			_centerX = x - Lxx;
		else if (x < 0)
			_centerX = Lxx + x;
		else
			_centerX = x;
		float Lyy = _model.systemSize.y;
		if (y >= Lyy)
			_centerY = y - Lyy;
		else if (y < 0)
			_centerY = Lyy + y;
		else
			_centerY = y;
		float Lzz = _model.systemSize.z;
		if (Lzz == 0)
			// 2D case
			_centerZ = 0;
		else if (z >= Lzz)
			_centerZ = z - Lzz;
		else if (z < 0)
			_centerZ = Lzz + z;
		else
			_centerZ = z;
		setShovingMapPositionAndAddToShovingMap();
	}

	
	/**
	 * Turns cell shedding off (no granule limit)
	 */
	public static void turnSheddingOn() {
		_sheedingOn = true;
	}
}