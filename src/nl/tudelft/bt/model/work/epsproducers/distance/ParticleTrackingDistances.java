package nl.tudelft.bt.model.work.epsproducers.distance;

import java.awt.Color;
import java.util.Collection;
import java.util.Iterator;

import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.util.ColorMaps;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Implements the mutator particle. A mutation will occur at a given time. The
 * change in trait will be defined by the biomass species. The mutante will
 * change its trait and will tag the remaining particles with the distance at
 * which they are to the mutant.
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Oct 18, 2006
 */
public class ParticleTrackingDistances extends BiomassParticle {
	private BiomassSpeciesTrackDistances _species;

	private float _distanceFocal;

	private static float _maximumDistanceToFocal = 0;;

	public ParticleTrackingDistances(BiomassSpeciesTrackDistances s) {
		super(s);
		_species = s;
	}

	/**
	 * @return the distance to the focal of ancestor cell
	 */
	public float getDistanceToFocal() {
		return _distanceFocal;
	}

	// @Override
	// public Color getColorCore() {
	// return ColorMaps.getBluescaleColor(_distanceFocal
	// / _maximumDistanceToFocal);
	// }

	/**
	 * Sets the distance to focal according to the distance to p. Also updated
	 * the maximum distance.
	 * 
	 * @param p
	 *            focal particle
	 */
	public void setDistanceTo(BiomassParticle p) {
		_distanceFocal = distanceTo(p);
		_maximumDistanceToFocal = ExtraMath.max(_distanceFocal,
				_maximumDistanceToFocal);
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		// TODO Auto-generated method stub
		return super.clone();
	}
}
