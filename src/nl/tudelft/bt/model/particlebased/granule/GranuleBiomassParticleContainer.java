/*
 * Created on Sep 15, 2004
 */
package nl.tudelft.bt.model.particlebased.granule;

import java.util.Iterator;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.detachment.cvf.ConnectedToCenterCvf;
import nl.tudelft.bt.model.detachment.cvf.ConnectedToTopCvf;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.BiomassParticleContainer;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Particle container for granular geometry
 * 
 * @author jxavier
 */
public class GranuleBiomassParticleContainer extends BiomassParticleContainer {
	private ContinuousCoordinate _centerOfComputationalVolume;

	private ContinuousCoordinate _movementForCenterOfMass;

	/**
	 * @param maximumRadius
	 * @param minimumRadius
	 * @param k
	 * @param f
	 */
	public GranuleBiomassParticleContainer(float maximumRadius,
			float minimumRadius, float k, float f) {
		super(maximumRadius, minimumRadius, k, f);
		// determine center coordinates
		_centerOfComputationalVolume = new ContinuousCoordinate(
				_model.systemSize.x / 2f, _model.systemSize.y / 2f,
				_model.systemSize.z / 2f);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seenl.tudelft.bt.model.particlebased.BiomassParticleContainer#
	 * createConnectedVolumeFiltrators()
	 */
	protected void createConnectedVolumeFilters() {
		_cvfForBiomass = new ConnectedToCenterCvf(_n, _m, _l);
		_cvfForLiquid = new ConnectedToTopCvf(_n + 10, _m, _l);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.BiomassContainer#spreadByShoving()
	 */
	public void spreadByShoving() {
		super.spreadByShoving();
		if (_detachmentHandler != null)
			if (!_detachmentHandler.detachmentIsOff()) {
				// if detachment is being performed
				// move all particles so that center of mass moves
				// to the center of the computational volume
				_movementForCenterOfMass = computeMovementForCenterOfMass();
				for (Iterator iter = particleList.iterator(); iter.hasNext();) {
					BiomassParticle b = (BiomassParticle) iter.next();
					b.move(_movementForCenterOfMass.x,
							_movementForCenterOfMass.y,
							_movementForCenterOfMass.z);
				}
			}
	}

	/**
	 * @return The vector movement necessary to put center of mass in the center
	 *         of computational volume
	 */
	private ContinuousCoordinate computeMovementForCenterOfMass() {
		// determine granule center of mass
		float centerOfMassX = 0;
		float centerOfMassY = 0;
		float centerOfMassZ = 0;
		float totalMass = 0;
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle b = (BiomassParticle) iter.next();
			float m = b.getTotalMass();
			centerOfMassX += b.getCenterX() * m;
			centerOfMassY += b.getCenterY() * m;
			centerOfMassZ += b.getCenterZ() * m;
			totalMass += m;
		}
		centerOfMassX /= totalMass;
		centerOfMassY /= totalMass;
		centerOfMassZ /= totalMass;
		// determine movement vector
		float delx = centerOfMassX - _centerOfComputationalVolume.x;
		float dely = centerOfMassY - _centerOfComputationalVolume.y;
		float delz = centerOfMassZ - _centerOfComputationalVolume.z;
		return new ContinuousCoordinate(delx, dely, delz);
	}

	/**
	 * Place a new particle in a random position near the center of the
	 * computational volume
	 * 
	 * @param s
	 *            species of particle to place
	 */
	public void placeInoculumParticle(BiomassSpecies s) {
		placeInoculumParticleInsideRadius(s, -1);
	}

	/**
	 * @param s
	 *            species of particle to add
	 * @param r
	 *            the radius whithin which to place particle (if r < 0 place at
	 *            a distance from the center within the radius of the particle)
	 */
	public void placeInoculumParticleInsideRadius(BiomassSpecies s, float r) {
		BiomassParticle b = s.createBiomassParticle();
		// set center to center of computational volume
		b.setCenter(_centerOfComputationalVolume.x,
				_centerOfComputationalVolume.y, _centerOfComputationalVolume.z);
		// random placement
		//float phi = ExtraMath.PI * 0.5f;
		float phi = _model.getRandom() * ExtraMath.PI;
		float theta = _model.getRandom() * ExtraMath.PI2;
		float radius = (r > 0 ? _model.getRandom() * r : b.getRadius());
		// // create a movement
		float delx = radius * (float) (Math.cos(theta) * Math.sin(phi));
		float dely = radius * (float) (Math.sin(theta) * Math.sin(phi));
		float delz = radius * (float) (Math.cos(phi));
		// // //distribute randomly inside a square
		// float delx = ((_model.getRandom() * 2 * r) - r);
		// float dely = ((_model.getRandom() * 2 * r) - r);
		// float delz = ((_model.getRandom() * 2 * r) - r);
		// move the particle to final position
		b.move(delx, dely, delz);
		// add to the bacteria list
		particleList.add(b);
		// shove to avoid overlapping
		spreadByShoving();
	}
	
	
	
	/**
	 * Places a BiomassParticle at distance r from center of space
	 * 
	 * @param s
	 * 			species fo particle to add
	 * @param r
	 * 			distance from center at which to add particle
	 */
	public void placeInoculumParticleAtRadius(BiomassSpecies s, float r) {
		// create the BiomassParticle
		BiomassParticle b = s.createBiomassParticle();
		
		// set the particle center to simulation space center
		b.setCenter(_centerOfComputationalVolume.x, _centerOfComputationalVolume.y, _centerOfComputationalVolume.z);
		
		// generate random theta, phi
		float theta = _model.getRandom() * ExtraMath.PI2;
		float phi = (_model.getDimensionality() == 3 ? 
				(float) _model.getRandom() * ExtraMath.PI : (float) 0.5 * ExtraMath.PI);
		
		// calculate movement components and move
		float delx = r * (float) (Math.cos(theta) * Math.sin(phi));
		float dely = r * (float) (Math.sin(theta) * Math.sin(phi));
		float delz = r * (float) (Math.cos(phi));
		b.move(delx, dely, delz);
		
		// add to the particle list
		particleList.add(b);
		
		// shove
		spreadByShoving();
	}

	/**
	 * Places a BiomassParticle at distance r and angle a from center of space
	 * 
	 * @param s
	 * 			species fo particle to add
	 * @param r
	 * 			distance from center at which to add particle
	 */
	public void placeInoculumParticleAtRadiusandAngle(BiomassSpecies s, float r, float a) {
		// create the BiomassParticle
		BiomassParticle b = s.createBiomassParticle();
		
		// set the particle center to simulation space center
		b.setCenter(_centerOfComputationalVolume.x, _centerOfComputationalVolume.y, _centerOfComputationalVolume.z);
		
		// generate random theta, phi
		float theta = (float) ((_model.getRandom()*b.getRadius()*1e-2 +a)*ExtraMath.PI2);
		float phi = (float) (_model.getDimensionality() == 3 ? 
				(float) ((_model.getRandom()*b.getRadius()*1e-2 +a)*ExtraMath.PI) : (float) 0.5 * ExtraMath.PI);
		
		// calculate movement components and move
		float delx = r * (float) (Math.cos(theta) * Math.sin(phi));
		float dely = r * (float) (Math.sin(theta) * Math.sin(phi));
		float delz = r * (float) (Math.cos(phi));
		b.move(delx, dely, delz);
		
		// add to the particle list
		particleList.add(b);
		
		// shove
		spreadByShoving();
	}
	/**
	 * @return Returns the _centerOfComputationalVolume.
	 */
	public ContinuousCoordinate get_centerOfComputationalVolume() {
		return _centerOfComputationalVolume;
	}

	/**
	 * Compute the granule run length along the x direction
	 * 
	 * @return the run length of the granule
	 */
	public float computeRunLengthX() {
		// iterate particles and find the highest and the lowest
		float highest = 0;
		float lowest = Float.POSITIVE_INFINITY;
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle b = (BiomassParticle) iter.next();
			float x = b.getCenterX();
			float r = b.getRadius();
			highest = (highest < (x + r) ? (x + r) : highest);
			lowest = (lowest > (x - r) ? (x - r) : lowest);
		}
		return highest - lowest;
	}

	/**
	 * Compute the granule run length along the y direction
	 * 
	 * @return the run length of the granule
	 */
	public float computeRunLengthY() {
		// iterate particles and find the highest and the lowest
		float highest = 0;
		float lowest = Float.POSITIVE_INFINITY;
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle b = (BiomassParticle) iter.next();
			float y = b.getCenterY();
			float r = b.getRadius();
			highest = (highest < (y + r) ? (y + r) : highest);
			lowest = (lowest > (y - r) ? (y - r) : lowest);
		}
		return highest - lowest;
	}

	/**
	 * Compute the granule run length along the z direction
	 * 
	 * @return the run length of the granule
	 */
	public float computeRunLengthZ() {
		// iterate particles and find the highest and the lowest
		float highest = 0;
		float lowest = Float.POSITIVE_INFINITY;
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle b = (BiomassParticle) iter.next();
			float z = b.getCenterZ();
			float r = b.getRadius();
			highest = (highest < (z + r) ? (z + r) : highest);
			lowest = (lowest > (z - r) ? (z - r) : lowest);
		}
		return highest - lowest;
	}

	public float getMaximumRadiusOfGranule() {
		// iterate particles and find the highest and the lowest
		float highestRadius = 0;
		float centerX = Model.model().systemSize.x * 0.5f;
		float centerY = Model.model().systemSize.y * 0.5f;
		float centerZ = Model.model().systemSize.z * 0.5f;
		for (Iterator<BiomassParticle> iter = particleList.iterator(); iter
				.hasNext();) {
			BiomassParticle b = iter.next();
			float x = b.getCenterX() - centerX;
			float y = b.getCenterY() - centerY;
			float z = b.getCenterZ() - centerZ;
			float r = b.getRadius()
					+ ExtraMath.sqrt(ExtraMath.sq(x) + ExtraMath.sq(y)
							+ ExtraMath.sq(z));
			highestRadius = (highestRadius < (z + r) ? (z + r) : highestRadius);
		}
		return highestRadius;
	}

}