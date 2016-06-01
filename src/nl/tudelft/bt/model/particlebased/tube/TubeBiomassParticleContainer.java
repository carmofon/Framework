/*
 * Created on Sep 15, 2004
 */
package nl.tudelft.bt.model.particlebased.tube;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.DiscreteCoordinate;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.detachment.cvf.ConnectedToCenterCvf;
import nl.tudelft.bt.model.detachment.cvf.ConnectedToTubeCvf;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.multigrid.BiomassBooleanMatrix;
import nl.tudelft.bt.model.multigrid.TubeVelocityField;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.BiomassParticleContainer;
import nl.tudelft.bt.model.particlebased.BiomassParticleSet;
import nl.tudelft.bt.model.particlebased.tube.detachment.DetachmentInducedByShear;
import nl.tudelft.bt.model.particlebased.tube.detachment.ParticleBasedDetachmentInducedByShear;
import nl.tudelft.bt.model.util.ColorMaps;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Particle container for tube reactor geometry. Used for simulating biofilm
 * growth in the interior of cyclindrical tubes (2D only)
 * 
 * @author jxavier
 */
public class TubeBiomassParticleContainer extends BiomassParticleContainer {
	private float _radiusOfTubularReactor;

	// this is used for convertion betweem radial coordinates and polar
	// coordinates
	private float _halfSystemSide;

	/**
	 * @param maximumRadius
	 * @param minimumRadius
	 * @param k
	 * @param f
	 */
	public TubeBiomassParticleContainer(float maximumRadius,
			float minimumRadius, float k, float f, float radiusOfTube) {
		super(maximumRadius, minimumRadius, k, f);
		_radiusOfTubularReactor = radiusOfTube;
		_halfSystemSide = _referenceSystemSide * 0.5f;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.particlebased.BiomassParticleContainer#createConnectedVolumeFiltrators()
	 */
	protected void createConnectedVolumeFilters() {
		_cvfForBiomass = new ConnectedToTubeCvf(_n, _m, this);
		_cvfForLiquid = new ConnectedToCenterCvf(_n, _m, _l);
	}

	/**
	 * Place a new particle in a random position near the center of the
	 * computational volume
	 * 
	 * @param s
	 *            species of particle to place
	 */
	public void placeInoculumParticle(BiomassSpecies s) {
		TubeBiomassParticle b = new TubeBiomassParticle(s);
		// coordinates to place the new particle (random
		float theta = _model.getRandom() * 2 * 3.1415f;
		float radius = _radiusOfTubularReactor - b.getRadius();

		// set center to center of computational volume
		b.setCenter(xConvertToCartesian(radius, theta), yConvertToCartesian(
				radius, theta), 0);
		// add to the bacteria list
		particleList.add(b);
		// shove to avoid overlapping
		spreadByShoving();
	}

	/**
	 * @return Returns the _radiusOfTubularReactor.
	 */
	public float getRadiusOfTubularReactor() {
		return _radiusOfTubularReactor;
	}

	/**
	 * Convert x,y coordinate to polar
	 * 
	 * @param x
	 * @param y
	 * @return radius coordinate
	 */
	public float rConvertToPolar(float x, float y) {
		x -= _halfSystemSide;
		y -= _halfSystemSide;
		// convert cartesian to polar
		return (float) Math.sqrt(ExtraMath.sq(x) + ExtraMath.sq(y));
	}

	/**
	 * Convert x,y coordinate to polar
	 * 
	 * @param x
	 * @param y
	 * @return theta coordinate
	 */
	public float thetaConvertToPolar(float x, float y) {
		x -= _halfSystemSide;
		y -= _halfSystemSide;
		// convert cartesian to polar, with theta in the range [0, 2 *pi]
		double theta = Math.atan2(y, x);
		theta = (theta + 2 * Math.PI) % (2 * Math.PI);
		return (float) theta;
	}

	/**
	 * Return half of the system size, to conver x,y coordinates into
	 * coordinates relatice to the tube center (0, 0)
	 * 
	 * @return the value of the system size *0.5
	 */
	public float getHalfSystemSide() {
		return _halfSystemSide;
	}

	/**
	 * Convert r, theta coordinate to cartesian
	 * 
	 * @param r
	 * @param theta
	 * @return the y coordinate
	 */
	public float xConvertToCartesian(float r, float theta) {
		return r * (float) Math.cos(theta) + _halfSystemSide;
	}

	/**
	 * Convert r, theta coordinate to cartesian
	 * 
	 * @param r
	 * @param theta
	 * @return the y coordinate
	 */
	public float yConvertToCartesian(float r, float theta) {
		return r * (float) Math.sin(theta) + _halfSystemSide;
	}

	/**
	 * Make sure that the particle with radius r is placed entirily inside the
	 * reactor, and if not move it until it is inside
	 * 
	 * @param rCoord
	 *            coordinate to place the particle
	 * @param r
	 *            radius of particle
	 * @return the new radius coordinate (updated or not)
	 */
	public float putInsideReactor(float rCoord, float r) {
		return (rCoord > _radiusOfTubularReactor - r ? _radiusOfTubularReactor
				- r : rCoord);
	}

	/**
	 * Checks if a grid node belongs to the border of the tube for
	 * initialization of the CVF (ConnectedToTubeCvf)
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return true if the point is located at the border, false otherwise
	 */
	public boolean belongsToTubeBorder(int i, int j, int k) {
		ContinuousCoordinate c = getGridElementCenter(i, j, k);
		if (rConvertToPolar(c.x, c.y) > _radiusOfTubularReactor)
			return false;
		// check all the grid element corners./
		// If any of those is outside the tube, then this is a border point
		// auxiliary variables
		float halfSide = shovingGridSide * 0.5f;
		float x, y;
		// NW corner
		x = c.x + halfSide;
		y = c.y - halfSide;
		if (rConvertToPolar(x, y) > _radiusOfTubularReactor)
			return true;
		// NE corner
		x = c.x + halfSide;
		y = c.y + halfSide;
		if (rConvertToPolar(x, y) > _radiusOfTubularReactor)
			return true;
		// SE corner
		x = c.x - halfSide;
		y = c.y + halfSide;
		if (rConvertToPolar(x, y) > _radiusOfTubularReactor)
			return true;
		// SW corner
		x = c.x - halfSide;
		y = c.y - halfSide;
		if (rConvertToPolar(x, y) > _radiusOfTubularReactor)
			return true;
		// none of the corners are outside,
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.particlebased.BiomassParticleContainer#isCarrier(int,
	 *      int, int)
	 */
	protected boolean isCarrier(int i, int j, int k) {
		ContinuousCoordinate c = getGridElementCenter(i, j, k);
		if (rConvertToPolar(c.x, c.y) > _radiusOfTubularReactor)
			return true;
		return false;
	}

	// /**
	// * This method overloads the removeDetachedBiomass of
	// * BiomassParticleContainer and does not use level set detachment (compute
	// * only if the kdet value is >0)
	// */
	// public void removeDetachedBiomass() throws ModelException {
	// if (_kdet > 0) {
	// // Run through the biomass particle list and erode each one
	// // find the border points
	// float timeStep = Model.model().getTimeStep();
	// for (int i = 0; i < _n; i++) {
	// for (int j = 0; j < _m; j++) {
	// for (int k = 0; k < _l; k++) {
	// // The time necessary for full erosion of a volume
	// // element at the border based on its mass and exposed
	// // area
	// // get the value of the shear stress
	// ContinuousCoordinate center = getGridElementCenter(i,
	// j, k);
	// float shear = TubeVelocityField.getInstance()
	// .getShearStress().getValueAt(center);
	// // compute the detachment speed
	// float rdet = _kdet * shear;
	// // time untill full detachment of voxel
	// float t = shovingGridSide / rdet;
	// BiomassParticleSet ps = shovingGrid[i][j][k];
	// if (ps != null) {
	// ps.resetIndex();
	// while (ps.hasNext()) {
	// //
	// float f = timeStep / t;
	// // erode particle and add detached biomass to
	// // container
	// Model.model().detachedBiomassContainer()
	// .addToErodedBiomassSeries(
	// ps.next().erode(f > 1 ? 1 : f));
	// }
	// }
	// }
	// }
	// }
	// }
	// }

	/**
	 * This method overloads the removeDetachedBiomass of
	 * BiomassParticleContainer and does not use level set detachment (compute
	 * only if the _detachmentShearStress value is >0)
	 */
	public void removeDetachedBiomass() throws ModelException {
		if (_detachmentHandler != null) {
			ParticleBasedDetachmentInducedByShear d = (ParticleBasedDetachmentInducedByShear) _detachmentHandler;
			if (!d.detachmentIsOff()) {
				// if (_detachmentHandler == null) {
				// // initialize the detachment handler on first call
				// _detachmentHandler = new
				// ParticleBasedDetachmentInducedByShear(
				// _detachmentShearStress);
				// }
				d.updateShearOfBorderPoints();
				// Don't iterate using iterator because list is changed
				// each time an element is removed
				for (int i = 0; i < particleList.size(); i++) {
					// get the next particle from the list
					BiomassParticle b = (BiomassParticle) particleList.get(i);
					// compute the shear stress acting on it
					float x = b.getCenterX();
					float y = b.getCenterY();

					// TODO also chck if particle is already on the border
					// (known
					// shear)
					float shear = d.computeShear(x, y);
					// float actualShear = shear
					// / ExtraMath.sqrt(b.getComposition()
					// .getFractionOfInert());
					// change color of particle, with color showing the shear
					b.overideColor(ColorMaps.getJetColor(shear
							/ d.getBreakingShearStress()));
					if (shear > d.getBreakingShearStress()) {
						b.setToDetachBySloughing();
					}
				}
				// detachMarkedAndOutOfBoundsBiomass();
				mapBiofilmElements();
				markParticlesForSloughing();
				detachMarkedAndOutOfBoundsBiomass();
			}
		}
	}

	// public void removeDetachedBiomass() throws ModelException {
	// if (_detachmentHandler != null) {
	// DetachmentInducedByShear d = (DetachmentInducedByShear)
	// _detachmentHandler;
	// if (!d.detachmentIsOff()) {
	// d.detach();
	// for (int i = 0; i < particleList.size(); i++) {
	// // get the next particle from the list
	// BiomassParticle b = (BiomassParticle) particleList.get(i);
	// if (d.checkIfBacteriaShouldDetach(b)) {
	// b.setToDetachBySloughing();
	// }
	// }
	// }
	// }
	// detachMarkedAndOutOfBoundsBiomass();
	// }

	/**
	 * Writes out text files with information about the detachmet algorithm
	 * 
	 * @param dir
	 */
	public void writeDetachmentInfo(String dir) {
		// write out several files (only if the handler has been initialized)
		if (_detachmentHandler != null)
			if (!_detachmentHandler.detachmentIsOff()) {
				// ((ParticleBasedDetachmentInducedByShear) _detachmentHandler)
				// .writeBorderWithValue(dir);
				// ((ParticleBasedDetachmentInducedByShear) _detachmentHandler)
				// .writeSegmentsToDetach(dir);
				// ((ParticleBasedDetachmentInducedByShear) _detachmentHandler)
				// .writeBiomassMatrix(dir);
			}
	}
}
