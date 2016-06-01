package nl.tudelft.bt.model.particlebased;

import java.awt.Color;
import java.io.Serializable;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.BiomassSpecies.Composition;
import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.exceptions.BacteriaNotInSetException;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * Hard sphere biomass component for shoving spreading mechanism (Individual
 * Based Modeling). Implements all behavior of a cell modeled as a hard sphere.
 * 
 * @author João Xavier (j.xavier@tnw.tudelft.nl)
 */
public class BiomassParticle implements Cloneable, Serializable {
	protected BiomassSpecies _biomassSpecies;

	protected BiomassSpecies.Composition _composition;

	protected float _r;

	protected float _centerX;

	protected float _centerY;

	protected float _centerZ;

	private int _shovL;

	private int _shovM;

	private int _shovN;

	private float _divisionTheta;

	private float _divisionPhi;

	protected float _radiusChangeInThisIteration;

	private ContinuousCoordinate _moveByPressureVector;

	private boolean _willDetachByErosion;

	private boolean _willDetachBySloughing;

	protected Model _model = Model.model();

	// handle to the particle container, initialized uppon creation of each
	// particle
	protected BiomassParticleContainer _particleContainer;

	private boolean _overideColor;

	private Color _color;

	/**
	 * Sorts the bacteria list so that objects in first planes are drawn later.
	 * 
	 * @author TU Delft
	 */
	public static class DepthComparator implements java.util.Comparator {
		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		public int compare(Object b1, Object b2) {
			float z1 = ((BiomassParticle) b1)._centerZ;
			float z2 = ((BiomassParticle) b2)._centerZ;
			return (z1 < z2 ? 1 : -1);
		}
	}

	/**
	 * Sorts the bacteria list so Bacterium in outer layers ares shoved first.
	 * 
	 * @author TU Delft
	 */
	public static class HeigthComparator implements java.util.Comparator {
		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		public int compare(Object b1, Object b2) {
			return (((BiomassParticle) b1)._centerX > ((BiomassParticle) b2)._centerX ? 1
					: -1);
		}
	}

	/**
	 * Sorts the bacteria list so Bacterium in outer layers ares shoved last.
	 * 
	 * @author TU Delft
	 */
	public static class InvertedHeigthComparator implements
			java.util.Comparator {
		/*
		 * (non-Javadoc)
		 * 
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		public int compare(Object b1, Object b2) {
			return (((BiomassParticle) b1)._centerX > ((BiomassParticle) b2)._centerX ? -1
					: 1);
		}
	}

	/**
	 * Check is bacteria is of a certain species
	 * 
	 * @param s
	 *            species
	 * @return true if bacteria is of species s, false otherwise
	 */
	public boolean isOfSpecies(BiomassSpecies s) {
		return _biomassSpecies.equals(s);
	}

	/**
	 * @return Returns the biomass species of the partricle.
	 */
	public BiomassSpecies getBiomassSpecies() {
		return _biomassSpecies;
	}

	/**
	 * Create a new bacteria with default compostion mass defined in the
	 * BiomassSpecies instance
	 * 
	 * @param s
	 *            the biomass species of particle to create
	 */
	public BiomassParticle(BiomassSpecies s) {
		_biomassSpecies = s;
		// create composition
		_composition = _biomassSpecies.newParticleComposition();
		// initialize position stuff
		_centerX = 0;
		_centerY = 0;
		_centerZ = 0;
		_shovL = 0;
		_shovM = 0;
		_shovN = 0;
		_radiusChangeInThisIteration = 0;
		_willDetachByErosion = false;
		_willDetachBySloughing = false;
		// initialize the radius
		_r = _composition.computeRadius();
		// get the reference to the particle container
		_particleContainer = (BiomassParticleContainer) Model.model().biomassContainer;
	}

	/**
	 * Grow bacteria acording to the determined local growth rate. Update the
	 * bacterium composition and radius
	 * 
	 * @param t
	 *            time step
	 * @return the total biomass produced/consumed
	 */
	public float grow(float t) {
		// get - volume before growth
		_radiusChangeInThisIteration = -_r;
		// create a copy of center
		ContinuousCoordinate c = getCenter();
		// update composition
		float m = _composition.grow(t, getCenter());
		// allow reseting of center position for flexibility in the
		// growth (e.g. to implement motility)
		_centerX = c.x;
		_centerY = c.y;
		_centerZ = c.z;
		// float m = _composition.grow(t, getCenter());
		// update radius
		_r = _composition.computeRadius();
		// subtract radius to determine radius change
		_radiusChangeInThisIteration += _r;
		return m;
	}

	/**
	 * Decrease mass of every component by factor f and returns a new instance
	 * of BiomassSpecies.Composition with the detached biomass
	 * 
	 * @param f
	 *            must be 0 <f <1
	 * @return detached biomass
	 */
	public BiomassSpecies.Composition erode(float f) {
		// get - volume before growth
		float radiusBeforeErosion = _r;
		// create a copy of center
		ContinuousCoordinate c = getCenter();
		// update composition
		// criterium changed to 1.1 instead of 1 as value 1 can occur
		// (although chances of this happening are very low) due to precision
		// of float type
		if ((f < 0) | (f >= 1.1) | (Float.isNaN(f)))
			throw new ModelRuntimeException("illegal value f = " + f
					+ ": f must be 0<f<1");
		// in case f is 1 or (slightly) higher, put value 0.999999f
		f = (f >= 1 ? 0.999999f : f);
		Composition cEroded =  _composition.divideMasses(1.0f - f);
		// allow reseting of center position for flexibility in the
		// growth (e.g. to implement motility)
		_centerX = c.x;
		_centerY = c.y;
		_centerZ = c.z;
		// float m = _composition.grow(t, getCenter());
		// update radius
		_r = _composition.computeRadius();
		// subtract radius to determine radius change
		_radiusChangeInThisIteration += (_r - radiusBeforeErosion);
		return cEroded;
	}
/*	public BiomassSpecies.Composition erode(float f) {
		// criterium changed to 1.1 instead of 1 as value 1 can occur
		// (although chances of this happening are very low) due to precision
		// of float type
		if ((f < 0) | (f >= 1.1) | (Float.isNaN(f)))
			throw new ModelRuntimeException("illegal value f = " + f
					+ ": f must be 0<f<1");
		// in case f is 1 or (slightly) higher, put value 0.999999f
		f = (f >= 1 ? 0.999999f : f);
		return _composition.divideMasses(1.0f - f);
	}*/

	/**
	 * Get the timeconstrint for this particle. The growth rates used for
	 * particle growth are precomputed in this step
	 * 
	 * @return volumetric growth rate for this particle [um3/h]
	 */
	public float getMaximumTimeStep() {
		_composition.computeGrowthRates(getCenter());
		return _composition.getMaximumTimeConstraint();
	}

	/**
	 * Get the volumetric growth rate for this particle.
	 * 
	 * @return volumetric growth rate for this particle [um3/h]
	 */
	public float getVolumetricGrowthRate() {
		return _composition.getVolumetricGrowthRate();
	}

	/**
	 * Check if cell should be removed from system.
	 * 
	 * @return true if _radius is lower than minimum radius
	 */
	public boolean isDead() {
		return _r < _particleContainer.getMinimumRadius();
	}

	/**
	 * Check if particle is to be removed fot being out of bounds ( which is
	 * true if it is above the system limit for plannar geometry).
	 * 
	 * @return true if particle is out of bounds
	 */
	public boolean outOfBounds() {
		return _centerX > _model.getMaximumBiofilmHeight();
	}

	/**
	 * @return true if particle is marked for detachment by erosion
	 */
	public boolean willDetachByErosion() {
		return _willDetachByErosion;
	}

	/**
	 * @return true if particle is marked for detachment by sloughing
	 */
	public boolean willDetachBySloughing() {
		return _willDetachBySloughing;
	}

	/**
	 * Mark this particle for detachment by erosion
	 */
	public void setToDetachByErosion() {
		_willDetachByErosion = true;
	}

	/**
	 * Mark this particle for detachment by sloughing
	 */
	public void setToDetachBySloughing() {
		_willDetachBySloughing = true;
	}

	/**
	 * Check if particle is ready for division.
	 * 
	 * @return true if radius is higher than maximum radius
	 */
	public boolean willDivide() {
		return _r > _particleContainer.getMaximumRadius();
	}

	/**
	 * Perform particle division. Reduce mass of this particle to 1/2 and return
	 * an identical particle. Both particle are shifted slightly following a
	 * random chosen direction.
	 * 
	 * @return The "daughter" biomass particle
	 */
	public BiomassParticle divide() {
		// add 0.01 for slightly assimetrical division
		// must not be too high, otherwise cumulative effects can lead to
		// particles that are too big.
		// This is an effect not accounted for in the definition of maximum
		// iteration step
		float divisionFraction = 0.5f + _model.getRandom() * 0.2f;
		// random division direction
		_divisionTheta = _model.getRandom() * ExtraMath.PI2;
		_divisionPhi = _model.getRandom() * ExtraMath.PI;
		// save the volume for later;
		float volumeBeforeDivision = computeVolume();
		float radiusBeforeDivision = _r;
		try {
			// cloning mother cell to new cell
			BiomassParticle baby = (BiomassParticle) this.clone();
			// put baby in shoving matrix
			baby.setShovingMapPositionAndAddToShovingMap();
			
			// divide masses asymmetrically among mother and daughter (mass
			// balance)
			baby._composition = _composition.divideMasses(divisionFraction);
			// Update radius of both mother and daughter
			_r = _composition.computeRadius();
			baby._r = baby._composition.computeRadius();
			// because baby had no volume before being created
			// update
			_radiusChangeInThisIteration = _r - radiusBeforeDivision;
			baby._radiusChangeInThisIteration = baby._r;
			//
			// create a movement direction for baby
			float delx = _r
					* (float) (Math.cos(_divisionTheta) * Math
							.sin(_divisionPhi));
			float dely = _r
					* (float) (Math.sin(_divisionTheta) * Math
							.sin(_divisionPhi));
			float delz = _r * (float) (Math.cos(_divisionPhi));
			// move baby cell
			// TODO remove try catch block after testing
			//try {
				baby.move(delx, dely, delz);
				// u pdate position of mother cell
				move(-delx, -dely, -delz);
				// return reference to baby
			/*} catch (ArrayIndexOutOfBoundsException e) {
				System.out.println("divisionFraction = " + divisionFraction);
				System.out.println("volumeBeforeDivision = "
						+ volumeBeforeDivision);
				System.out.println("radiusBeforeDivision = "
						+ radiusBeforeDivision);
				System.out.println("_r = " + _r);
				System.out.println("_composition = " + _composition);
				System.out.println("baby._r = " + baby._r);
				System.out.println("baby._composition = " + baby._composition);

				throw new ModelRuntimeException(e.toString());
			}*/
			return baby;
		} catch (CloneNotSupportedException e) {
			// Do nothing, the exception will never be thrown
		}
		// just so that it works, this code is never reached
		return null;
	}

	/**
	 * Check if particle is ready for excreting Eps.
	 * 
	 * @return true if volumetric fraction of eps is higher than 0.5
	 */
	public boolean willExcreteEps() {
		float f = _composition.getEpsVolumetricFraction();
		return ((f > 0.7f) && (f < 1.0f));
	}

	/**
	 * Excrete eps as a particle of the same species but with only eps in its
	 * compoition. The Eps only particle will contain all the EPS from the
	 * original particle
	 * 
	 * @return The "eps only" biomass particle
	 */
	public BiomassParticle excreteEps() {
		// random division direction
		float epsPlacementTheta = _model.getRandom() * ExtraMath.PI2;
		float epsPlacementPhi = _model.getRandom() * ExtraMath.PI;
		// save the radius for later;
		float radiusBeforeDivision = _r;
		try {
			// cloning mother cell to new cell
			BiomassParticle excretedEps = (BiomassParticle) this.clone();
			// put new particle in shoving matrix
			excretedEps.setShovingMapPositionAndAddToShovingMap();
			// divide masses assymetrically among mother and daughter (mass
			// balance)
			excretedEps._composition = _composition.removeCapsule();
			// udpate radius of both mother and daughter
			_r = _composition.computeRadius();
			excretedEps._r = excretedEps._composition.computeRadius();
			// because baby had no volume before being created
			// update
			_radiusChangeInThisIteration = _r - radiusBeforeDivision;
			excretedEps._radiusChangeInThisIteration = excretedEps._r;
			//
			// create a movement direction for baby
			float delx = _r
					* (float) (Math.cos(epsPlacementTheta) * Math
							.sin(epsPlacementPhi));
			float dely = _r
					* (float) (Math.sin(epsPlacementTheta) * Math
							.sin(epsPlacementPhi));
			float delz = _r * (float) (Math.cos(epsPlacementPhi));
			// move baby cell
			excretedEps.move(delx, dely, delz);
			// update position of mother cell
			move(-delx, -dely, -delz);
			// return reference to baby
			return excretedEps;
		} catch (CloneNotSupportedException e) {
			// Do nothing, the exception will never be thrown
		}
		// just so that it works, flow never reaches this point
		return null;
	}

	/**
	 * Set the position of this cell according to a discrete map. Add Bacterium
	 * to the shoving grid at the grid element defined by __shovN, _shovM and
	 * _shovL.
	 */
	protected final void setShovingMapPositionAndAddToShovingMap() {
		_particleContainer.addToShovingGrid(this, _centerX, _centerY, _centerZ);
	}

	/**
	 * Remove an entry from a shoving grid.
	 */
	protected void removeFromShovingMap() {
		try {
			_particleContainer.shovingGrid[_shovN][_shovM][_shovL].remove(this);
		} catch (BacteriaNotInSetException e) {
			throw new ModelRuntimeException("Tried to remove a non existing"
					+ " element");
		}
	}

	/**
	 * @return true if the location in the shoving map directly above is
	 *         occupied
	 */
	protected boolean hasNeighborOnTop() {
		int n = _shovN;
		while (n++ < _particleContainer._n)
			if (!shovingGridEmpty(n, _shovM, _shovL))
				return true;
		return false;
	}

	/**
	 * @return the approximate distance to the nearest neighbor directly below
	 *         this one
	 */
	protected float distanceToNearestNeighborBelow() {
		int ngrids = 0;
		int n = _shovN - 1;
		// return immediately if particle is placed in the substrate
		if (n == -1)
			return _centerX - _r;
		// otherwise, check is particle is floating
		while (shovingGridEmpty(n--, _shovM, _shovL)) {
			ngrids++;
			// if n == 1 it means substrate was reached
			// move particle all the way to the bottom
			if (n == -1)
				return _centerX - _r;
		}
		return ngrids * _particleContainer.shovingGridSide;
	}

	/**
	 * @return the approximate distance to the nearest empty grid node of the
	 *         shoving grid
	 */
	protected float distanceToNearestEmptyPlaceAbove() {
		int ngrids = 0;
		int n = _shovN + 1;
		// otherwise, check is particle is floating
		while (!shovingGridEmpty(n++, _shovM, _shovL)) {
			ngrids++;
		}
		return ngrids * _particleContainer.shovingGridSide;
	}

	/**
	 * @param n
	 * @param m
	 * @param l
	 * @return true if the shoving map grid is empty at that location
	 */
	private boolean shovingGridEmpty(int n, int m, int l) {
		if (_particleContainer.shovingGrid[n][m][l] == null)
			return true;
		if (_particleContainer.shovingGrid[n][m][l].getNumberOfParticles() == 0)
			return true;
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
		// implement solid substratum for planar geometry
		_centerX = (x < _r ? _r : x);
		// Implement periodic boundaries (side walls)
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
	 * Move center of bacteria by vector addition with m. Uses setCenter so
	 * implementation of system boundary conditions is implicit.
	 * 
	 * @param x
	 *            [microns]
	 * @param y
	 *            [microns]
	 * @param z
	 *            [microns]
	 */
	public void move(float x, float y, float z) {
		removeFromShovingMap();
		setCenter(_centerX - x / _biomassSpecies.getShovingHierarchy(),
				_centerY - y / _biomassSpecies.getShovingHierarchy(), _centerZ
						- z / _biomassSpecies.getShovingHierarchy());
	}

	/**
	 * Add contribution of a neighbor to the moving direction if there is
	 * overlapping.
	 * 
	 * @param neib
	 *            Neighboring cell that's potentially overlapping
	 * @param m
	 *            move direction to be updated in case of overlapping
	 * @return true if there is overlapping, false otherwise
	 */
	protected boolean addOverlappingMoveDirection(BiomassParticle neib,
			ContinuousCoordinate m) {
		// check if the neighbor tested is actually this cell
		if (neib == this)
			return false;
		// check shoving hierarchy and shove only if neighbor has higher
		/*if (neib._biomassSpecies.getShovingHierarchy() < _biomassSpecies
				.getShovingHierarchy())
			return false;*/
		// if not, then calculate if there is overlap
		float neibX = neib._centerX;
		float neibY = neib._centerY;
		float neibZ = neib._centerZ;
		float r = neib._r;
		float dx;
		float dy;
		float dz;
		float dxq = _particleContainer.shovingGridSide;
		// no periodicity for vertical direction
		dx = neibX - _centerX;
		// take periodicity into account in y-direction
		float Lyy = _model.systemSize.y;
		if (neibY < dxq && _centerY > Lyy - dxq)
			// neighbor on the left, this on right
			dy = Lyy + neibY - _centerY;
		else if (_centerY < dxq && neibY > Lyy - dxq)
			// this on the left, neighbor on right
			dy = neibY - Lyy - _centerY;
		else
			// the regular case
			dy = neibY - _centerY;
		// take periodicity into account in z-direction
		float Lzz = _model.systemSize.z;
		if (neibZ < dxq && _centerZ > Lzz - dxq)
			dz = Lzz + neibZ - _centerZ;
		else if (_centerZ < dxq && neibZ > Lzz - dxq)
			dz = neibZ - Lzz - _centerZ;
		else
			dz = neibZ - _centerZ;
		// length of difference vector between
		// the two cells
		float d = (float) Math.sqrt(dx * dx + dy * dy + dz * dz);
		// if d is zero (overlapping cells) put a
		if (d == 0)
			d = 1e-3f * _particleContainer.getMaximumRadius();
		// normalized overlapping distance overlapping length
		float nod = (_particleContainer.shovingParameter * (_r + r) - d) / d;
		if (nod > 0.1) {
			m.x += dx * nod;
			m.y += dy * nod;
			m.z += dz * nod;
			return true;
		}
		return false;
	}

	/**
	 * @param neib
	 * @return the distance to particle neib taking periodicity into account
	 */
	public float distanceTo(BiomassParticle neib) {
		float neibX = neib._centerX;
		float neibY = neib._centerY;
		float neibZ = neib._centerZ;
		float dx;
		float dy;
		float dz;
		// no periodicity for vertical direction
		dx = neibX - _centerX;
		// take periodicity into account in y-direction
		float Lyy = _model.systemSize.y;
		dy = ExtraMath.abs(neibY - _centerY);
		dy = ExtraMath.min(dy, Lyy - dy);
		// take periodicity into account in z-direction
		float Lzz = _model.systemSize.z;
		dz = ExtraMath.abs(neibZ - _centerZ);
		dz = ExtraMath.min(dz, Lzz - dz);
		// length of difference vector between
		// the two cells
		return (float) Math.sqrt(dx * dx + dy * dy + dz * dz);
	}

	/**
	 * Returns the _shovingMapPosition. NOTE: this practice is not very good, as
	 * the _shovingMapPosition may be changed
	 * 
	 * @return _shovN
	 */
	protected int getShovingMapN() {
		return _shovN;
	}

	/**
	 * @return _shovM
	 */
	protected int getShovingMapM() {
		return _shovM;
	}

	/**
	 * @return _shovL
	 */
	protected int getShovingMapL() {
		return _shovL;
	}

	/**
	 * Set the value of the shoving map coordinates
	 * 
	 * @param n
	 * @param m
	 * @param l
	 */
	protected void setShovingMapPosition(int n, int m, int l) {
		_shovN = n;
		_shovL = l;
		_shovM = m;
	}

	/**
	 * @return the cell center coordinates
	 */
	public ContinuousCoordinate getCenter() {
		return new ContinuousCoordinate(_centerX, _centerY, _centerZ);
	}

	public float getCenterX() {
		return _centerX;
	}

	public float getCenterY() {
		return _centerY;
	}

	public float getCenterZ() {
		return _centerZ;
	}

	/**
	 * @return the cell radius
	 */
	public float getRadius() {
		return _r;
	}

	/**
	 * Computes the volume of the biomass particle from the radius and the
	 * dimensionality of the system: if 2D the volume is compute assuming the
	 * particle is a cylinder, if 3D volume is computed considering spherical
	 * particle
	 * 
	 * @return the volume of the particle [um3]
	 */
	public float computeVolume() {
		return computeVolume(_r);
	}

	/**
	 * Computes the volume of a biomass particle with radius r from the
	 * dimensionality of the system: if 2D the volume is compute assuming the
	 * particle is a cylinder, if 3D volume is computed considering spherical
	 * particle
	 * 
	 * @param r
	 *            radius of particle
	 * @return the volume of the particle [um3]
	 */
	public static float computeVolume(float r) {
		float v = 0;
		// Volume of particle is computed either from sphere (3D) or
		// cylinder (2D)
		v = (Model.model().getDimensionality() == 3 ? ExtraMath
				.volumeOfASphere(r) : ExtraMath.volumeOfACylinder(r, Model
				.model().get2DSystem3rdDimension()));
		return v;
	}

	/**
	 * Determine the radius of spherical particle core (all masses except eps)
	 * 
	 * @return radius of spherical particle core
	 */
	public float getCoreRadius() {
		return _composition.getCoreRadius();
	}

	/**
	 * Returns the color for this bacterium, based on its growth rate
	 * 
	 * @return the color species of the bacteria (or _color if overidden)
	 */
	public Color getColorCore() {
		if (_overideColor)
			return _color;
		return _composition.getColorCore(getCenter());
	}

	/**
	 * @return true if particle has EPS capsule
	 */
	public boolean hasCapsule() {
		return _biomassSpecies.hasEpsCapsule();
	}

	/**
	 * @return the color to draw capsule
	 */
	public Color getColorCapsule() {
		return _composition.getColorCapsule();
	}

	/**
	 * Set the active mass of this bacterium (and update valu of radius)
	 * 
	 * @param s
	 *            fixed species to set value
	 * @param m
	 *            value to set mass to
	 */
	public void setMass(ParticulateSpecies s, float m) {
		_composition.setMass(s, m);
		// update the radius
		_r = _composition.computeRadius();
	}

	/**
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		Color colorCore = getColorCore();
		Color colorCapsule = getColorCapsule();
		return _centerX + "\t" + _centerY + "\t" + _centerZ + "\t"
				+ getCoreRadius() + "\t" + colorCore.getRed() + "\t"
				+ colorCore.getGreen() + "\t" + colorCore.getBlue() + "\t" + _r
				+ "\t" + colorCapsule.getRed() + "\t" + colorCapsule.getGreen()
				+ "\t" + colorCapsule.getBlue() + "\n";
	}

	/**
	 * Add contribution of this particle to the discrete data matrices
	 */
	protected void addContributionToBacteriaSpeciesDiscreteData() {
		_composition.addContributionToDiscreteData(getCenter());
	}

	/**
	 * Add contribution of a neighboring to the moving direction which is a
	 * consequence of pressures caused by particle growth/shricking
	 * 
	 * @param neib
	 *            neighbor that is acting uppon this cells
	 * @param m
	 *            move direction to be updated in case of overlapping
	 * @return true if there is overlapping, flase otherwise
	 */
	protected void addInfluenceToPressureMoveDirection(BiomassParticle neib) {
		float neibX = neib._centerX;
		float neibY = neib._centerY;
		float neibZ = neib._centerZ;
		float r = neib._r;
		float dx;
		float dy;
		float dz;
		float dxq = _particleContainer.shovingGridSide;
		// no periodicity for vertical direction
		dx = neibX - _centerX;
		// take periodicity into account in y-direction
		float Lyy = _model.systemSize.y;
		if (neibY < dxq && _centerY > Lyy - dxq)
			dy = Lyy + neibY - _centerY;
		else if (_centerY < dxq && neibY > Lyy - dxq)
			dy = neibY - Lyy - _centerY;
		else
			dy = neibY - _centerY;
		// take periodicity into account in z-direction
		float Lzz = _model.systemSize.z;
		if (neibZ < dxq && _centerZ > Lzz - dxq)
			dz = Lzz + neibZ - _centerZ;
		else if (_centerZ < dxq && neibZ > Lzz - dxq)
			dz = neibZ - Lzz - _centerZ;
		else
			dz = neibZ - _centerZ;
		// length of difference vector between
		// the two cells
		float d = (float) Math.sqrt(dx * dx + dy * dy + dz * dz);
		// if d is zero (overlapping cells) put a
		if (d == 0)
			d = 1e-3f * _particleContainer.getMaximumRadius();
		// only bring cells together if they are not overlapping
		float multiplier = neib._radiusChangeInThisIteration / d;
		if (multiplier < 0) {
			_moveByPressureVector.x += dx * multiplier;
			_moveByPressureVector.y += dy * multiplier;
			_moveByPressureVector.z += dz * multiplier;
		}
	}

	/**
	 * Reset the move by pressure vector
	 */
	public void resetMoveByPressureVector() {
		_moveByPressureVector = new ContinuousCoordinate();
	}

	/**
	 * Move the particle center according to the _moveByPressureVector
	 */
	public void moveByPressure() {
		move(_moveByPressureVector.x, _moveByPressureVector.y,
				_moveByPressureVector.z);
	}

	/**
	 * @return the total mass of a particle
	 */
	public float getTotalMass() {
		return _composition.getTotalMass();
	}

	/**
	 * @return the mass of Eps in a particle
	 */
	public float getTotalEps() {
		return _composition.getEpsMass();
	}

	/**
	 * @return true if particle is solely composed of EPS
	 */
	public boolean isEpsOnly() {
		return _composition.getTotalMass() == _composition.getEpsMass();
	}

	/**
	 * @return the toal mass after subtracting the EPS mass
	 */
	public float getActiveMass() {
		return _composition.getTotalMass() - _composition.getEpsMass();
	}

	/**
	 * @param s
	 *            species to get mass of
	 * @return the mass of species in a particle
	 */
	public float getMassOfSpecies(ParticulateSpecies s) {
		return _composition.getSpeciesMass(s);
	}

	/**
	 * @param c
	 *            The _color to set.
	 */
	public void overideColor(Color c) {
		_overideColor = true;
		_color = c;
	}

	/**
	 * @return Returns the _composition.
	 */
	public BiomassSpecies.Composition getComposition() {
		return _composition;
	}
}