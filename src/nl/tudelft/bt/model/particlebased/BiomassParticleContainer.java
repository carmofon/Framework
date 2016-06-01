package nl.tudelft.bt.model.particlebased;

import java.io.File;
import java.io.IOException;
import java.util.*;

import nl.tudelft.bt.model.BiomassContainer;
import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.ContinuousCoordinate;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.detachment.*;
import nl.tudelft.bt.model.detachment.cvf.ConnectedToBottomCvf;
import nl.tudelft.bt.model.detachment.cvf.ConnectedToTopCvf;
import nl.tudelft.bt.model.detachment.cvf.ConnectedVolumeFilter;
import nl.tudelft.bt.model.detachment.levelset.DensityLevelSet;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Height2MassDetachment;
import nl.tudelft.bt.model.detachment.levelset.functions.HeightVolumetricDetachment;
import nl.tudelft.bt.model.detachment.levelset.functions.MassDetachment;
import nl.tudelft.bt.model.detachment.levelset.functions.VolumetricDetachment;
import nl.tudelft.bt.model.exceptions.ModelException;
import nl.tudelft.bt.model.exceptions.ModelRuntimeException;
import nl.tudelft.bt.model.multigrid.MultigridUtils;
import nl.tudelft.bt.model.multigrid.ParticulateSpecies;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;
import nl.tudelft.bt.model.util.ContinuousCoordinateOps;
import nl.tudelft.bt.model.util.ExtraMath;
import nl.tudelft.bt.model.work.will.ChemotacticBiomassParticle;

/**
 * Container of hard spheres for modeling of biomass spreading through shoving
 * (individual based modeling). Stores all instances of bacteria in a map,
 * ordered by neighborhoods.
 * 
 * @author Jo?o Xavier (j.xavier@tnw.tudelft.nl)
 */
public class BiomassParticleContainer extends BiomassContainer {
	protected float shovingParameter;

	protected float shovingGridSide;

	protected DetachmentHandler _detachmentHandler;

	protected float _referenceSystemSide;

	private float _shovingFraction; // fraction of particles to consider

	// maximum and minimum radius of biomass particles
	private float _maximumRadius;

	private float _minimumRadius;

	private float _minimumMass;

	// relaxation process complete
	// sizes of shoving grid
	protected int _l;

	protected int _m;

	protected int _n;

	// The structures which hold the references to the particles
	protected BiomassParticleSet[][][] shovingGrid;

	protected List<BiomassParticle> particleList;

	// For abrasion operations
	private DensityLevelSet _densityLevelSet;

	//
	private float _biomassProducedInThisIteration;

	// refrence to model
	protected Model _model = Model.model();

	// an array for determining the potential shovers
	private BiomassParticleSet _ps;

	protected ConnectedVolumeFilter _cvfForBiomass;

	protected ConnectedVolumeFilter _cvfForLiquid;

	protected boolean[][][] _map;

	private boolean[][][] _totalLiquid;

	private boolean[][][] _carrierElements;

	private boolean _sloughinfIsOn = true;

	/**
	 * Creates a new BacteriaConainer with a shoving grid which is optimized for
	 * shoving. The grid is defined from a "reference system size".
	 * 
	 * @param maximumRadius
	 * @param minimumRadius
	 * @param kshov
	 * @param f
	 * @param rdetach
	 */
	public BiomassParticleContainer(float maximumRadius, float minimumRadius,
			float kshov, float f) {
		// create a bacteria set with enough space for 3*3*3 = 27
		// neighboring shoving grid elements
		_ps = new BiomassParticleSet(27);
		//
		_maximumRadius = maximumRadius;
		_minimumRadius = minimumRadius;
		// Value of minum mass will be set at first read of variable value
		// to ensure that species have been defined
		_minimumMass = 0;
		//
		shovingParameter = kshov;
		_shovingFraction = f;
		// determine the optimum value for grid size
		double o = _model.referenceSystemSide * 1.1 / kshov
				/ (2 * _maximumRadius);
		// and its floored integer
		int ofloor = (int) Math.round(o);
		ofloor = (ofloor < o ? ofloor : ofloor - 1);
		// shoving grid side
		_referenceSystemSide = _model.referenceSystemSide;
		shovingGridSide = _model.referenceSystemSide / ofloor;
		// shoving grid sizes
		// n includes one more grid level for cells that are about
		// to detach
		_n = Math.round(_model.systemSize.x / shovingGridSide);
		_m = Math.round(_model.systemSize.y / shovingGridSide);
		_l = Math.round(_model.systemSize.z / shovingGridSide);
		// in case of 2D (z = 0) put _n to 1)
		_l = (_l == 0 ? 1 : _l);
		// initialize shoving grid and bacteria list
		shovingGrid = new BiomassParticleSet[_n + 10][_m][_l];
		particleList = new ArrayList();
		//
		_biomassProducedInThisIteration = 0;
		// create connected volume filtrator
		createConnectedVolumeFilters();
		// initialize auxiliary matrices for detachment computation
		_map = new boolean[_n][_m][_l];
		_totalLiquid = new boolean[_n + 10][_m][_l];
	}

	/**
	 * Initialize the connected volume filters
	 */
	protected void createConnectedVolumeFilters() {
		_cvfForBiomass = new ConnectedToBottomCvf(_n, _m, _l);
		_cvfForLiquid = new ConnectedToTopCvf(_n + 10, _m, _l);
	}

	/**
	 * Grow each biomass particle according to a time step
	 * 
	 * @param t
	 *            time step
	 */
	public void grow(float t) {
		// update the values of maximum specific growth rate for ech fixed
		// species
		_biomassProducedInThisIteration = 0;
		_model.updateValuesOfCurrentMaximumSpecificGrowthRates();
		for (Iterator it = particleList.iterator(); it.hasNext();) {
			BiomassParticle b = (BiomassParticle) it.next();
			_biomassProducedInThisIteration += b.grow(t);
		}
	}

	/**
	 * Iterate through bacteria and perform division when necessary.
	 */
	public void divideAndExcreteEPS() {
		// total number of bacteria at this time
		int nbac = particleList.size();
		// index of new bacteria in bacteria List
		int inew = nbac;
		// reset the divided particle counter for each biomass species
		for (Iterator iter = _model.getBiomassSpecies().iterator(); iter
				.hasNext();) {
			BiomassSpecies s = (BiomassSpecies) iter.next();
			// get number of particles to attach
			s.resetCreatedParticleCounter();
		}
		// iterate through bacteria checking which to divide
		for (int i = 0; i < nbac; i++) {
			// get the current bacterium
			BiomassParticle b = (BiomassParticle) particleList.get(i);
			// excrete EPS if its the case
			if (b.willExcreteEps()) {
				BiomassParticle epsOnlyParticle = b.excreteEps();
				// add eps only particle to the list
				particleList.add(epsOnlyParticle);
				// check if eps only particle is also larger than maximum size
				// and if so divide it and add also this new one
				if (epsOnlyParticle.willDivide()) {
					BiomassParticle epsOnlyParticle2 = epsOnlyParticle.divide();
					// add baby eps particle to the list
					particleList.add(epsOnlyParticle2);
				}
			}
			// divide if its the case
			if (b.willDivide()) {
				BiomassParticle baby = b.divide();
				// add baby particle to the list
				particleList.add(baby);
				// Increment the divided particle counter
				b.getBiomassSpecies().incrementCreatedParticleCounter();
			}
		}
	}

	/**
	 * For each biomass species in the system perform attachment based on the
	 * number of particles created at this iteration
	 */
	public void attach() {
		// for each biomass species perform
		for (Iterator iter = _model.getBiomassSpecies().iterator(); iter
				.hasNext();) {
			BiomassSpecies s = (BiomassSpecies) iter.next();
			// get number of particles to attach
			int n = (int) (0.1 * s
					.getNumberOfParticlesCreatedInPresentIteration());
			// Create a n new particles of this species and attach to the
			// biofilm
			for (int i = 0; i < n; i++) {
				BiomassParticle b = s.createBiomassParticle();
				ContinuousCoordinate c = determineAttachmentPositionFromRandomWalk();
				b.setCenter(c.x, c.y, c.z);
				// add to the bacteria list
				particleList.add(b);
			}
		}
	}

	/**
	 * Relax completly by shoving
	 */
	public void spread() {
		// spread (or shrink) by internet pressure
		// spreadByShoving();
		// TODO: pulling is no
		performSpreadingByPulling();
		spreadByShoving();
	}

	public void spreadByShoving() {
		// Undo particle overlaps
		// int shovingIteration = 0;
		int shovlimit = (int) (particleList.size() * _shovingFraction);

		// Shoving until a completly stable system is achieved is:
		// while ( niter < _maximumShovingIteration && performShovingStep() );
		// here we use a fraction of shoved cells for increased efficiency
		while (performSpreadingStep() > shovlimit) {
			// do nothing
		}
	}

	/**
	 * Perform a single shoving iteration
	 * 
	 * @return number of shoved particles in this iteration
	 */
	public int performSpreadingStep() {
		// create checker for shoving overlap in any of the
		// cells
		int shoved = 0;
		// create a new movement vector
		ContinuousCoordinate m = new ContinuousCoordinate();
		// iterate through bacteria list
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle b = (BiomassParticle) iter.next();
			if (true) {
				// create checker for occurence of shoving in this
				// bacteria
				boolean thisShoved = false;
				// set to zeros the movement vector
				m.reset();
				// get an array of potential shovers from the
				// shoving grid
				BiomassParticleSet potentialShovers = getPotentialShovers(b);
				int numberOfOverlappingNeighbors = 0;
				while (potentialShovers.hasNext()) {
					BiomassParticle neib = potentialShovers.next();
					if (b.addOverlappingMoveDirection(neib, m)) {
						thisShoved = true;
						numberOfOverlappingNeighbors++;
					}
				}
				// update cell position
				if (thisShoved) {
					shoved++;
					// move the cell
					b.move(m.x, m.y, m.z);
				}
			} // end if(b.willDetach())
		}
		return shoved;
	}

	/**
	 * Get the neighbors from shoving map that may be overlapping
	 * 
	 * @param c
	 * @return set of barticles in neighboring grids of the shoving map
	 */
	public BiomassParticleSet getPotentialShovers(BiomassParticle b) {
		// Pre-allocate
		int m, n, l;
		//
		int i = b.getShovingMapN();
		int j = b.getShovingMapM();
		int k = b.getShovingMapL();
		// reset the array that holds the potential shoving neighbors
		_ps.reset();

		// check 8-neighbors for Bacteria
		// check neighbors in vertical direction
		// (no periodicity)
		for (n = (i == 0 ? 0 : i - 1); n <= (i == _n - 1 ? _n - 1 : i + 1); n++) {
			// Horizontal directions (with periodicity):
			for (m = j - 1; m <= j + 1; m++) {
				// wrapping of edges - wrapped m
				int mwrap = (m < 0 ? _m - 1 : (m == _m ? 0 : m));
				if (_l > 1) {
					for (l = k - 1; l <= k + 1; l++) {
						// 3D case
						int lwrap = (l < 0 ? _l - 1 : (l == _l ? 0 : l));
						_ps.addAll(shovingGrid[n][mwrap][lwrap]);
					}
				} else {
					// 2D case
					_ps.addAll(shovingGrid[n][mwrap][0]);
				}
			}
		}
		return _ps;
	}

	/**
	 * Remove all cells to be detached by erosion, sloughing and that reached
	 * out of bounds
	 */
	public void removeDetachedBiomass() throws ModelException {
		// execute detachMarkedAndOutOfBoundsBiomass() 3 times!
		// 1st - to remove out of bonds biomass
		detachMarkedAndOutOfBoundsBiomass();
		// if a detachment function is not defined or the detachment is set off
		// do not go through with detachment routines
		if ((_detachmentHandler != null)
				& (!_detachmentHandler.detachmentIsOff())) {
			// 2nd - to remove eroded biomass
			markParticlesForErosion();
			detachMarkedAndOutOfBoundsBiomass();
			// 3rd - to remove eroded biomass
			if (_sloughinfIsOn) {
				mapBiofilmElements();
				markParticlesForSloughing();
				detachMarkedAndOutOfBoundsBiomass();
			}
			// Erode border particles to complete detachment
			// update the liquid elements matrix after having removed particles
			// by other mechanisms
			_densityLevelSet.setBulkLiquidElements(getLiquidElements());
			_densityLevelSet.erodeBorder();
		} else {
			// skip erosion but still perform detachment
			mapBiofilmElements();
			detachMarkedAndOutOfBoundsBiomass();
		}
	}

	/**
	 * Remove all cells marked for detachment and include them in the eroded
	 * category
	 */
	protected void detachMarkedAndOutOfBoundsBiomass() {
		// Don't iterate using iterator because list is changed
		// each time an element is removed
		for (int i = 0; i < particleList.size(); i++) {
			// get the current bacterium
			BiomassParticle b = (BiomassParticle) particleList.get(i);
			if (b.outOfBounds() | b.willDetachByErosion()
					| b.willDetachBySloughing()) {
				remove(b);
				// account for detached biomass
				Model.model().detachedBiomassContainer()
						.addToDetachedBiomass(b);
				// decrement index to repeat this value
				i--;
			}
		}
	}

	/**
	 * Detach all particles that should be removed in this iteration and compute
	 * the density level set which can be used later for eroding remaining
	 * border particles
	 */
	public void markParticlesForErosion() {
		// also initializes _carrierElements
		boolean[][][] liq = getLiquidElements();
		_densityLevelSet = new DensityLevelSet(shovingGrid, liq,
				_carrierElements, shovingGridSide,
				(DetachmentSpeedFunction) _detachmentHandler,
				_model.getTimeStep());
	}

	/**
	 * Method used for testing phase
	 * 
	 * @throws IOException
	 */
	public void writeLevelSetUsingDiverseDetachmentFunctions()
			throws IOException {
		DetachmentSpeedFunction[] functions = {
				new VolumetricDetachment(10f / 70.0f), new MassDetachment(10f),
				new HeightVolumetricDetachment(10f / 70.0f / 50f),
				new Height2MassDetachment(1.0e-3f) };
		for (int i = 0; i < functions.length; i++) {
			boolean[][][] liq = getLiquidElements();
			_densityLevelSet = new DensityLevelSet(shovingGrid, liq,
					_carrierElements, shovingGridSide, functions[i],
					_model.getTimeStep());
			writeDetachmentLevelSet("D:\\temp\\"
					+ functions[i].getClass().getName() + ".txt");
		}
	}

	/**
	 * @return a list with all elements
	 */
	private boolean[][][] getLiquidElements() {
		// map total liquid
		for (int i = 0; i < _n + 10; i++)
			for (int j = 0; j < _m; j++)
				for (int k = 0; k < _l; k++) {
					if (getCarrierElement(i, j, k))
						// set to false if it belongs to the carrier
						_totalLiquid[i][j][k] = false;
					else
						// set to true if entry has no particles, false
						// otherwise
						_totalLiquid[i][j][k] = ((shovingGrid[i][j][k] == null) ? true
								: (shovingGrid[i][j][k].getNumberOfParticles() == 0));
				}
		boolean[][][] cvfContainer = _cvfForLiquid.computeCvf(_totalLiquid);
		return cvfContainer;
	}

	/**
	 * Update the map array with all elements in the biofilm
	 */
	protected void mapBiofilmElements() {
		// map all biofilm elements
		for (int i = 0; i < _n; i++)
			for (int j = 0; j < _m; j++)
				for (int k = 0; k < _l; k++)
					// set map entry to true if there is particle in the
					// position
					// false, otherwise
					_map[i][j][k] = (((shovingGrid[i][j][k] != null) ? (shovingGrid[i][j][k]
							.getNumberOfParticles() > 0) : false));
	}

	/**
	 * Determine which cells will detach via sloughing and mark them
	 */
	protected void markParticlesForSloughing() {
		// perform connected volume filtration (connected to bottom
		boolean[][][] cvf = _cvfForBiomass.computeCvf(_map);
		// mark as detachable all particles in nonvalid map positions
		for (int i = 0; i < _n; i++) {
			for (int j = 0; j < _m; j++) {
				for (int k = 0; k < _l; k++) {
					if (_map[i][j][k] & (!cvf[i][j][k])) {
						// for all particles in a map entry not connected to
						// the carrier
						shovingGrid[i][j][k].resetIndex();
						while (shovingGrid[i][j][k].hasNext())
							shovingGrid[i][j][k].next()
									.setToDetachBySloughing();
					}
				}
			}
		}
	}

	/**
	 * Iterate through all bacteria in the system removing "dead" cells.
	 */
	public void removeDeadBiomass() {
		// Get the current number of bacteria in the system, 'ncl'
		for (int i = 0; i < particleList.size(); i++) {
			// get the current bacterium
			BiomassParticle b = (BiomassParticle) particleList.get(i);
			if (b.isDead()) {
				remove(b);
				// decrement index to repeat this value
				i--;
			}
		}
	}

	/**
	 * Remove bacteria from the bacteria list and from the shoving grid
	 * 
	 * @param b
	 *            bacterium to remove
	 */
	public void remove(BiomassParticle b) {
		particleList.remove(b);
		b.removeFromShovingMap();
	}

	/**
	 * Snap a continuous coordinate to shoving grid indexes
	 * 
	 * @param x
	 *            spatial coordinate (continuous)
	 * @return shoving grid index (discrete)
	 */
	protected void addToShovingGrid(BiomassParticle b, float x, float y, float z) {
		// snap continuous coordinates
		int i = snapToGridI(x);
		int j = snapToGridJ(y);
		int k = snapToGridK(z);
		// initialize the particle set in this grid element if it is null
		if (shovingGrid[i][j][k] == null)
			shovingGrid[i][j][k] = new BiomassParticleSet();
		// add the particle
		shovingGrid[i][j][k].add(b);
		// and retrun the grid coordinates to update particle attributes
		b.setShovingMapPosition(i, j, k);
	}

	/**
	 * Convert a continuous coordinate into shoving grid coordinates
	 * 
	 * @param x
	 * @return shoving grid coordinates
	 */
	private int snapToGridI(float x) {
		int i = (int) (x / shovingGridSide);
		return i < _n ? i : 0;
	}

	/**
	 * Convert a continuous coordinate into shoving grid coordinates
	 * 
	 * @param y
	 * @return shoving grid coordinates
	 */
	private int snapToGridJ(float y) {
		int i = (int) (y / shovingGridSide);
		return i < _m ? i : 0;
	}

	/**
	 * Convert a continuous coordinate into shoving grid coordinates
	 * 
	 * @param z
	 * @return shoving grid coordinates
	 */
	private int snapToGridK(float z) {
		int i = (int) (z / shovingGridSide);
		return i < _l ? i : 0;
	}

	/**
	 * @param i
	 * @param j
	 * @param k
	 * @return the center of a grid element in continuous space
	 */
	public ContinuousCoordinate getGridElementCenter(int i, int j, int k) {
		return new ContinuousCoordinate(((float) i + 0.5f) * shovingGridSide,
				((float) j + 0.5f) * shovingGridSide, ((float) k + 0.5f)
						* shovingGridSide);
	}

	/**
	 * Return the number of cells of a species
	 * 
	 * @param s
	 *            species
	 * @return number of particles of this species presently in the container
	 */
	public int getNumberOfParticles(BiomassSpecies s) {
		int n = 0;
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle element = (BiomassParticle) iter.next();
			if (element.isOfSpecies(s))
				n++;
		}
		return n;
	}

	/**
	 * Get the biofilm height
	 * 
	 * @return biofilm height
	 */
	public float getTopBiofilmHeight() {
		// iterate through array of bacteria and increase
		// concentration value in corresponding BacteriaSpecies
		float h = 0;
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			// get bacterium height
			float x = ((BiomassParticle) iter.next()).getCenter().x;
			h = (x > h ? x : h);
		}
		return h;
	}

	/**
	 * Add contribution of biomass particles to the biomass discrete data
	 * matrices
	 */
	public void addContributionsToSpeciesDiscreteData() {
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle b = (BiomassParticle) iter.next();
			b.addContributionToBacteriaSpeciesDiscreteData();
		}
	}

	/**
	 * For a monospecies system, add a bacterium at a defined position with a
	 * defined mass
	 * 
	 * @param x
	 *            height from solid surface [um]
	 * @param y
	 *            horizontal coordinate [um]
	 * @param z
	 *            horizontal coordinate [um]
	 * @param m
	 *            mass
	 * @return reference to particle placed
	 */
	public BiomassParticle placeBiomassAt(BiomassSpecies s, float x, float y,
			float z) {
		BiomassParticle b = s.createBiomassParticle();
		b.setCenter(x, y, z);
		// add to the particle list
		particleList.add(b);
		return b;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.photobiofilms.model.BiomassContainer#getBiomassAsBacteriaCollection()
	 */
	public Collection getBiomassAsBiomassParticleCollection() {
		return Collections.unmodifiableCollection(particleList);
	}

	/**
	 * computes a time constrint that prevents the biofilm to increase in size
	 * relative to a function of the maximum particle radius
	 * 
	 * @return the biofilm matrix quantity contraint
	 */
	public TimeStepConstraint getGlobalSizeConstraint() {
		float biofilmSizeRate = 0; // area if 2D volume if 3D
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle b = (BiomassParticle) iter.next();
			biofilmSizeRate += Math.abs(b.getVolumetricGrowthRate());
		}
		float constraint = _l * _m
				* BiomassParticle.computeVolume(_maximumRadius)
				/ biofilmSizeRate * 3;
		// throw exception if constraint is NaN
		if (Float.isNaN(constraint))
			throw new ModelRuntimeException(
					"getGlobalSizeConstraint returned a NaN value");
		// otherwise, return the value
		return new TimeStepConstraint(constraint, "GlobalBiofilmSize");
	}

	/**
	 * The time constraint at the particle scale
	 * 
	 * @return the time step constraint from particle growth
	 */
	public TimeStepConstraint getParticleTimeConstraint() {
		float t = Float.POSITIVE_INFINITY;
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle b = (BiomassParticle) iter.next();

			// keep the lowest value
			t = Math.min(b.getMaximumTimeStep(), t);
			if (Float.isNaN(t)) {
				throw new ModelRuntimeException(
						"getParticleTimeConstraint produced a NaN value");
			}
		}
		TimeStepConstraint c = new TimeStepConstraint(t, "ParticleGrowth");
		return c;
	}

	/**
	 * Return the maximum radius that bacteria can assume in the system. This is
	 * used for purposes of defining size of elements in the shoving grid and
	 * the particle division criterium.
	 * 
	 * @return maximum radius of Bacteria in the system
	 */
	public float getMaximumRadius() {
		return _maximumRadius;
	}

	/**
	 * Return the minimum radius that bacteria can assume in the system. This is
	 * used for purposes of determining if a particle will be removed.
	 * 
	 * @return minimum radius of Bacteria in the system
	 */
	protected float getMinimumRadius() {
		return _minimumRadius;
	}

	/**
	 * Perform the spreading of particles by pressure created by neighboring
	 * particles growing/shrinking
	 */
	public void performSpreadingByPulling() {
		// iterate through bacteria list and compute movement vectors
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle b = (BiomassParticle) iter.next();
			b.resetMoveByPressureVector();
			// get an array of potential shovers from the
			// shoving grid
			BiomassParticleSet neighborsInShovingMap = getPotentialShovers(b);
			while (neighborsInShovingMap.hasNext()) {
				BiomassParticle neib = neighborsInShovingMap.next();
				if (neib != b)
					b.addInfluenceToPressureMoveDirection(neib);
			}
		}
		// iterate through bacteria list again and perform the movements
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle b = (BiomassParticle) iter.next();
			b.moveByPressure();
		}
	}

	/**
	 * @param function
	 *            The _detachmentFunction to set.
	 */
	public void setDetachmentFunction(DetachmentHandler function) {
		_detachmentHandler = function;
	}

	/**
	 * @return Returns the _detachmentFunction.
	 */
	public DetachmentHandler getDetachmentFunction() {
		return _detachmentHandler;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.BiomassContainer#getTotalProducedBiomass()
	 */
	public float getTotalProducedBiomass() {
		return _biomassProducedInThisIteration;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.photobiofilms.model.BiomassContainer#getCurrentAverageBiofilmHeight()
	 */
	public float getTotalBiomass() {
		float m = 0;
		for (Iterator iter = particleList.iterator(); iter.hasNext();) {
			BiomassParticle b = (BiomassParticle) iter.next();
			m += b.getTotalMass();
		}
		return m;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.photobiofilms.model.BiomassContainer#writeDetachmentLevelSet(java
	 * .lang.String)
	 */
	public void writeDetachmentLevelSet(String f) throws IOException {
		if (_detachmentHandler != null) {
			_densityLevelSet.writeToFile(f);
		} else {
			System.out.println("Warning! "
					+ "trying to write level set without detachment function");
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.photobiofilms.model.BiomassContainer#writeBiofilmFront(java.lang.
	 * String)
	 */
	public void writeBiofilmFront(String s) throws IOException {
		// get the liquid boolean matrix
		boolean[][][] le = getLiquidElements();
		// convert to biofilm boolean matrix (connected volume
		int n = le.length;
		int m = le[0].length;
		int l = le[0][0].length;
		float[][][] be = new float[le.length][le[0].length][le[0][0].length];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				for (int k = 0; k < l; k++) {
					//
					be[i][j][k] = (le[i][j][k] ? 0
							: (getCarrierElement(i, j, k) ? 0 : 1.0f));
				}
		// create the file and write matrix to it
		File f = new File(s + ".txt");
		java.io.FileWriter fr = new java.io.FileWriter(f);
		fr.write(MultigridUtils.matrixToString(be));
		fr.close();
	}

	/**
	 * method for debugging
	 */
	private void writeMatrices() throws IOException {
		// get the liquid boolean matrix
		boolean[][][] le = getLiquidElements();
		// convert to biofilm boolean matrix (connected volume
		int n = le.length;
		int m = le[0].length;
		int l = le[0][0].length;
		float[][][] be = new float[n][m][l];
		float[][][] carrier = new float[n][m][l];
		float[][][] liquid = new float[n][m][l];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				for (int k = 0; k < l; k++) {
					//
					be[i][j][k] = (le[i][j][k] ? 0
							: (getCarrierElement(i, j, k) ? 0 : 1.0f));
					// carrier
					carrier[i][j][k] = (_carrierElements[i][j][k] ? 1.0f : 0);
					// carrier
					liquid[i][j][k] = (le[i][j][k] ? 1.0f : 0);
				}
		// create the file and write matrix to it
		File f = new File("/Users/jxavier/results/tube/biofilm.txt");
		java.io.FileWriter fr = new java.io.FileWriter(f);
		fr.write(MultigridUtils.matrixToString(be));
		fr.close();
		// create the file and write matrix to it
		f = new File("/Users/jxavier/results/tube/carrier.txt");
		fr = new java.io.FileWriter(f);
		fr.write(MultigridUtils.matrixToString(carrier));
		fr.close();
		// create the file and write matrix to it
		f = new File("/Users/jxavier/results/tube/liquid.txt");
		fr = new java.io.FileWriter(f);
		fr.write(MultigridUtils.matrixToString(liquid));
		fr.close();
	}

	/**
	 * Auzxiliary method to write out a boolean matrix to file (for debuging
	 * purposes)
	 * 
	 * @param bm
	 * @param filename
	 */
	static public void writeBooleanMatrixToFile(boolean[][][] bm,
			String filename) {
		// convert to biofilm boolean matrix (connected volume
		int n = bm.length;
		int m = bm[0].length;
		int l = bm[0][0].length;
		float[][][] be = new float[n][m][l];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				for (int k = 0; k < l; k++) {
					// carrier
					be[i][j][k] = (bm[i][j][k] ? 1.0f : 0);
				}
		// create the file and write matrix to it
		try {
			File f = new File(filename);
			java.io.FileWriter fr = new java.io.FileWriter(f);
			fr.write(MultigridUtils.matrixToString(be));
			fr.close();
		} catch (IOException e) {
			throw new ModelRuntimeException("Error trying to write " + filename);
		}
	}

	/**
	 * Compute the value of minimum mass if it is set to zero, and return it. So
	 * far this method is not very well designed as it casts the
	 * BiomassContainer to BiomassParticleContainer. Minimum mass value is the
	 * value that a particle with _minimumRadius would have if it was composed
	 * entirely by the particulate species that has the lowest density of all
	 * particulates in the system
	 * 
	 * @return a lower value for the biomass (for levelset determination)
	 */
	private float getMinimumMassValue() {
		if (_minimumMass == 0) {
			// if minimum mass is not yet determined, determine it here
			float minDensity = Float.POSITIVE_INFINITY;
			for (Iterator iter = _model.getParticulateSpecies().iterator(); iter
					.hasNext();) {
				ParticulateSpecies s = (ParticulateSpecies) iter.next();
				minDensity = Math.min(s.getDensity(), minDensity);
			}
			_minimumMass = minDensity
					* BiomassParticle.computeVolume(_minimumRadius);
		}
		return _minimumMass;
	}

	/**
	 * Check if a grid element is empty
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return true if element is empty, false otherwise
	 */
	private boolean elementIsEmpty(int i, int j, int k) {
		BiomassParticleSet ps = shovingGrid[i][j][k];
		if (ps == null)
			return true;
		else if (ps.getNumberOfParticles() == 0)
			return true;
		return false;
	}

	/**
	 * The density that results from the sums of all particle masses devided by
	 * the volume of a shoving grid element
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return density
	 */
	public float getElementDensity(int i, int j, int k) {
		float m = 0;
		BiomassParticleSet ps = shovingGrid[i][j][k];
		// use minimum mass defined in system in case set is empty (or has not
		// been initialized)
		if (ps == null)
			m = getMinimumMassValue();
		else if (ps.getNumberOfParticles() == 0)
			m = getMinimumMassValue();
		else {
			ps.resetIndex();
			while (ps.hasNext())
				m += ps.next().getTotalMass();
		}
		return m / ExtraMath.cube(shovingGridSide);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * org.photobiofilms.model.BiomassContainer#getElementDensity(org.photobiofilms
	 * .model.ContinuousCoordinate)
	 */
	public float getElementDensity(ContinuousCoordinate c) {
		return getElementDensity(snapToGridI(c.x), snapToGridJ(c.y),
				snapToGridK(c.z));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * nl.tudelft.bt.model.BiomassContainer#getConcentrationInElement(nl.tudelft
	 * .bt.model.ContinuousCoordinate,
	 * nl.tudelft.bt.model.multigrid.ParticulateSpecies)
	 */
	public float getConcentrationInElement(ContinuousCoordinate c,
			ParticulateSpecies s) {
		int i = snapToGridI(c.x);
		int j = snapToGridI(c.y);
		int k = snapToGridI(c.z);
		float m = 0;
		BiomassParticleSet ps = shovingGrid[i][j][k];
		// use minimum mass defined in system in case set is empty (or has not
		// been initialized)
		if (ps == null)
			m = getMinimumMassValue();
		else if (ps.getNumberOfParticles() == 0)
			m = getMinimumMassValue();
		else {
			ps.resetIndex();
			while (ps.hasNext())
				m += ps.next().getMassOfSpecies(s);
		}
		return m / ExtraMath.cube(shovingGridSide);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * nl.tudelft.bt.model.BiomassContainer#getEpsInElement(nl.tudelft.bt.model
	 * .ContinuousCoordinate)
	 */
	public float getEpsInElement(ContinuousCoordinate c) {
		int i = snapToGridI(c.x);
		int j = snapToGridI(c.y);
		int k = snapToGridI(c.z);
		float m = 0;
		BiomassParticleSet ps = shovingGrid[i][j][k];
		// use minimum mass defined in system in case set is empty (or has not
		// been initialized)
		if (ps == null)
			m = getMinimumMassValue();
		else if (ps.getNumberOfParticles() == 0)
			m = getMinimumMassValue();
		else {
			ps.resetIndex();
			while (ps.hasNext())
				m += ps.next().getTotalEps();
		}
		return m / ExtraMath.cube(shovingGridSide);
	}

	/**
	 * Determines the position to attach a biomass particle comming from the top
	 * boundary of the system performing a random walk. Random walk is performed
	 * on the coordinates of the shoving grid
	 * 
	 * @return the position to attach the particle
	 */
	private ContinuousCoordinate determineAttachmentPositionFromRandomWalk() {
		// choose starting point on top
		int n = _n - 1;
		int m = (int) (_model.getRandom() * _m);
		int l = (int) (_model.getRandom() * _l);
		// Do the Random-walk until a grid element containing biofilm or close
		// to the solid surface is reached
		do {
			// a direction for movement is chosen from 5 options (all
			// directions except top)
			int dir = (int) (_model.getRandom() % 5);
			// Now move ...
			switch (dir) {
			case 0:
				// down
				n--;
				break;
			case 1:
				// right (periodic boundary)
				m = (m == _m - 1 ? 0 : m - 1);
				break;
			case 2:
				// left (periodic boundary)
				m = (m == 0 ? _m - 1 : m + 1);
				break;
			case 3:
				// backward (periodic boundary)
				l = (l == _l - 1 ? 0 : l - 1);
				break;
			case 4:
				// forward (periodic boundary)
				l = (l == 0 ? _l - 1 : l + 1);
				break;
			}
			// Check if there is biofilm in this grid element
		} while ((n > 0) & (elementIsEmpty(n, m, l)));
		// convert grid coordinates to continuous coordinates
		return getGridElementCenter(n, m, l);
	}

	/**
	 * Compute the biovolume from the mapping representation of the biomass
	 * 
	 * @return the total volume occupied by the biofilm (biovolume)
	 */
	public float computeBiovolume() {
		// count the voxels occupied by the biofilm
		int n = 0;
		// get the information concerning all the biomass
		for (int i = 0; i < _n; i++)
			for (int j = 0; j < _m; j++)
				for (int k = 0; k < _l; k++)
					if (_map[i][j][k])
						n++;
		// convert to volume
		return ExtraMath.cube(shovingGridSide) * ((float) n);
	}

	/**
	 * Tests if a grid node belongs to the solid carrier material. This method
	 * should be overloaded when new carrier shapes are implemented, e.g. the
	 * tubular reactor, in which the carrier is located all arround the biofilm
	 * and liquid in a circular shape. For the default BiomassParticleContainer,
	 * which implements a planar geometry, this method always returns false, as
	 * the carrier is located below the computational volume, and no entries in
	 * the matrix actually represent the carrier. The same for the granule
	 * geometry, meaning that GranuleBiomassContainer does not need to override
	 * this method
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return true if i, j, k is the carrier, false otherwise
	 */
	protected boolean isCarrier(int i, int j, int k) {
		return false;
	}

	/**
	 * get the value of carrier elements matrix and initialize it if it is not
	 * initialized yet
	 * 
	 * @param i
	 * @param j
	 * @param k
	 * @return the value at i, j, k
	 */
	private boolean getCarrierElement(int i, int j, int k) {
		if (_carrierElements == null) {
			// create and initialize the carrier elements matrix
			_carrierElements = new boolean[_n + 10][_m][_l];
			for (int i2 = 0; i2 < _n + 10; i2++)
				for (int j2 = 0; j2 < _m; j2++)
					for (int k2 = 0; k2 < _l; k2++)
						_carrierElements[i2][j2][k2] = isCarrier(i2, j2, k2);
		}
		return _carrierElements[i][j][k];
	}
	
	/**
	 * Turns off the removal of disconnected cells
	 */
	public void turnSloughingOff() {
		_sloughinfIsOn = false;
	}
}