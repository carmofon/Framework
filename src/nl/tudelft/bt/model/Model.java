package nl.tudelft.bt.model;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;

import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.multigrid.boundary_layers.BoundaryLayer;
import nl.tudelft.bt.model.parsers.SimulationResultParser;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.BiomassParticleContainer;
import nl.tudelft.bt.model.particlebased.granule.GranuleBiomassParticleContainer;
import nl.tudelft.bt.model.particlebased.tube.TubeBiomassParticleContainer;
import nl.tudelft.bt.model.reaction.Reaction;
import nl.tudelft.bt.model.timeconstraint.TimeStepConstraint;
import nl.tudelft.bt.model.util.ExtraMath;
import nl.tudelft.bt.model.util.UnlimitedFloatArray;
import nl.tudelft.bt.model.apps.output.SimulationResultsWriter;
import nl.tudelft.bt.model.apps.output.SolidsConcentrationWriter;
import nl.tudelft.bt.model.bulkconcentrations.BulkConcentration;
import nl.tudelft.bt.model.bulkconcentrations.ConstantBulkConcentration;
import nl.tudelft.bt.model.detachment.levelset.DetachedBiomassContainer;
import nl.tudelft.bt.model.exceptions.*;

/**
 * Main model class. Implements singleton design pattern to impose a unique
 * instance of model. 22042004 new change
 * 
 * @author Jo?o Xavier (j.xavier@tnw.tudelft.nl)
 */
/**
 * @author xavierj
 * 
 */
public class Model implements Serializable {
	// the singleton instance of model (exists 1 and only 1)
	private static Model _model = new Model();

	private int _seed; // random number generator seed

	private Random _rng; // random number generator

	// computation
	private int _dimensionality;

	private int _iterationCounter;

	private float _time;

	private long _realInitialTime;

	private float _timeStep;

	private float _compulsoryTimeStep;

	private float _maximumTimeStep;

	private float _minimumTimeStep;

	private boolean _writeTimedWriters; // flag for timed writers

	private boolean _overrideTimeStep;

	private ArrayList<TimeStepConstraint> _timeStepConstraints;

	// attributes for simulation finishing criterium
	private float _finishIterationTime;

	private float _finishIterationCount;

	private float _maxBiovolume;

	private float _maxRunLength;

	private boolean _finishSimulation;

	private UnlimitedFloatArray _timeSeries;

	public ContinuousCoordinate systemSize;

	public float referenceSystemSide;

	private float _systemGridSide;

	private int _numberOfGridNodes;

	public BiomassContainer biomassContainer;

	private BiomassSpecies[] _biomassSpecies;

	private SoluteSpecies[] _soluteSpecies;

	private ParticulateSpecies[] _particulateSpecies;

	private ArrayList<Reaction> _reactions;

	private BoundaryLayer _boundaryLayer;

	private float _maxBiofilmHeight;

	private float _finalBiofilmHeight;

	// detached biomass
	private DetachedBiomassContainer _detachedBiomass;

	// global system scale
	private float _carrierArea;

	private float _residenceTime;

	private float _reactorVolume;

	/**
	 * Model constructor (private to implement Singleton Model), initializes
	 * random seed to its default value.
	 */
	private Model() {
		// initialize the random number generator
		// to default seed value
		// predifined
		// _seed = 1007;
		// random
		_seed = (int) System.currentTimeMillis();
		reset();
	}

	/**
	 * Returns the Singleton instance of model
	 * 
	 * @return the singleton instace of model
	 */
	public static Model model() {
		return _model;
	}

	/**
	 * Resets the model.
	 */
	public void reset() {
		// initialize the random number generator
		_rng = new java.util.Random(_seed); // use this for fixed seed
		// TODO delete:
		System.out.println("seed " + _seed);
		// _rng = new java.util.Random(); // use this for variable seed
		_dimensionality = 0;
		_iterationCounter = 0;
		_realInitialTime = System.currentTimeMillis(); // save the initial time
		// for benchmark purposes
		_time = 0;
		_timeStep = 0;
		_compulsoryTimeStep = 0;
		_maximumTimeStep = 0;
		_minimumTimeStep = 0;
		_writeTimedWriters = true;
		_overrideTimeStep = false;
		_finishIterationTime = Float.POSITIVE_INFINITY;
		_finishIterationCount = Float.POSITIVE_INFINITY;
		_maxBiovolume = Float.POSITIVE_INFINITY;
		_maxRunLength = Float.POSITIVE_INFINITY;
		_finishSimulation = false;
		_timeSeries = new UnlimitedFloatArray();
		_timeSeries.add(_time);
		referenceSystemSide = 0;
		_maxBiofilmHeight = 0;
		_finalBiofilmHeight = Float.POSITIVE_INFINITY;
		systemSize = null;
		_timeStepConstraints = new ArrayList<TimeStepConstraint>();
		// delete all references to bacteria and the shoving grid:
		biomassContainer = null;
		// removes all discrete multigrid data on species
		_biomassSpecies = null;
		_soluteSpecies = null;
		// the biofilm reactor propetries
		_residenceTime = 0;
		_carrierArea = 0;
		_reactorVolume = 0;
		_detachedBiomass = new DetachedBiomassContainer(_timeSeries);
		// Initialize the reactions container
		_reactions = new ArrayList<Reaction>();
	}

	/**
	 * Set the value for the random number generator seed
	 * 
	 * @param i
	 *            value to set seed to
	 */
	public void setSeed(int i) {
		_seed = i;
	}

	/**
	 * The only method o provides a random number that should be used anywhere
	 * where random number generation is require in the application.
	 * 
	 * @return random float between 0 and 1
	 * @see Random
	 */
	public float getRandom() {
		return _rng.nextFloat();
	}

	/**
	 * The only method o provides a random number that should be used to
	 * generate a normally distributed random variavle.
	 * 
	 * @return random with average 0 and standard deviation 1
	 * @see Random
	 */
	public float getRandomFromNormalDistribution() {
		return (float) (_rng.nextGaussian());
	}

	/**
	 * Builds the spatial system
	 * 
	 * @param d
	 *            dimensionality (2D or 3D)
	 * @param s
	 *            system size
	 * @param g
	 *            grid size for solute concentration field
	 * @throws InvalidValueException
	 *             id ((d != 2) && (d != 3)) and if g is not of form 2^n + 1
	 */
	public void buildSystem(int d, float s, int g) throws InvalidValueException {
		if ((d != 2) && (d != 3))
			throw new InvalidValueException("dimensionality value (" + d
					+ ") not valid");
		_dimensionality = d;
		systemSize = new ContinuousCoordinate(s, s, (d == 3 ? s : 0));
		referenceSystemSide = s;
		_systemGridSide = s / g;
		_numberOfGridNodes = g;
		_maxBiofilmHeight = s;
		// define the grid resolution
		MultigridVariable.setGrid(g, g, (d == 3 ? g : 1));
	}

	/**
	 * Expands the system along the vertical direction. If expansionOrder is 0
	 * it keeps the same size. If expansionOrder is 1 it doubles the size. If
	 * expansionOrder is 2 it quadruples the size. etc.
	 * 
	 * @param expansionOrder
	 */
	public void expandSystemVertically(int expansionOrder) {
		int newVerticalGridSize = ExtraMath.exp2(MultigridUtils
				.order(_numberOfGridNodes) + expansionOrder) + 1;
		float newVerticalSize = newVerticalGridSize * _systemGridSide;
		systemSize = new ContinuousCoordinate(newVerticalSize,
				referenceSystemSide,
				(_dimensionality == 3 ? referenceSystemSide : 0));
		// define the grid resolution
		MultigridVariable.setGrid(newVerticalGridSize, _numberOfGridNodes,
				(_dimensionality == 3 ? _numberOfGridNodes : 1));
	}

	/**
	 * Get value of time step
	 * 
	 * @return time step
	 */
	public float getTimeStep() {
		return _timeStep;
	}
	
	/**
	 * Return the maximum time step
	 * 
	 * @return _maximumTimeStep
	 */
	public float getMaximumTimeStep() {
		return _maximumTimeStep;
	}

	/**
	 * Set the size for biofilm detachment (thickness if flate geometry biofilm,
	 * radius if granule or colony
	 * 
	 * @param h
	 *            maximum biofilm size (units of length)
	 */
	public void setVerticalCutoffSize(float h) {
		_maxBiofilmHeight = h;
	}

	/**
	 * Get the biofilm detachment height
	 * 
	 * @return height biofilm detachment height
	 */
	public float getMaximumBiofilmHeight() {
		return _maxBiofilmHeight;
	}

	/**
	 * Get the current biofilm height
	 * 
	 * @return current biofilm height
	 */
	public float getCurrentBiofilmHeight() {
		return biomassContainer.getTopBiofilmHeight();
	}

	/**
	 * Get the current biovolume. Biovolume is computed from the multigrid
	 * 
	 * @return current biovolume
	 */
	public float getCurrentBiovolume() {
		return biomassContainer.computeBiovolume();
		// return MultigridVariable.computeBiovolume(_particulateSpecies);
	}

	/**
	 * @return the total fixed biomass [kg/m2]
	 */
	public float getCurrentTotalBiomass() {
		return biomassContainer.getTotalBiomass();
	}

	/**
	 * Set the array of bacteria species, chemical species and the boundary
	 * layer object that will be used by the multigrid solver.
	 * 
	 * @param b
	 *            array of bacteria species
	 * @param c
	 *            array of bacteria species
	 * @param bl
	 *            boundary layer object
	 */
	public void setupDiffusionReactionSystem(BiomassSpecies[] b,
			SoluteSpecies[] c, BoundaryLayer bl) {
		_biomassSpecies = b;
		_soluteSpecies = c;
		_boundaryLayer = bl;
		// get the array of fixed substances from the biomass species array
		Collection<ParticulateSpecies> col = new ArrayList<ParticulateSpecies>();
		for (int i = 0; i < _biomassSpecies.length; i++) {
			Collection<ParticulateSpecies> clocal = _biomassSpecies[i]
					.getFixedSpeciesAsCollection();
			col.removeAll(clocal);
			col.addAll(clocal);
		}
		_particulateSpecies = new ParticulateSpecies[col.size()];
		int i = 0;
		for (Iterator<ParticulateSpecies> iter = col.iterator(); iter.hasNext();) {
			_particulateSpecies[i++] = (ParticulateSpecies) iter.next();
		}
	}

	/**
	 * Add a reaction to the reaction container. This method is called directly
	 * by the Reaction contructor, so that each reaction created is
	 * automatically stored
	 * 
	 * @param r
	 *            reactioin to add to model
	 */
	public void addReaction(Reaction r) {
		_reactions.add(r);
	}

	/**
	 * Creates the bacteria container as a spherical particle container.
	 * 
	 * @param maximumRadius
	 *            maximum radius of biomass particles (for division)
	 * @param minimumRadius
	 *            minimum radius of biomass particles (for removal)
	 * @param kshov
	 *            shoving multiplier
	 * @param fshov
	 *            fraction of remaining particles
	 * @param rdetach
	 *            detachment rate [10^-15g/um^2/h]
	 */
	public void buildBacteriaContainer(float maximumRadius,
			float minimumRadius, float kshov, float fshov) {
		biomassContainer = new BiomassParticleContainer(maximumRadius,
				minimumRadius, kshov, fshov);
	}

	/**
	 * Creates a granule biomass container as a spherical particle container.
	 * 
	 * @param maximumRadius
	 *            maximum radius of biomass particles (for division)
	 * @param minimumRadius
	 *            minimum radius of biomass particles (for removal)
	 * @param kshov
	 *            shoving multiplier
	 * @param fshov
	 *            fraction of remaining particles
	 * @param rdetach
	 *            detachment rate [10^-15g/um^2/h]
	 */
	public void buildGranuleBiomassContainer(float maximumRadius,
			float minimumRadius, float kshov, float fshov) {
		biomassContainer = new GranuleBiomassParticleContainer(maximumRadius,
				minimumRadius, kshov, fshov);
	}

	/**
	 * Create a new container for tube geometry
	 * 
	 * @param maximumRadius
	 * @param minimumRadius
	 * @param kshov
	 * @param fshov
	 */
	public void buildTubeBiomassContainer(float maximumRadius,
			float minimumRadius, float kshov, float fshov, float tubeRadius) {
		// check validity of size of tube radius
		if (referenceSystemSide < tubeRadius
				/ (0.5f - 1f / (float) _numberOfGridNodes))
			throw new ModelRuntimeException("Radius of tube (" + tubeRadius
					+ ") not agreeing with system size (" + referenceSystemSide
					+ ")");
		// create the new instance of the biomass tube container
		biomassContainer = new TubeBiomassParticleContainer(maximumRadius,
				minimumRadius, kshov, fshov, tubeRadius);
	}

	/**
	 * Inoculate system at random positions with arbitrary number of cells of
	 * each defined species. Mass of cell is the average mass defined in species
	 * 
	 * @param n
	 *            number of cells to inoculate for each species
	 * @throws NonMatchingNumberException
	 *             if length of n does not match number of defined bacteria
	 *             species
	 */
	public void inoculateRandomlyMultispecies(int[] n)
			throws NonMatchingNumberException {
		if (n.length != _biomassSpecies.length)
			throw new NonMatchingNumberException(
					"trying to inoculate system with " + n.length
							+ " when it has " + _biomassSpecies.length
							+ " species");
		// Inoculate cells randomly (shhuffling order)
		int nspecies = n.length;
		int[] numberOfCellsToPlace = new int[n.length];
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = n[i];
		while (true) {
			// get the random value
			int r = (int) Math.floor(getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					// set the center to a random position at the surface
					float y = getRandom() * systemSize.y;
					float z = getRandom() * systemSize.z;
					// add to the bacteria list
					biomassContainer
							.placeBiomassAt(_biomassSpecies[r], 0, y, z);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}
		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();
	}

	/**
	 * Place a single cell in center of system with masses defined by the
	 * default array of masses of the biomass species
	 * 
	 * @param s
	 *            species of biomass to add
	 * @param x
	 *            height from solid surface [um]
	 * @param y
	 *            horizontal coordinate [um]
	 * @param z
	 *            horizontal coordinate [um]
	 * @return reference to particle added
	 */
	public BiomassParticle placeSingleCellInCenter(BiomassSpecies s) {
		float x = 0;
		float y = systemSize.y * 0.5f;
		float z = systemSize.z * 0.5f;
		return biomassContainer.placeBiomassAt(s, x, y, z);
	}

	/**
	 * Place biomass in the system with masses defined by the default array of
	 * masses of the biomass species
	 * 
	 * @param s
	 *            species of biomass to add
	 * @param x
	 *            height from solid surface [um]
	 * @param y
	 *            horizontal coordinate [um]
	 * @param z
	 *            horizontal coordinate [um]
	 * @return reference to particle added
	 */
	public BiomassParticle placeBiomass(BiomassSpecies s, float x, float y,
			float z) {
		return biomassContainer.placeBiomassAt(s, x, y, z);
	}

	/**
	 * @param n
	 * @throws NonMatchingNumberException
	 */
	public void inoculateRandomlyMultispeciesInGranule(int[] n)
			throws NonMatchingNumberException {
		int nspecies = n.length;
		if (nspecies != _biomassSpecies.length)
			throw new NonMatchingNumberException(
					"trying to inoculate system with " + n.length
							+ " when it has " + _biomassSpecies.length
							+ " species");
		// copy the n array
		int[] numberOfCellsToPlace = new int[n.length];
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = n[i];
		while (true) {
			// get the random value
			int r = (int) Math.floor(getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					((GranuleBiomassParticleContainer) biomassContainer)
							.placeInoculumParticle(_biomassSpecies[r]);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}

		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();
	}
	/**
	 * Inoculates particles in a granule at a random location of given distance 
	 * from center
	 * 
	 * @param n
	 * @param radius
	 * @throws NonMatchingNumberException
	 */
	public void inoculateRandomlyMultispeciesInGranuleAtRadius(int[] n, float radius)
		throws NonMatchingNumberException {
		int nspecies = n.length;
		if (nspecies != _biomassSpecies.length)
			throw new NonMatchingNumberException(
					"trying to inoculate system with " + n.length
							+ " when it has " + _biomassSpecies.length
							+ " species");
		// calculate a random coordinate at distance r from center at which to inoculate
		float phi =  (_dimensionality == 3) ?
			_model.getRandom() * ExtraMath.PI : 0.5f * ExtraMath.PI;
		float theta = _model.getRandom() * ExtraMath.PI2;
		float bx = radius * (float) (Math.cos(theta) * Math.sin(phi));
		float by = radius * (float) (Math.sin(theta) * Math.sin(phi));
		float bz = radius * (float) (Math.cos(phi));
		// copy the n array
		int[] numberOfCellsToPlace = new int[n.length];
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = n[i];
		while (true) {
			int r = (int) Math.floor(getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					BiomassSpecies species = _biomassSpecies[r];
					// place a new cell in the center
					ContinuousCoordinate center = 
						((GranuleBiomassParticleContainer) biomassContainer).get_centerOfComputationalVolume();
					BiomassParticle b = placeBiomass(species, center.x, center.y, center.z);
					// compute a random displacement as long as the cell's radius
					float d = b.getRadius();
					float dphi = (_dimensionality == 3) ? 
							_model.getRandom() * ExtraMath.PI : 0.5f * ExtraMath.PI;
					float dtheta = _model.getRandom() * ExtraMath.PI2;
					float dx = d * (float) (Math.cos(dtheta) * Math.sin(dphi));
					float dy = d * (float) (Math.sin(dtheta) * Math.sin(dphi));
					float dz = d * (float) (Math.cos(dphi));
					// move to the location of the inoculum + random offset
					b.move(bx + dx, by + dy, bz + dz);
				}
			}
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}
		spreadByShovingCompletely();
	}

	/**
	 * Place randonmly cells of each species inside a central region of radius r
	 * 
	 * @param n
	 * @param rad
	 * @throws NonMatchingNumberException
	 */
	public void inoculateRandomlyMultispeciesInColony(int[] n, float rad)
			throws NonMatchingNumberException {
		int nspecies = n.length;
		if (nspecies != _biomassSpecies.length)
			throw new NonMatchingNumberException(
					"trying to inoculate system with " + n.length
							+ " when it has " + _biomassSpecies.length
							+ " species");
		// copy the n array
		int[] numberOfCellsToPlace = new int[n.length];
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = n[i];
		while (true) {
			// get the random value
			int r = (int) Math.floor(getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					((GranuleBiomassParticleContainer) biomassContainer)
							.placeInoculumParticleInsideRadius(
									_biomassSpecies[r], rad);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}

		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();

	}

	/**
	 * Inoculates multiple species at fixed distance from center of simulation
	 * space
	 * 
	 * @param n
	 *            array of number of cells to place
	 * @param rad
	 *            fixed radius at which to place the cells
	 * @throws NonMatchingNumberException
	 *             if number of cell number does not match number of species in
	 *             system
	 */
	public void inoculateRandomlyMultispeciesAtRadius(int[] n, float rad)
			throws NonMatchingNumberException {
		int nspecies = n.length;
		if (nspecies != _biomassSpecies.length)
			throw new NonMatchingNumberException(
					"trying to inoculate system with " + n.length
							+ " when it has " + _biomassSpecies.length
							+ " species");
		// copy the n array
		int[] numberOfCellsToPlace = new int[n.length];
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = n[i];
		while (true) {
			// get the random value
			int r = (int) Math.floor(getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					((GranuleBiomassParticleContainer) biomassContainer)
							.placeInoculumParticleAtRadius(_biomassSpecies[r],
									rad);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}

		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();

	}

	/**
	 * Place randonmly cells of each species everythere in space
	 * 
	 * @param n
	 * @param rad
	 * @throws NonMatchingNumberException
	 */
	public void inoculateRandomlyMultispeciesEverywhere(int[] n)
			throws NonMatchingNumberException {
		int nspecies = n.length;
		if (nspecies != _biomassSpecies.length)
			throw new NonMatchingNumberException(
					"trying to inoculate system with " + n.length
							+ " when it has " + _biomassSpecies.length
							+ " species");
		// copy the n array
		int[] numberOfCellsToPlace = new int[n.length];
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = n[i];
		while (true) {
			// get the random value
			int r = (int) Math.floor(getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					// set the center to a random position at the surface
					float x = getRandom() * systemSize.x;
					float y = getRandom() * systemSize.y;
					float z = getRandom() * systemSize.z;
					// add to the bacteria list
					biomassContainer
							.placeBiomassAt(_biomassSpecies[r], x, y, z);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}

		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();

	}

	/**
	 * @param n
	 * @throws NonMatchingNumberException
	 */
	public void inoculateRandomlyMultispeciesInTube(int[] n)
			throws NonMatchingNumberException {
		int nspecies = n.length;
		if (nspecies != _biomassSpecies.length)
			throw new NonMatchingNumberException(
					"trying to inoculate system with " + n.length
							+ " when it has " + _biomassSpecies.length
							+ " species");
		// copy the n array
		int[] numberOfCellsToPlace = new int[n.length];
		for (int i = 0; i < nspecies; i++)
			numberOfCellsToPlace[i] = n[i];
		while (true) {
			// get the random value
			int r = (int) Math.floor(getRandom() * nspecies);
			if (r < nspecies) {
				if (numberOfCellsToPlace[r] > 0) {
					numberOfCellsToPlace[r]--;
					((TubeBiomassParticleContainer) biomassContainer)
							.placeInoculumParticle(_biomassSpecies[r]);
				}
			}
			// do not break if any species has not reached 0 yet
			boolean done = true;
			for (int i = 0; i < nspecies; i++)
				done &= (numberOfCellsToPlace[i] == 0);
			if (done)
				break;
		}
		// After innoculating, spread the biomass completelly for homogenizing
		// the innoculum
		spreadByShovingCompletely();
	}

	/**
	 * Update the concentration of biomass in the diffusion-reaction data
	 * matrices.
	 */
	public void updateBioDiscreteData() {
		// reset values of all BacteriaSpecies
		for (int i = 0; i < _particulateSpecies.length; i++) {
			_particulateSpecies[i].resetDiscreteMatrix();
		}
		// iterate through array of bacteria and increase
		// concentration value in corresponding BacteriaSpecies
		biomassContainer.addContributionsToSpeciesDiscreteData();
		// update boundary values of all BacteriaSpecies
		for (int i = 0; i < _particulateSpecies.length; i++) {
			_particulateSpecies[i].refreshBoundaryConditions();
		}
	}

	/**
	 * Compute the initial concentration fields and initializaes the global
	 * rates time series
	 * 
	 * @throws MultigridSystemNotSetException
	 */
	public void initializeSoluteConcentrations()
			throws MultigridSystemNotSetException {
		updateBioDiscreteData();
		//TODO put multigrid back on
		MultigridVariable.solve(_soluteSpecies, _particulateSpecies,
				_boundaryLayer);
//		MultigridVariable.initializeSolveByRelax(_soluteSpecies, _particulateSpecies,
//		_boundaryLayer);
		// perform the global reactor mass balances
		for (int i = 0; i < _soluteSpecies.length; i++) {
			_soluteSpecies[i].initializeGlobalRateSeries();
		}
	}

	/**
	 * Updates the bulk concentrations of solutes and the discrete
	 * representation of the solids concentration fields and then determines the
	 * concentration fields of all solutes by solving the diffusion/reaction
	 * PDE's using multigrid
	 * 
	 * @throws MultigridSystemNotSetException
	 */
	public void updateSoluteConcentrations()
			throws MultigridSystemNotSetException {
		// perform the global reactor mass balances
		for (int i = 0; i < _soluteSpecies.length; i++) {
			_soluteSpecies[i].updateBulkConcentrationAndRateSeries(_timeStep);
		}
		updateBioDiscreteData();
		solveWithMultigrid();
		// NOTE to solve by relaxation use:
		// MultigridVariable.solveByRelax(
		// _chemicalSpecies,
		// _fixedSpecies,
		// _boundaryLayer);
	}

	/**
	 * Determines the concentration fields of all solutes by solving the
	 * diffusion/reaction PDE's using multigrid
	 * 
	 * @throws MultigridSystemNotSetException
	 * 
	 */
	public void solveWithMultigrid() throws MultigridSystemNotSetException {
//		MultigridVariable.solve(_soluteSpecies, _particulateSpecies,
//				_boundaryLayer);
		//TODO put multigrid back on
//		MultigridVariable.solve(_soluteSpecies, _particulateSpecies,
//				_boundaryLayer);
		MultigridVariable.solveByRelax(_soluteSpecies, _particulateSpecies,
		_boundaryLayer);

	}

	/**
	 * Compute the time step for growth of biomass and global mass balances.
	 * Uses TimeStepConstraint in some operations to compare timesteps and keep
	 * track of which is the constraint setting the time step.
	 */
	private void computeTimeStep() throws GrowthProblemException {
		TimeStepConstraint presentTimeStepConstraint = new TimeStepConstraint();
		if (!_overrideTimeStep) {
			// ///////determine the variable time step
			// Constraint 1: minimum particle doubling time
			// Has to be the first to be computed!!!
			TimeStepConstraint particleConstraint = biomassContainer
					.getParticleTimeConstraint();

			// Constraint 2: global volume change
			TimeStepConstraint biofilmSizeConstraint = biomassContainer
					.getGlobalSizeConstraint();
			// NOTE comment previous line and uncomment following if no biofilm
			// growth constraint is required
			// float biofilmSizeConstraint = Float.POSITIVE_INFINITY;
			// Constraint 3: the computation of bulk concentration time
			// constraints
			TimeStepConstraint bulkConcentrationConstraint = new TimeStepConstraint();
			for (int i = 0; i < _soluteSpecies.length; i++) {
				bulkConcentrationConstraint = TimeStepConstraint.getMinimum(
						bulkConcentrationConstraint,
						_soluteSpecies[i].getBulkSpeciesTimeStepConstraint());
			}
			// other timestep limiting criteria may be added here...
			presentTimeStepConstraint = TimeStepConstraint.getMinimum(
					particleConstraint, biofilmSizeConstraint,
					bulkConcentrationConstraint);
			_timeStep = presentTimeStepConstraint.getTimeStep();
			// if time step is infinity it means biofilm is dead!!!
			// break iteration
			if (_timeStep == Float.POSITIVE_INFINITY)
				throw new GrowthProblemException();
			// If a value for the compulsory time step is defined, determine if
			// compulsory time step should be preformed instead of using the
			// variable time step just determined
			if (_compulsoryTimeStep != 0) {
				// compute the number of compulsory steps already passed
				double nCompulsories = Math.floor(_time / _compulsoryTimeStep);
				// compute the nexte compulsory step
				double nextCompulsory = (nCompulsories + 1)
						* _compulsoryTimeStep;
				// compute the time remaining to next compulsory
				float timeToCompulsory = (float) nextCompulsory - _time;
				if (timeToCompulsory < _timeStep * 1e-3f) {
					nextCompulsory += _compulsoryTimeStep;
					timeToCompulsory += _compulsoryTimeStep;
				}
				// Due to float precision, this operation may cause the
				// simulation to get stuck if timeToCompulsory returns a very
				// small value. Therefore, a check is introduced:
				// _timeStep is used not only if:
				// timeToCompulsory < _timeStep
				// but also if:
				// timeToCompulsory > _timeStep * 1e-3
				// or
				// _writeTimedWriters = false
				if ((timeToCompulsory <= _timeStep)
						| (nextCompulsory == (_time + _timeStep))) {
					// the check (nextCompulsory == (_time + _timeStep))
					// was also included because of precision issues
					_timeStep = timeToCompulsory;
					presentTimeStepConstraint.setTimeStep(_timeStep);
					presentTimeStepConstraint.setName("Compulsory time step");
					// set the flag to write in this iteration
					_writeTimedWriters = true;
				} else {
					// set the flag NOT to write in this iteration
					_writeTimedWriters = false;
				}
			}
		}
		// If a maximum time step is defined, truncate at that value
		if (_maximumTimeStep > 0) {
			if (_maximumTimeStep < _timeStep) {
				_timeStep = _maximumTimeStep;
				presentTimeStepConstraint.setTimeStep(_maximumTimeStep);
				presentTimeStepConstraint.setName("Maximum time step");
			}
		}
		// If a minimum time step is defined, lift to at least that value
		if (_minimumTimeStep > 0) {
			_timeStep = Math.max(_timeStep, _minimumTimeStep);
			if (_minimumTimeStep > _timeStep) {
				_timeStep = _minimumTimeStep;
				presentTimeStepConstraint.setTimeStep(_minimumTimeStep);
				presentTimeStepConstraint.setName("Minimum time step");
			}
		}
		// update counters
		incrementIteration();
		_time += _timeStep;
		_timeSeries.add(_time);
		_timeStepConstraints.add(presentTimeStepConstraint);
	}

	/**
	 * Increment the iteration counter by 1
	 */
	public void incrementIteration() {
		_iterationCounter++;
	}

	/**
	 * Perform the growth and division part of the iterative step
	 */
	public void performGrowthAndDivision() throws ModelException {
		// reset the reaction rates for all the reactions
		for (Iterator<Reaction> iter = _reactions.iterator(); iter.hasNext();) {
			Reaction r = iter.next();
			r.resetGlobalReactionRate();
		}
		// compute the timestep for present iteration
		computeTimeStep();
		// Here are the growth steps
		// 1. GROWTH step: only updates masses of all bacteria
		biomassContainer.grow(_timeStep);

		// TODO the following is preliminary code for the attachment of biomass
		// dependent o n the gorwth rate of the bacteria species:
		// ((BiomassParticleContainer) biomassContainer).attach();
		// 2. DEATH step: remove cells with mass below a critical level
		biomassContainer.removeDeadBiomass();
		// 3. DIVISION step: splits bacteria over critial mass
		biomassContainer.divideAndExcreteEPS();
		// 4. DETACHMENT step: performs erosion, sloughing and removes
		// out-of-bounds particles
		//detach();
	}

	/**
	 * Spread biomass after growth and division
	 */
	public void spread() {
		biomassContainer.spread();
	}

	public void spreadByShovingCompletely() {
		biomassContainer.spreadByShoving();
	}

	/**
	 * Perform a single spreading step (for testing or illustration)
	 */
	public void performSpreadingStep() {
		biomassContainer.performSpreadingStep();
	}

	/**
	 * Perform the biomass shrinking by pressure
	 */
	public void performSpreadingByPulling() {
		biomassContainer.performSpreadingByPulling();
	}

	/**
	 * Remove all biomass to be dettached
	 */
	public void detach() throws ModelException {
		biomassContainer.removeDetachedBiomass();
	}

	/**
	 * Get the biomass at present state as collection of Particles
	 * 
	 * @return collection of biomass particles
	 */
	public Collection getBiomassAsBiomassParticleCollection() {
		return biomassContainer.getBiomassAsBiomassParticleCollection();
	}

	/**
	 * Write the concetrations of all solutes to text files in a subdirectory of
	 * directory dirName. The subdirectory 'iteration*' is created, where * is a
	 * variable number numbering of the iteration. Files are written using
	 * method finestGridToFile.
	 * 
	 * @param dirName
	 *            directory to write files to
	 * @throws Exception
	 *             if subdircannot be created or in case of failure to write
	 *             files
	 */
	public void writeSoluteConcentrationsToFile(String dirName)
			throws Exception {
		// create directory for this iteration
		File dir = new java.io.File(dirName + "/iteration"
				+ getFormatedIterationNumber());
		try {
			dir.mkdir();
		} catch (SecurityException e) {
			throw new SystemEditViolationException("directory " + dir
					+ " could not be created for security reasons");
		}
		for (int i = 0; i < _soluteSpecies.length; i++) {
			_soluteSpecies[i].finestGridToFile(dir.toString() + "/");
		}
	}

	/**
	 * Write the concetrations of all solids to text files in a subdirectory of
	 * directory dirName. The subdirectory 'iteration*' is created, where * is a
	 * variable number numbering of the iteration. Files are written using
	 * method finestGridToFile.
	 * 
	 * @param dirName
	 *            directory to write files to
	 * @throws Exception
	 *             if subdircannot be created or in case of failure to write
	 *             files
	 */
	public void writeSolidConcentrationsToFile(String dirName) throws Exception {
		// create directory for this iteration
		File dir = new java.io.File(dirName + "/iteration"
				+ getFormatedIterationNumber());
		try {
			dir.mkdir();
		} catch (SecurityException e) {
			throw new SystemEditViolationException("directory " + dir
					+ " could not be created for security reasons");
		}
		for (int i = 0; i < _particulateSpecies.length; i++) {
			_particulateSpecies[i].finestGridToFile(dir.toString() + "/");
		}
	}

	/**
	 * Write the values in a multigrid variable to text files in a subdirectory
	 * of directory dirName. The subdirectory 'iteration*' is created, where *
	 * is a variable number numbering of the iteration. Files are written using
	 * method finestGridToFile.
	 * 
	 * @param dirName
	 *            directory to write files to
	 * @throws Exception
	 *             if subdircannot be created or in case of failure to write
	 *             files
	 */
	public void writeMultigridVariableToFile(String dirName,
			MultigridVariable var) throws Exception {
		// create directory for this iteration
		File dir = new java.io.File(dirName + "/iteration"
				+ getFormatedIterationNumber());
		try {
			dir.mkdir();
		} catch (SecurityException e) {
			throw new SystemEditViolationException("directory " + dir
					+ " could not be created for security reasons");
		}
		var.finestGridToFile(dir.toString() + "/");
	}

	/**
	 * Read the discrete concentrations of all solutes from file
	 * 
	 * @param inputDirectory
	 *            the directory where solute concentrations are
	 * @param inputIterationNumber
	 *            the iteration number to input
	 */
	public void readBiomassMatricesFromDisk(String inputDirectory,
			int inputIterationNumber) {
		String dir = inputDirectory + "/" + SolidsConcentrationWriter.SOLIDDIR
				+ "/iteration" + getFormatedNumber(inputIterationNumber);
		for (int i = 0; i < _particulateSpecies.length; i++) {
			_particulateSpecies[i].readDiscreteConcentrationFromDir(dir,
					inputIterationNumber);
		}
	}

	public void readBulkConcentrationsFromDisk(String inputDirectory,
			int inputIterationNumber) {
		String resultFile = inputDirectory + "/"
				+ SimulationResultsWriter.FILENAME;
		SimulationResultParser p = new SimulationResultParser(resultFile);
		for (int i = 0; i < _soluteSpecies.length; i++) {
			float v = p.getValueForBulkConcentration(_soluteSpecies[i]);
			BulkConcentration bc = new ConstantBulkConcentration(v);
			_soluteSpecies[i].setBulkConcentration(bc);
		}

	}

	/**
	 * Writes the detachment level set matrix to file, before detachment is
	 * carried out
	 * 
	 * @param dirName
	 * @throws Exception
	 */
	public void writeDetachmentLevelSet(String dirName) throws Exception {
		biomassContainer.writeDetachmentLevelSet(dirName + "/iteration"
				+ getFormatedIterationNumber());
	}

	/**
	 * Write the biofilm front boolean matrix to file
	 * 
	 * @param dirName
	 * @throws IOException
	 *             if write to file fails
	 */
	public void writeBiofilmFront(String dirName) throws IOException {
		biomassContainer.writeBiofilmFront(dirName + "/iteration"
				+ getFormatedIterationNumber());
	}

	/**
	 * Writes the positions of bacteria as particles
	 * 
	 * @param dirName
	 *            the directory to write text file to
	 * @throws IOException
	 *             in case of failure to write file
	 */
	public void writeParticlePositionsToFile(String dirName) throws IOException {
		biomassContainer.bacteriaToFile(dirName + "/iteration"
				+ getFormatedIterationNumber() + ".txt");
	}

	/**
	 * @return the present time
	 */
	public float getTime() {
		return _time;
	}

	/**
	 * @return the present iteration number
	 */
	public int getIterationCounter() {
		return _iterationCounter;
	}

	/**
	 * @return the present iteration number formated to a String
	 */
	public String getFormatedIterationNumber() {
		return getFormatedNumber(_iterationCounter);
	}

	/**
	 * @param n
	 *            number to return in formatted string
	 * @return a number formated to a String
	 */
	public String getFormatedNumber(int n) {
		NumberFormat numberFormat = new DecimalFormat("########");
		numberFormat.setMinimumIntegerDigits(8);
		return numberFormat.format(n);
	}

	/**
	 * Get the current number of particles of a given species
	 * 
	 * @param s
	 * @return number of particles of a given species
	 */
	public int getNumberOfParticles(BiomassSpecies s) {
		return biomassContainer.getNumberOfParticles(s);
	}

	/**
	 * Return a tab separated string with some parameters concerning the present
	 * iteration: iteration number, time
	 * 
	 * @return string with some parameters
	 */
	public String getIterationParameters() {
		long t = (System.currentTimeMillis() - _realInitialTime);
		// The next line is to include real time of simulations
		String s = _iterationCounter + "\t" + _time + "\t" + t;
		return s;
	}

	/**
	 * @return the header for a table with iteration parameters
	 */
	public String getIterationParametersHeader() {
		// The next line is to include real time of simulations
		return "Iteration\tTime [h]\tCurrent real time [ms]";
	}

	/**
	 * @return Returns the dimensionality of the system (2 for 2D, 3 for 3D)
	 */
	public int getDimensionality() {
		return _dimensionality;
	}

	/**
	 * @return Returns the carrierArea.
	 */
	public float getComputationalVolumeMultiplier() {
		return _carrierArea / getComputationalVolumeCarrierArea();
	}

	/**
	 * Set the ratio between carrier area in the reactor and the computational
	 * volume carrier area
	 * 
	 * @param m
	 *            The computational volume multiplier.
	 */
	public void setCarrierArea(float m) {
		this._carrierArea = m;
	}

	/**
	 * @return Returns the residenceTime.
	 */
	public final float getResidenceTime() {
		return _residenceTime;
	}

	/**
	 * @param residenceTime
	 *            The residenceTime to set.
	 */
	public void setResidenceTime(float residenceTime) {
		this._residenceTime = residenceTime;
	}

	/**
	 * @return Returns the _reactorVolume.
	 */
	public float getReactorVolume() {
		return _reactorVolume;
	}

	/**
	 * @param volume
	 *            The _reactorVolume to set.
	 */
	public void setReactorVolume(float volume) {
		_reactorVolume = volume;
	}

	/**
	 * @return Returns the time series of simulation.
	 */
	public UnlimitedFloatArray getTimeSeries() {
		return _timeSeries;
	}

	/**
	 * @return Returns an unmodifiable collection with solute species.
	 */
	public Collection<SoluteSpecies> getSoluteSpecies() {
		return Collections
				.unmodifiableCollection(Arrays.asList(_soluteSpecies));
	}

	/**
	 * @return the number of solute species in the system
	 */
	public int getNumberOfSolutes() {
		return _soluteSpecies.length;
	}

	/**
	 * Solute species at index i in the _soluteSpecies array
	 * 
	 * @param i
	 * @return the solute species number i
	 */
	public SoluteSpecies getSoluteSpecies(int i) {
		return _soluteSpecies[i];
	}

	/**
	 * @return Returns an unmodifiable collection with particulate species.
	 */
	public Collection<ParticulateSpecies> getParticulateSpecies() {
		return Collections.unmodifiableCollection(Arrays
				.asList(_particulateSpecies));
	}

	/**
	 * @return Returns an unmodifiable collection with biomass species.
	 */
	public Collection<BiomassSpecies> getBiomassSpecies() {
		return Collections.unmodifiableCollection(Arrays
				.asList(_biomassSpecies));
	}

	/**
	 * @return Returns the _detachedBiomass.
	 */
	public DetachedBiomassContainer detachedBiomassContainer() {
		return _model._detachedBiomass;
	}

	/**
	 * @return Returns the biomassContainer.
	 */
	public BiomassContainer getBiomassContainer() {
		return biomassContainer;
	}

	/**
	 * @return true if (_time >= _finishIterationTime) or _finishSimulation
	 */
	public boolean endSimulation() {
		// if model is in granular geometry check if maximum diameter
		// was reached
		if (biomassContainer instanceof GranuleBiomassParticleContainer) {
			float diameter = ((GranuleBiomassParticleContainer) biomassContainer)
					.getMaximumRadiusOfGranule() * 2.0f;
			// TODO delete trace
			System.out.println("diameter " + diameter + " of " + _maxRunLength);
			if (diameter >= _maxRunLength)
				return true;
		}
		// otherwise, check other criteria
		// return (_time >= _finishIterationTime) | _finishSimulation
		// | (getCurrentBiovolume() >= _maxBiovolume);
		return (_iterationCounter >= _finishIterationCount || _time >= _finishIterationTime)
				| _finishSimulation
				| (getCurrentBiovolume() >= _maxBiovolume)
				| (getCurrentBiofilmHeight() > _finalBiofilmHeight);
	}

	/**
	 * Update values of maximum specific growth rates for each fixed species
	 */
	public void updateValuesOfCurrentMaximumSpecificGrowthRates() {
		for (int i = 0; i < _particulateSpecies.length; i++) {
			_particulateSpecies[i].computeMaximumSpecificGrowthRate();
		}
	}

	/**
	 * @return the area represented in the computational volume [um^2}
	 */
	public float getComputationalVolumeCarrierArea() {
		return (_dimensionality == 2 ? systemSize.y * _systemGridSide
				: systemSize.y * systemSize.z);
	}

	/**
	 * @return Returns the the depth of the system for 2D computations.
	 */
	public float get2DSystem3rdDimension() {
		if (_dimensionality == 3)
			throw new ModelRuntimeException(
					"Trying to get System depth in a 3D system");
		return _systemGridSide;
	}

	/**
	 * @return the side of a grid element
	 */
	public float getGridElementSide() {
		return _systemGridSide;
	}

	/**
	 * Set a time criterium for finishing simulations. Default is
	 * Float.POSITIVE_INFINITY, meaning that simulations must be stopped
	 * manually
	 * 
	 * @param iterationTime
	 *            The _finishIterationTime to set.
	 */
	public void setFinishIterationTime(float iterationTime) {
		_finishIterationTime = iterationTime;
	}

	/**
	 * Alternatively set a maximum number of iterations before finish
	 * 
	 * @param iterationCount
	 *            sets _finishIterationCount.
	 */

	public void setFinishIterationCount(float iterationCount) {
		_finishIterationCount = iterationCount;
	}

	/**
	 * Set the iteration time step to a fixed value, overriding the models
	 * variable time step feature
	 * 
	 * @param t
	 *            time step
	 */
	public void overrideTimeStep(float t) {
		_timeStep = t;
		_overrideTimeStep = true;
		_writeTimedWriters = true;
	}

	/**
	 * Turn off the overriding of time step
	 */
	public void setTimeStepBackToManual() {
		_timeStep = 0;
		_overrideTimeStep = false;
		_writeTimedWriters = true;
	}

	/**
	 * Set a value for the compulsory time step. This means that at each time
	 * interval defined by attribute _compulsoryTimeStep an iteration will be
	 * performed, inspite of the variable time step. This feature may be used so
	 * that output is presented at regular time intervals
	 * 
	 * @param f
	 */
	public void setCompulsoryTimeStep(float f) {
		_compulsoryTimeStep = f;
		_overrideTimeStep = false;
	}

	/**
	 * Sets a maximum time step, writing outputs at each iteration independently
	 * of the time step value
	 * 
	 * @param f
	 */
	public void setMaximumTimeStep(float f) {
		_maximumTimeStep = f;
	}

	/**
	 * @return Returns true if timed state writers should be written at this
	 *         iteration.
	 */
	public boolean writeTimedWriters() {
		return _writeTimedWriters;
	}

	/**
	 * set the simulation to finish
	 */
	public void setFinishSimulation() {
		_finishSimulation = true;
	}

	/**
	 * @return Returns the _boundaryLayer.
	 */
	public BoundaryLayer getBoundaryLayer() {
		return _boundaryLayer;
	}

	/**
	 * @param t
	 *            The value for the minimum time step in the model simulation.
	 */
	public void setMinimumTimeStep(float t) {
		_minimumTimeStep = t;
	}

	/**
	 * @return the last time constraint or a default time contraint
	 */
	public TimeStepConstraint getLastTimeConstraint() {
		int l = _timeStepConstraints.size();
		if (l > 0) {
			TimeStepConstraint c = (_timeStepConstraints.get(l - 1));
			return c;
		}
		// if iterations have not started yet, return a new TimeConstraint
		return new TimeStepConstraint(0, "Initial time point");
	}

	/**
	 * @return the Singleton instance of model for dessirialization
	 */
	protected Object readResolve() {
		_model = this;
		return _model;
	}

	private BiomassBooleanMatrix _biomassBooleanMatrix;

	/**
	 * Returns a boolean matrix with the locations of the biomass
	 * 
	 * @return the present instance of BiomassBooleanMatrix
	 */
	public BiomassBooleanMatrix getBiomassBooleanMatrix() {
		if (_biomassBooleanMatrix == null)
			_biomassBooleanMatrix = new BiomassBooleanMatrix();
		_biomassBooleanMatrix.upateBiomassConcentrations(_particulateSpecies);
		return _biomassBooleanMatrix;
	}

	/**
	 * This method is for debugging purposes. It will only work in the case the
	 * model used is a 2D tube reactor
	 * 
	 * @param dirName
	 *            the directory to write out to
	 */
	public void writeTubeDetachmentAlgorithmOutput(String dirName) {
		// create directory for this iteration
		File dir = new java.io.File(dirName + "/iteration"
				+ getFormatedIterationNumber());
		try {
			dir.mkdir();
		} catch (SecurityException e) {
			throw new SystemEditViolationException("directory " + dir
					+ " could not be created");
		}
		((TubeBiomassParticleContainer) biomassContainer)
				.writeDetachmentInfo(dir.toString() + "/");
	}

	/**
	 * Set a maximum biovolume as a criterium to end simulations
	 * 
	 * @param biovolume
	 */
	public void setMaxBiovolume(float biovolume) {
		_maxBiovolume = biovolume;
	}

	/**
	 * Set a maximum runlength as a criterium to end simulations. Only works for
	 * granule
	 * 
	 * @param r
	 *            maximum run length
	 */
	public void setMaxRunLength(float r) {
		_maxRunLength = r;
	}

	/**
	 * Set maximum biofilm height as criterium to end simulaiton.
	 * 
	 * @param maxHeight
	 */
	public void setMaxHeight(float maxHeight) {
		_finalBiofilmHeight = maxHeight;
	}
}