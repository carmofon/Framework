package nl.tudelft.bt.model.examples.monospecies2D;

import java.awt.Color;
import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.bulkconcentrations.BulkConcentration;
import nl.tudelft.bt.model.bulkconcentrations.ConstantBulkConcentration;
import nl.tudelft.bt.model.detachment.*;
import nl.tudelft.bt.model.detachment.levelset.functions.Height2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.util.Semaphore;

/**
 * Handler for a monospcies model with a single solute species
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Monospecies2DModel extends ModelHandler {
	// VARIABLES
	// attibutes that are variable
	private float _qMax;

	private float _kO;

	private float _cO;

	private float _kd;

	private float _boundaryLayerThickness;

	private int _initialParticleNumber;

	// CONSTANTS
	// geometry (2D or 3D)
	protected static int GEOMETRY = 2;

	// Oxygen
	private static float DIFFUSIVITY = 8e6f; // [um2/h]

	// Particulate species
	protected static float DENSITY = 200f; // [gCOD_X/l]

	// Yield coefficients
	private static float YSX = 0.495f; // [gCOD_X/gCOS_S]

	// Computation parameters
	protected static float SYSTEMSIZE = 2000; // [um]

	protected static float RELATIVEMAXIMUMRADIUS = 0.006f;

	protected static float RELATIVEMINIMUMRADIUS = RELATIVEMAXIMUMRADIUS * 0.0001f;

	// multigrid settings
	protected static int NPRE = 5;

	protected static int NPOS = 50;

	// other model parameters
	protected static int GRIDSIZE = 33; // multigrid grid side

	protected static float KSHOV = 1.0f; // shoving parameter[dim/less]

	//
	// references to components with variable properties
	private Reaction _substrateUptake; // Substrate uptake rate reaction

	private Saturation _mO; // Mono saturetion factor for oxygen

	private BulkConcentration _oxygenBulkConcentration; // oxygen solute species

	private Height2MassDetachment _detachmentFunction; // detachment speed

	// function

	// the Semaphores
	private Semaphore _endSemaphore;

	private Semaphore _pauseSemaphore;

	// reference to the applet to get the user input
	private Monospecies2DGUI _gui;

	/**
	 * 
	 */
	public Monospecies2DModel(Monospecies2DGUI a) throws ModelException {
		_gui = a;
		// set initial values for parameters
		_qMax = _gui.getUserQMax();
		_kO = _gui.getUserKO(); // [gO/l]
		_initialParticleNumber = _gui.getUserInoculum(); // [particles]
		_cO = _gui.getUserCO(); // [gO/l]
		_boundaryLayerThickness = _gui.getUserBoundaryLayerThickness(); // [um]
		_kd = _gui.getUserDetachment();
		// intialize the semaphores
		_endSemaphore = new Semaphore();
		_pauseSemaphore = new Semaphore();
		// set pre and pos steps for multigrid
		MultigridVariable.setSteps(NPRE, NPOS);
		// intialize system
		initializeSystem();
		updateParameters();
	}

	/**
	 * Initializes simulation
	 * 
	 * @throws Exception
	 */
	public void start() throws Exception {
		// start a new thread
		_endSemaphore.close();
		while (!Model.model().endSimulation()) {
			_pauseSemaphore.waitIfClosed();
			performFullIteration();
			updateParameters();
			_gui.updateParticleVisualizer();
		}
		// simulation is finished, release semafore
		System.out.println("done running, open end semaphore");
		_endSemaphore.open();
		System.out.println("opened");
	}

	/**
	 * pauses simulation if it is running or continues simulation if it is
	 * paused
	 */
	public void pause() {
		_pauseSemaphore.switchOpenClose();
	}

	/**
	 * Restarts the model
	 * 
	 * @throws ModelException
	 */
	public void restart() throws ModelException {
		setFinishSimulationSwitch();
		_pauseSemaphore.open();
		_endSemaphore.waitIfClosed();
		updateParameters();
		initializeSystem();
	}

	/**
	 * Updates the runtime variables
	 */
	private void updateParameters() {
		_qMax = _gui.getUserQMax();
		_kO = _gui.getUserKO();
		_initialParticleNumber = _gui.getUserInoculum();
		_cO = _gui.getUserCO();
		_boundaryLayerThickness = _gui.getUserBoundaryLayerThickness();
		_kd = _gui.getUserDetachment();
		_substrateUptake.setConstant(_qMax);
		_mO.setK(_kO);
		_oxygenBulkConcentration.setValue(_cO);
		_detachmentFunction.setDetachmentRateConstant(_kd);
		_boundaryLayer.setThickness(_boundaryLayerThickness);
	}

	/**
	 * Resets and then initializes the whole system, including the species,
	 * reactions and inoculation
	 * 
	 * @throws ModelException
	 */
	private void initializeSystem() throws ModelException {
		// reset the model properties;
		Model.model().setSeed((int) System.currentTimeMillis());
		reset();
		// create the system
		setSystemSpaceParameters(GEOMETRY, SYSTEMSIZE, RELATIVEMAXIMUMRADIUS,
				RELATIVEMINIMUMRADIUS, _boundaryLayerThickness / SYSTEMSIZE,
				GRIDSIZE, KSHOV);
		// initialize
		initializeSystemSpace();
		// create the solute species and particulate species
		initializeDiffusionReactionSystem();
		initializeDetachmentFunction();
	}

	/**
	 * Initialize the species and reactions
	 */
	public void initializeDiffusionReactionSystem() throws ModelException {
		// oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", DIFFUSIVITY);
		_oxygenBulkConcentration = new ConstantBulkConcentration(_cO);
		oxygen.setBulkConcentration(_oxygenBulkConcentration);
		// create the particulate species
		// active mass
		ParticulateSpecies activeH1 = new ParticulateSpecies("activeMass",
				DENSITY, Color.ORANGE);
		// array of fixed species that constitute H1
		ParticulateSpecies[] spH1 = { activeH1 };
		float[] fractionalVolumeCompositionH1 = { 1.0f };
		// create the biomass species
		BiomassSpecies speciesH1 = new BiomassSpecies("heterotroph", spH1,
				fractionalVolumeCompositionH1);
		speciesH1.setActiveMass(activeH1);
		// Create the Reaction factors, Monod and inhibition coefficients
		_mO = new Saturation(oxygen, _kO);
		// create the reactions
		// substrate uptake
		_substrateUptake = new Reaction("substrateUptake", activeH1, _qMax, 1);
		_substrateUptake.addFactor(_mO);
		// assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsH1active = new NetReaction(1);
		rsH1active.addReaction(_substrateUptake, YSX);
		activeH1.setProcesses(rsH1active);
		// assign reaction stoichiometry to the solutes
		// oxygen
		NetReaction rsOxygen = new NetReaction(1);
		rsOxygen.addReaction(_substrateUptake, -(1 - YSX));
		oxygen.setProcesses(rsOxygen);
		// add the species to system
		addBiomassSpecies(speciesH1);
		addSoluteSpecies(oxygen);
		//
		super.initializeDiffusionReactionSystem();
	}

	/**
	 * Inocultes the system, this is called in
	 * super.initializeDiffusionReactionSystem which in turn is called in
	 * this.initializeDiffusionReactionSystem
	 */
	protected void inoculate() {
		int[] nCells = { _initialParticleNumber };
		inoculateRandomly(nCells);
	}

	/**
	 * Initializes the detachment function
	 */
	public void initializeDetachmentFunction() {
		_detachmentFunction = new Height2MassDetachment(_kd);
		setDetachmentHandler(_detachmentFunction);
	}

	/**
	 * @param layerThickness
	 *            The _boundaryLayerThickness to set.
	 */
	public void setBoundaryLayerThickness(float layerThickness) {
		_boundaryLayerThickness = layerThickness;
	}

	/**
	 * @param _co
	 *            The _cO to set.
	 */
	public void setOxygenBulkConcentration(float _co) {
		_cO = _co;
	}

	/**
	 * @param _kd
	 *            The _kd to set.
	 */
	public void setDetachmentRateConstant(float _kd) {
		this._kd = _kd;
	}

	/**
	 * @param _ko
	 *            The _kO to set.
	 */
	public void setOxygenMonodConstant(float _ko) {
		_kO = _ko;
	}

	/**
	 * @param max
	 *            The _qMax to set.
	 */
	public void setSubstrateUptakeRate(float max) {
		_qMax = max;
	}
}