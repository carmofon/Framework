/*
 * Created on Sep 8, 2004
 */
package nl.tudelft.bt.model.particlebased.granule.model2d;

import java.awt.Color;
import java.util.Iterator;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.output.ParticlePositionWriter;
import nl.tudelft.bt.model.apps.output.SimulationResultsWriter;
import nl.tudelft.bt.model.apps.output.SolidsConcentrationWriter;
import nl.tudelft.bt.model.apps.output.SoluteConcentrationWriter;
import nl.tudelft.bt.model.bulkconcentrations.BulkConcentration;
import nl.tudelft.bt.model.bulkconcentrations.ConstantBulkConcentration;
import nl.tudelft.bt.model.detachment.*;
import nl.tudelft.bt.model.detachment.levelset.functions.Radius2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.multigrid.boundary_layers.SphericalDilationBoundaryLayer;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.util.Semaphore;

/**
 * Handler for a monospcies model with a single solute species
 * 
 * @author jxavier
 */
public class Monospecies2DGranuleModel extends GranuleModelHandler {
	// Set the output directory
	private String outputDirectory = "/Users/jxavier/model/lixo";

	// VARIABLES
	// attibutes that are variable
	private float _uMax;

	private float _bDecay;

	private float _cO;

	private float _kd;

	private float _boundaryLayerThickness;

	private int _initialParticleNumber;

	// CONSTANTS
	// geometry (2D or 3D)
	protected static int GEOMETRY = 2;

	// Oxygen
	private static float DIFFUSIVITY_OXYGEN = 8.3e6f; // [um2/h]

	// Particulate species
	protected static float DENSITY_ACTIVE = 200f / 0.7f; // [gCOD-H/L]

	protected static float DENSITY_INERT = DENSITY_ACTIVE; // [gCOD-I/L]

	// Yield coefficients
	private static float YOH = 0.505f; // [gCOD-H/gCOD-S]

	private static float KO = 3.5e-4f; // [gO/l]

	// Computation parameters
	protected static float SYSTEMSIZE = 3000; // [um]

	protected static float RELMAXRAD = 0.0038f;

	protected static float RELMINRAD = RELMAXRAD * 0.0001f;

	// multigrid settings
	protected static int NPRE = 5;

	protected static int NPOS = 50;

	// other model parameters
	protected static int GRIDSIZE = 33; // multigrid grid side

	protected static float KSHOV = 1.0f; // shoving parameter[dim/less]

	//
	// references to components with variable properties
	private Reaction _growthH;

	private Reaction _heterotrophDecay; // Substrate uptake rate reaction

	private BulkConcentration _oxygenBulkConcentration; // oxygen solute species

	private BulkConcentration _substrateBulkConcentration; // substrate solute

	// species

	private Radius2MassDetachment _detachmentFunction; // detachment speed

	// function

	// the Semaphores
	private Semaphore _endSemaphore;

	private Semaphore _pauseSemaphore;

	// reference to the applet to get the user input
	private Monospecies2DGranuleGUI _gui;

	/**
	 * 
	 */
	public Monospecies2DGranuleModel(Monospecies2DGranuleGUI a)
			throws ModelException {
		_gui = a;
		// set initial values for parameters
		_uMax = _gui.getUserUMax();
		_bDecay = _gui.getUserBDecay(); // [gO/l]
		_initialParticleNumber = _gui.getUserInoculum(); // [particles]
		_cO = _gui.getUserCO(); // [gO/l]
		_boundaryLayerThickness = _gui.getUserBoundaryLayerThickness(); // [um]
		_kd = _gui.getUserDetachment() / DENSITY_ACTIVE; // [1e-15
		// gCOD-H/um^4/h]
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
		updateParameters();
		while (!Model.model().endSimulation()) {
			_pauseSemaphore.waitIfClosed();
			performFullIteration();
			// write the state to file
			writeState();
			//
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
		_uMax = _gui.getUserUMax(); // [h-1]
		_bDecay = _gui.getUserBDecay(); // [h-1]
		_initialParticleNumber = _gui.getUserInoculum(); // [particles]
		_cO = _gui.getUserCO(); // [gO/l]
		_boundaryLayerThickness = _gui.getUserBoundaryLayerThickness(); // [um]
		_kd = _gui.getUserDetachment() / DENSITY_ACTIVE; // [1e-15
		// gCOD-H/um^4/h]
		// update the growth rate
		_growthH.setConstant(_uMax);
		//
		_heterotrophDecay.setConstant(_bDecay);
		_oxygenBulkConcentration.setValue(_cO);
		_detachmentFunction.setDetachmentRateConstant(_kd);
		((SphericalDilationBoundaryLayer) _boundaryLayer)
				.setThickness(_boundaryLayerThickness);
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
		setSystemSpaceParameters(GEOMETRY, SYSTEMSIZE, RELMAXRAD, RELMINRAD,
				_boundaryLayerThickness / SYSTEMSIZE, GRIDSIZE, KSHOV);
		// initialize
		initializeSystemSpace();
		/*
		 * // add writers here intializeStateWriters(outputDirectory);
		 * addTimedStateWriter(new SoluteConcentrationWriter());
		 * addTimedStateWriter(new SolidsConcentrationWriter());
		 * addTimedStateWriter(new ParticlePositionWriter()); // the simulation
		 * parameters writter SimulationParametersWriter spw = new
		 * SimulationParametersWriter(); // add the state writer to the
		 * appllication addStateWriter(spw); // add series to the state writer
		 * spw.addSeries(Model.detachedBiomassContainer()
		 * .getTotalDetachedBiomassSeries()); spw
		 * .addSeries(Model.detachedBiomassContainer()
		 * .getErodedBiomassSeries());
		 * spw.addSeries(Model.detachedBiomassContainer()
		 * .getSloughedBiomassSeries());
		 */
		// create the solute species and particulate species
		initializeDiffusionReactionSystem();
		// create the detachment function*/
		initializeDetachmentFunction();
		/*
		 * // add bulk concentrations of all solutes as variable series for
		 * (Iterator i = Model.model().getSoluteSpecies().iterator(); i
		 * .hasNext();) { SoluteSpecies s = (SoluteSpecies) i.next();
		 * spw.addSeries(s.getBulkConcentrationSeries());
		 * spw.addSeries(s.getRateTimeSeries()); } // add particulate global
		 * masses to write for (Iterator i =
		 * Model.model().getParticulateSpecies().iterator(); i .hasNext();) {
		 * ParticulateSpecies s = (ParticulateSpecies) i.next();
		 * spw.addSeries(s.getTotalMassSeries()); } // write the header in
		 * parameters.txt file spw.initializeParametersWriting();
		 */
	}

	/**
	 * Initialize the species and reactions
	 */
	public void initializeDiffusionReactionSystem() throws ModelException {
		// oxygen
		SoluteSpecies oxygen;
		// create the particulate species
		// active mass
		ParticulateSpecies activeH;
		// inert mass
		ParticulateSpecies inert;
		try {
			oxygen = new SoluteSpecies("oxygen", DIFFUSIVITY_OXYGEN);
			_oxygenBulkConcentration = new ConstantBulkConcentration(_cO);
			oxygen.setBulkConcentration(_oxygenBulkConcentration);
			activeH = new ParticulateSpecies("activeMass", DENSITY_ACTIVE,
					Color.orange);
			inert = new ParticulateSpecies("inert", DENSITY_INERT, Color.white);
			// array of fixed species that constitute the heterotrophic species
			ParticulateSpecies[] spH = { activeH, inert };
			float[] fractionalVolumeCompositionH = { 1.0f, 0 };
			// create the biomass species
			BiomassSpecies heterotroph = new BiomassSpecies("heterotroph", spH,
					fractionalVolumeCompositionH);
			heterotroph.setActiveMass(activeH);
			heterotroph.setInertMass(inert);
			// Create the Reaction factors, Monod and inhibition coefficients
			// for heterotrophs
			ProcessFactor mOH = new Saturation(oxygen, KO);
			// create the reactions
			// biomass growth
			_growthH = new Reaction("growth H", activeH, _uMax, 1);
			_growthH.addFactor(mOH);
			// heterotroph decay
			_heterotrophDecay = new Reaction("heterotrophDecay", activeH,
					_bDecay, 0);
			// assign reaction to the species through ReactionStoichiometries
			// active mass
			NetReaction rsHactive = new NetReaction(2);
			rsHactive.addReaction(_growthH, 1);
			rsHactive.addReaction(_heterotrophDecay, -1);
			activeH.setProcesses(rsHactive);
			// inert mass
			NetReaction rsInert = new NetReaction(1);
			rsInert.addReaction(_heterotrophDecay, 0.4f);
			inert.setProcesses(rsInert);
			// assign reaction stoichiometry to the solutes
			// oxygen
			NetReaction rsOxygen = new NetReaction(5);
			rsOxygen.addReaction(_growthH, -YOH);
			oxygen.setProcesses(rsOxygen);
			// add the species to system
			addBiomassSpecies(heterotroph);
			addSoluteSpecies(oxygen);
			//
			super.initializeDiffusionReactionSystem();
		} catch (MultigridSystemNotSetException e) {
			throw new ModelRuntimeException(e.toString());
		}

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
		_detachmentFunction = new Radius2MassDetachment(_kd);
		setDetachmentHandler(_detachmentFunction);
	}

}