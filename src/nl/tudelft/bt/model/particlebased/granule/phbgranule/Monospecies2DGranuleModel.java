/*
 * Created on Sep 8, 2004
 */
package nl.tudelft.bt.model.particlebased.granule.phbgranule;

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
	private String outputDirectory = "E:/results/SBR/lixo/";

	// VARIABLES
	// attibutes that are variable
	private float _qMax;

	private float _bDecay = 3.3e-3f;

	private float _cO;

	private float _cS = 0.6f; // gCOD-S/L

	private float _kd;

	private float _boundaryLayerThickness;

	private int _initialParticleNumber;

	// PAO'S kinetics
	protected static float QSMAX_PAO = 0.4167f; // [gCOD-S/gCOD-PAO/h]

	protected static float MPAOANAEROBIC = 0.0021f; // [gP/gCOD-PAO/h]

	protected static float QPHBPAO = 0.3146f; // [gCOD-PHB/gCOD-PAO/h]

	protected static float QPPPAO = 0.0187f; // [gP/gCOD-PAO/h]

	protected static float QGLYPAO = 0.0454f; // [gCOD-GLY/gCOD-PAO/h]

	// CONSTANTS
	// geometry (2D or 3D)
	protected static int GEOMETRY = 2;

	// Oxygen
	private static float DIFFUSIVITY_OXYGEN = 8e6f; // [um2/h]

	// Substrate
	private static float DIFFUSIVITY_SUBSTRATE = 4e6f; // [um2/h]

	// Particulate species
	protected static float DENSITY_ACTIVE = 200f; // [gCOD-H/L]

	protected static float DENSITY_INERT = DENSITY_ACTIVE; // [gCOD-I/L]

	protected static float DENSITY_PHB = 1000; // [gCOD-PHB/L]

	// Yield coefficients
	private static float YHS = 0.495f; // [gCOD-H/gCOD-S]

	private static float YPS = 0.667f; // [gCOD-P/gCOD-S]

	private static float YHP = 0.668f; // [gCOD-H/gCOD-P]

	private static float KO = 3.5e-4f; // [gO/l]

	private static float KO_star = 7E-04f; // [gO/l]

	private static float KS = 4e-3f; // [gCOD-S/l]

	private static float KP = 1; // [gCOD-P/gCOD-H]

	// PAO's properties
	private static float KS_PAO_ANAEROBIC = 4e-3f; // [gCOD-S/L]

	private static float KGLY_PAO_ANAEROBIC = 0.01e-3f; // [gCOD-GLY/L]

	private static float KPP_PAO_ANAEROBIC = 0.01e-3f; // [gP/L]

	private static float KO_PAO = 2e-4f; // [gO/l]

	private static float KPHP_FRACTION_PAO = 0.33f; // [gCOD-PHB/gCOD-PAO]

	private static float KPP_PAO_AEROBIC = 0.35f; // [gP/gCOD-PAO]

	private static float KGLY_PAO_AEROBIC = 0.45f; // [gCOD-GLY/gCOD-PAO]

	// PHB consumption
	private static float KPHB = 0.15f; // [gCOD_PHB/gCOD_X/h]

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
		_qMax = _gui.getUserQMax();
		_bDecay = _gui.getUserBDecay(); // [gO/l]
		_initialParticleNumber = _gui.getUserInoculum(); // [particles]
		_cO = _gui.getUserCO(); // [gO/l]
		_boundaryLayerThickness = _gui.getUserBoundaryLayerThickness(); // [um]
		_kd = _gui.getUserDetachment(); // 
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
		_qMax = _gui.getUserQMax(); // [h-1]
		_bDecay = _gui.getUserBDecay(); // [h-1]
		_initialParticleNumber = _gui.getUserInoculum(); // [particles]
		_cO = _gui.getUserCO(); // [gO/l]
		_boundaryLayerThickness = _gui.getUserBoundaryLayerThickness(); // [um]
		_kd = _gui.getUserDetachment(); // 
		// update the growth rate
		_substrateUptake.setConstant(_qMax);
		_growthH.setConstant(_qMax * YHS);
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
		setSystemSpaceParameters(GEOMETRY, SYSTEMSIZE, RELATIVEMAXIMUMRADIUS,
				RELATIVEMINIMUMRADIUS, _boundaryLayerThickness / SYSTEMSIZE,
				GRIDSIZE, KSHOV);
		// initialize
		initializeSystemSpace();
		// add writers here
		intializeStateWriters(outputDirectory);
		addTimedStateWriter(new SoluteConcentrationWriter());
		addTimedStateWriter(new SolidsConcentrationWriter());
		addTimedStateWriter(new ParticlePositionWriter());
		// the simulation parameters writter
		SimulationResultsWriter spw = new SimulationResultsWriter();
		// add the state writer to the appllication
		addStateWriter(spw);
		// add series to the state writer
		spw.addSeries(Model.model().detachedBiomassContainer()
				.getTotalDetachedBiomassSeries());
		spw
				.addSeries(Model.model().detachedBiomassContainer()
						.getErodedBiomassSeries());
		spw.addSeries(Model.model().detachedBiomassContainer()
				.getSloughedBiomassSeries());
		// create the solute species and particulate species
		initializeDiffusionReactionSystem();
		// create the detachment function
		initializeDetachmentFunction();
		// add bulk concentrations of all solutes as variable series
		for (Iterator i = Model.model().getSoluteSpecies().iterator(); i
				.hasNext();) {
			SoluteSpecies s = (SoluteSpecies) i.next();
			spw.addSeries(s.getBulkConcentrationSeries());
			spw.addSeries(s.getRateTimeSeries());
		}
		// add particulate global masses to write
		for (Iterator i = Model.model().getParticulateSpecies().iterator(); i
				.hasNext();) {
			ParticulateSpecies s = (ParticulateSpecies) i.next();
			spw.addSeries(s.getTotalMassSeries());
		}
	}

	/**
	 * Initialize the species and reactions
	 */
	public void initializeDiffusionReactionSystem() throws ModelException {
		// oxygen
		SoluteSpecies oxygen;
		// substrate
		SoluteSpecies substrate;
		// create the particulate species
		// active mass
		ParticulateSpecies activeH;
		// inert mass
		ParticulateSpecies inert;
		// PHB
		ParticulateSpecies phb;
		// create the particulate species used for PAO
		// active mass
		ParticulateSpecies activePao;
		// PHB
		ParticulateSpecies phbPao;
		// glycogen
		ParticulateSpecies glycogen;
		// poly-p
		ParticulateSpecies polyP;
		try {
			oxygen = new SoluteSpecies("oxygen", DIFFUSIVITY_OXYGEN);
			_oxygenBulkConcentration = new ConstantBulkConcentration(_cO);
			oxygen.setBulkConcentration(_oxygenBulkConcentration);
			substrate = new SoluteSpecies("substrate", DIFFUSIVITY_SUBSTRATE);
			_substrateBulkConcentration = new ConstantBulkConcentration(_cS);
			substrate.setBulkConcentration(_substrateBulkConcentration);
			activeH = new ParticulateSpecies("activeMass", DENSITY_ACTIVE,
					Color.blue);
			inert = new ParticulateSpecies("inert", DENSITY_INERT, Color.white);
			phb = new ParticulateSpecies("phb", DENSITY_PHB, Color.yellow);
			activePao = new ParticulateSpecies("activePAO", DENSITY_ACTIVE,
					Color.red);
			phbPao = new ParticulateSpecies("PHB PAO", DENSITY_PHB,
					Color.yellow);
			glycogen = new ParticulateSpecies("glycogen", DENSITY_PHB,
					Color.green);
			polyP = new ParticulateSpecies("polyP", DENSITY_PHB, Color.magenta);
		} catch (MultigridSystemNotSetException e) {
			throw new ModelRuntimeException(e.toString());
		}
		// array of fixed species that constitute the heterotrophic species
		ParticulateSpecies[] spH = { activeH, inert, phb };
		float[] fractionalVolumeCompositionH = { 1.0f, 0, 0 };
		// create the biomass species
		BiomassSpecies heterotroph = new BiomassSpecies("heterotroph", spH,
				fractionalVolumeCompositionH);
		heterotroph.setActiveMass(activeH);
		heterotroph.setInertMass(inert);
		// array of fixed species that constitute H1
		ParticulateSpecies[] spPao = { activePao, phbPao, glycogen, polyP };
		float[] fractionalVolumeCompositionPao = { 1.0f, 0, 0, 0 };
		// create the biomass species
		BiomassSpecies pao = new BiomassSpecies("pao", spPao,
				fractionalVolumeCompositionPao);
		pao.setActiveMass(activePao);
		// Create the Reaction factors, Monod and inhibition coefficients
		// for heterotrophs
		ProcessFactor mOH = new Saturation(oxygen, KO);
		ProcessFactor mOH_star = new Saturation(oxygen, KO_star);
		ProcessFactor mS = new Saturation(substrate, KS);
		ProcessFactor iS = new Inhibition(substrate, KS);
		ProcessFactor iPhb = new InhibitionFromFraction(phb, activeH, KP);
		// for PAO's
		// anaerobic
		ProcessFactor mSPao = new Saturation(substrate, KS_PAO_ANAEROBIC);
		ProcessFactor mGlyPao = new Saturation(glycogen, KGLY_PAO_ANAEROBIC);
		ProcessFactor mPolyPPao = new Saturation(polyP, KPP_PAO_ANAEROBIC);
		ProcessFactor iOPao = new Inhibition(oxygen, KO_PAO);
		// aerobic
		ProcessFactor mOPao = new Saturation(oxygen, KO_PAO);
		ProcessFactor mPhbPaoFraction = new SaturationFromFraction(phbPao,
				activePao, KPHP_FRACTION_PAO);
		ProcessFactor iPpPao = new InhibitionFromFraction(polyP, activePao,
				KPP_PAO_AEROBIC);
		ProcessFactor iGlyPao = new InhibitionFromFraction(glycogen, activePao,
				KGLY_PAO_AEROBIC);
		// create the reactions
		// Heterotrophic:
		// substrate uptake
		_substrateUptake = new Reaction("substrateUptake", activeH, _qMax, 3);
		_substrateUptake.addFactor(mOH);
		_substrateUptake.addFactor(mS);
		_substrateUptake.addFactor(iPhb);
		// biomass growth
		_growthH = new Reaction("growth H", activeH, _qMax * YHS, 2);
		_growthH.addFactor(mS);
		_growthH.addFactor(mOH_star);
		// PHB consumption
		Reaction phbConsumption = new Reaction("PHB consumption", phb, KPHB, 2);
		phbConsumption.addFactor(iS);
		phbConsumption.addFactor(mOH);
		// heterotroph decay
		_heterotrophDecay = new Reaction("heterotrophDecay", activeH, _bDecay,
				0);
		// PAO'S
		// storage of PHB:
		Reaction paoPhbStorage = new Reaction("paoPhbStorage", activePao,
				QSMAX_PAO, 3);
		paoPhbStorage.addFactor(mSPao);
		paoPhbStorage.addFactor(mGlyPao);
		paoPhbStorage.addFactor(mPolyPPao);
		// anaerobic maintenance
		Reaction paoAnaerobicMaintenance = new Reaction(
				"paoAnaerobicMaintenance", activePao, MPAOANAEROBIC, 2);
		paoAnaerobicMaintenance.addFactor(iOPao);
		paoAnaerobicMaintenance.addFactor(mPolyPPao);
		// aerobic consumption of PHB
		Reaction paoPhbConsumption = new Reaction("paoPhbConsumption",
				activePao, QPHBPAO, 2);
		paoPhbConsumption.addFactor(mPhbPaoFraction);
		paoPhbConsumption.addFactor(mOPao);
		// aerobic storage of Poly-P
		Reaction paoPpStorage = new Reaction("paoPpStorage", activePao, QPPPAO,
				2);
		paoPpStorage.addFactor(iPpPao);
		paoPpStorage.addFactor(mOPao);
		// aerobic storage of glycogen
		Reaction paoGlyStorage = new Reaction("paoGlyStorage", activePao,
				QGLYPAO, 3);
		paoGlyStorage.addFactor(iGlyPao);
		paoGlyStorage.addFactor(mPhbPaoFraction);
		paoGlyStorage.addFactor(mOPao);

		// Anaerobic processes:

		// assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsHactive = new NetReaction(3);
		rsHactive.addReaction(_growthH, 1);
		rsHactive.addReaction(phbConsumption, YHP);
		rsHactive.addReaction(_heterotrophDecay, -1);
		activeH.setProcesses(rsHactive);
		// phb
		NetReaction rsPhb = new NetReaction(3);
		rsPhb.addReaction(_substrateUptake, YPS);
		rsPhb.addReaction(_growthH, -YPS / YHS);
		rsPhb.addReaction(phbConsumption, -1);
		phb.setProcesses(rsPhb);
		// inert mass
		NetReaction rsInert = new NetReaction(1);
		rsInert.addReaction(_heterotrophDecay, 0.4f);
		inert.setProcesses(rsInert);
		// PAO active mass
		NetReaction rsActivePao = new NetReaction(3);
		rsActivePao.addReaction(paoPhbConsumption, 0.667f);
		rsActivePao.addReaction(paoPpStorage, -0.229f);
		rsActivePao.addReaction(paoGlyStorage, -0.9f);
		activePao.setProcesses(rsActivePao);
		// PAO phb
		NetReaction rsPhbPao = new NetReaction(2);
		rsPhbPao.addReaction(paoPhbStorage, 1.5f);
		rsPhbPao.addReaction(paoPhbConsumption, -1f);
		phbPao.setProcesses(rsPhbPao);
		// PolyP
		NetReaction rsPp = new NetReaction(3);
		rsPp.addReaction(paoPhbStorage, 0.3f);
		rsPp.addReaction(paoAnaerobicMaintenance, -1f);
		rsPp.addReaction(paoPpStorage, 1);
		polyP.setProcesses(rsPp);
		// Glycogen
		NetReaction rsGly = new NetReaction(2);
		rsGly.addReaction(paoPhbStorage, -0.5f);
		rsGly.addReaction(paoGlyStorage, 1);
		glycogen.setProcesses(rsGly);
		// assign reaction stoichiometry to the solutes
		// oxygen
		NetReaction rsOxygen = new NetReaction(5);
		rsOxygen.addReaction(_substrateUptake, -(1 - YPS));
		rsOxygen.addReaction(_growthH, -(YPS / YHS - 1));
		rsOxygen.addReaction(paoPhbConsumption, -0.333f);
		rsOxygen.addReaction(paoPpStorage, -0.229f);
		rsOxygen.addReaction(paoGlyStorage, 0.1f);
		oxygen.setProcesses(rsOxygen);
		// substrate
		NetReaction rsSubstrate = new NetReaction(2);
		rsSubstrate.addReaction(_substrateUptake, -1);
		rsSubstrate.addReaction(paoPhbStorage, -1);
		substrate.setProcesses(rsSubstrate);
		// add the species to system
		addBiomassSpecies(heterotroph);
		addBiomassSpecies(pao);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(substrate);
		//
		super.initializeDiffusionReactionSystem();
	}

	/**
	 * Inocultes the system, this is called in
	 * super.initializeDiffusionReactionSystem which in turn is called in
	 * this.initializeDiffusionReactionSystem
	 */
	protected void inoculate() {
		int[] nCells = { _initialParticleNumber, _initialParticleNumber };
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