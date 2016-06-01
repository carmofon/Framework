package nl.tudelft.bt.model.work.benzoate;

import java.awt.Color;
import java.util.Iterator;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Height2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;

/**
 * Two species system with Acinetobacter sp. and Pseudomonas putida
 * -- added bulk concentrations of benzyl alcohol and benzoate accumulate
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu)
 */
public class TwoSpecies3 extends ModelHandler {
	// output directory name
	private static String outputDirectory = "C:/results/twospecies3/";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.
	// The output directory is were the program will store all the results.
	// Choose a path to an existing folder in your system.
	// EXAMPLE: if you choose "C:\results\twospecies1" directory "c:\results"
	// must
	// exist in your computer. The directory "C:\results\twospecies1" will be
	// created
	// if it is non-existant. If it exists, its contents will be deleted
	// during the program initialization

	// geometry of simulation (2D or 3D)
	private static int geometry = 3;

	// Parameters for solute species
	// Benzyl alcohol
	private static float benzylOHBulkConcentration = 0.042f; // [gC/L]

	private static float benzylOHDiffusivity = 4e6f; // [um2/h]

	private static float KBenzylOH = 0.010f; // [gO/l]

	// Benzoate
	private static float benzoateBulkConcentration = 0; // [gO/l]

	private static float benzoateDiffusivity = 4e6f; // [um2/h]

	private static float KBenzoate = 0.020f; // [gC/L]

	// Oxygen
	private static float oxygenBulkConcentration = 0.001f; // [gO/L]

	private static float oxygenDiffusivity = 8e6f; // [um2/h]

	private static float KOxygen = 0.004e-3f; // [gO/l]

	// Parameters for particulate species
	// density of active biomass
	private static float densityX = 130f; // [gC/L]

	// Yield coefficients
	// biomass on substrate
	private static float YXBenzylOH = 0.03f; // [gC/gC]

	private static float YXBenzoate = 0.8f; // [gCOD_X/gCOS_S]

	private static float YBenzoateProduction = 0.4f; // [gCOD_X/gCOS_S]

	// Processes
	// Benzyl alcohol uptake by Acinetobacter
	private static float qBenzylOHMax = 2f; // [gC/gC/h]

	// Benzoate uptake by P. putida
	private static float qBenzoate = 20f; // [gC/gC/h]

	// Computation parameters
	private static float systemSize = 200f; // [micron]

	// maximum radius of biomass particles, relative to the system size
	//private static float relativeMaximumRadius = 3f / systemSize;
	private static float relativeMaximumRadius = 0.5f / systemSize;

	// minimum radius of biomass particles, relative to the system size
	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// boundary layer thickness, relative to the system size
	private static float relativeBoundaryLayer = 0.1f;

	// other model parameters
	private static int gridSide = 33; // multigrid grid side

	private static float kShov = 1.0f; // shoving parameter

	// detachment rate coefficient
	private static float rdetach = 0;// detachment constant

	// private static float rdetach = 3e-3f;// detachment constant

	//reactor properties
	private static float reactorVolume = 0.160e-3f; // [L]
	private static float residenceTime = 0.053f; // [h]
	private static float carrierArea = 0.04f; // [dm2]

	private static float finishSimulationTime = 240f;
	
	// inoculation
	private static int initialCellNumber = 30;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// ---- Define the solutes ----
		// --- Benzyl alcohol
		SoluteSpecies benzylOH = new SoluteSpecies("benzylOH",
				benzylOHDiffusivity);
		benzylOH.setBulkConcentration(new ConstantBulkConcentration(
				benzylOHBulkConcentration));
		//benzylOH.setBulkConcentration(new DynamicBulkConcentrationImplicit(
		//		benzylOHBulkConcentration));
		// --- Benzoate
		SoluteSpecies benzoate = new SoluteSpecies("benzoate",
				benzoateDiffusivity);
		benzoate.setBulkConcentration(new ConstantBulkConcentration(
				benzoateBulkConcentration));
		//benzoate.setBulkConcentration(new DynamicBulkConcentrationImplicit(
		//		benzoateBulkConcentration));
		// --- Oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		// ---- Create the particulates ----
		// ---H-PHB active mass
		ParticulateSpecies activeAcinetobacter = new ParticulateSpecies(
				"Acinetobacter-active-mass", densityX, Color.red);
		// ---H-EPS active mass
		ParticulateSpecies activePseudomonas = new ParticulateSpecies(
				"Pseudomonas-active-mass", densityX, Color.green);
		// ---- Create the biomass species ----
		// ----PHB producer
		// array of fixed species that constitute H-PHB biomass species
		ParticulateSpecies[] spAcinetobacter = { activeAcinetobacter };
		float[] fractionalVolumeCompositionAcinetobacter = { 1.0f };
		BiomassSpecies acinetobacter = new BiomassSpecies("Acinetobacter",
				spAcinetobacter, fractionalVolumeCompositionAcinetobacter);
		acinetobacter.setActiveMass(activeAcinetobacter);
		// ----EPS producer
		// array of fixed species that constitute H-PHB biomass species
		ParticulateSpecies[] spPseudomonas = { activePseudomonas };
		float[] fractionalVolumeCompositionPseudomonas = { 1.0f };
		BiomassSpecies pseudomonas = new BiomassSpecies("Pseudomonas",
				spPseudomonas, fractionalVolumeCompositionPseudomonas);
		pseudomonas.setActiveMass(activePseudomonas);
		// as
		// capsule
		// ---- Create the process terms ----
		// Monod and inhibition coefficients
		ProcessFactor mBenzylOH = new Saturation(benzylOH, KBenzylOH);
		ProcessFactor mBenzoate = new Saturation(benzoate, KBenzoate);
		ProcessFactor mOxygen = new Saturation(oxygen, KOxygen);
		// ---- Create the reactions ----
		// --- uptake of benzyl alcohol by acinetobacter
		Reaction benzylOHUptakeAcinetobacter = new Reaction(
				"benzylOH uptake by Acinetobacter", activeAcinetobacter,
				qBenzylOHMax, 2);
		benzylOHUptakeAcinetobacter.addFactor(mBenzylOH);
		benzylOHUptakeAcinetobacter.addFactor(mOxygen);
		// --- benzoate uptake by Pseudomonas
		Reaction benzoateUptakePseudomonas = new Reaction(
				"benzoate uptake by Pseudomonas", activePseudomonas, qBenzoate,
				2);
		benzoateUptakePseudomonas.addFactor(mBenzoate);
		benzoateUptakePseudomonas.addFactor(mOxygen);
		// ---- Assign reaction to each species through NetRaction ----
		// --- Acinetobacter active mass
		NetReaction rsAcinetobacterActive = new NetReaction(1);
		rsAcinetobacterActive.addReaction(benzylOHUptakeAcinetobacter,
				YXBenzylOH);
		activeAcinetobacter.setProcesses(rsAcinetobacterActive);
		// --- Pseudomonas active mass
		NetReaction rsPseudomonasActive = new NetReaction(1);
		rsPseudomonasActive.addReaction(benzoateUptakePseudomonas, YXBenzoate);
		activePseudomonas.setProcesses(rsPseudomonasActive);
		// --- Benzyl alcohol
		NetReaction rsBenzylOH = new NetReaction(1);
		rsBenzylOH.addReaction(benzylOHUptakeAcinetobacter, -1);
		benzylOH.setProcesses(rsBenzylOH);
		// --- Benzoate
		NetReaction rsBenzoate = new NetReaction(2);
		rsBenzoate.addReaction(benzoateUptakePseudomonas, -1);
		rsBenzoate
				.addReaction(benzylOHUptakeAcinetobacter, YBenzoateProduction);
		benzoate.setProcesses(rsBenzoate);
		// --- Oxygen
		NetReaction rsOxygen = new NetReaction(2);
		rsOxygen.addReaction(benzylOHUptakeAcinetobacter,
				-(1 - YBenzoateProduction - YXBenzylOH));
		rsOxygen.addReaction(benzoateUptakePseudomonas, -(1 - YXBenzoate));
		oxygen.setProcesses(rsOxygen);
		// ---- Add the species to system ----
		addBiomassSpecies(acinetobacter);
		addBiomassSpecies(pseudomonas);
		addSoluteSpecies(benzoate);
		addSoluteSpecies(benzylOH);
		addSoluteSpecies(oxygen);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	/*
	 * (non-Javadoc)
	 */
	protected void inoculate() {
		int[] nCells = { initialCellNumber, initialCellNumber };
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation of 1000 iterative steps, storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// set preprocessing and postprocessing steps for multigrid algorithm
		MultigridVariable.setSteps(50, 500);
		//
		ApplicationComponent app = new TwoSpecies3();
		app = new BiomassVizualizer(app);
		
		// the biomass thickness visualizer
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		//app = new SeriesVizualizer(app, thickness);
		
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// uncomment the following line for plot of produced biomass
		 app = new SeriesVizualizer(app, prod);
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// uncomment the following line for plot of total biomass in biofilm
		app = new BulkConcentrationVizualizer(app);
		// add vizualizer for solutes rates
		app = new SoluteRateSeriesVizualizer(app);
		// detached biomass
		//app = new DetachedBiomassVizualizer(app);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// set reactor dimensions
			// set the global mass balance parameters
			// --- nothing to set in this case: constant bulk concentration
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			app.addStateWriter(new PovRayWriter());
			app.addStateWriter(new SoluteConcentrationWriter());
			app.addStateWriter(new SolidsConcentrationWriter());
			app.addStateWriter(new ParticlePositionWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			// add the time constraints writer
			app.addStateWriter(new TimeConstraintsWriter());
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
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
			// set the reactor dimensions
			app.setReactorParameters(residenceTime, carrierArea, reactorVolume);
			Model.model().setFinishIterationTime(finishSimulationTime);
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
		}
		try {
			// wait for user to press start iteration
			// app.waitForStartIteratingRequest();
			// start iterating cycle
			app.startIterating();
		} catch (Exception e1) {
			e1.printStackTrace();
			System.out.println("Allocated memory:"
					+ Runtime.getRuntime().totalMemory() + " byte");
		}
		System.out.println("Simulation finished.");
	}
}