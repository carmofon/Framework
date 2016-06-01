package nl.tudelft.bt.model.work.will;

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
import nl.tudelft.bt.model.detachment.levelset.functions.Radius2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.work.benzoate.MotileBiomassSpecies;;

/**
 * Tumor with two cell types
 * For now X will be producer and Y will be consumer
 * 
 * @author Will Chang
 */
public class TestChemotaxis extends GranuleModelHandler {
	// All the model parameters are defined here as static attributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory = "/Users/wkc/xlab/results/tests";
	
	// How long to run the simulation?
	protected static float simulationFinishTime = Float.POSITIVE_INFINITY;
	
	// How many iterations to run the simulation?
	protected static int simulationFinishIterationCount = (int) Float.POSITIVE_INFINITY;
	
	// Alternatively, use max tumor radius as end criterium
	protected static float simulationFinishRadius = Float.POSITIVE_INFINITY;

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.
	// The output directory is were the program will store all the results.
	// Choose a path to an existing folder in your system.
	// EXAMPLE: if you choose "e:/results/example1/" directory "e:\results" must
	// exist in your computer. The subdirectory "example1" will be created
	// if it is non-existant. If it exists, its contents will be deleted
	// during the program initialization

	// geometry (default is 2D) - change value to 3 for 3D
	protected static int geometry = 2;
	
	// default run with graphics?
	private static boolean runWithGraphics = true;
	
	// set a fixed time step
	private static float fixedTimeStep = 10f;

	// Solute species
	// Oxygen (O)
	protected static float oxygenBulkConcentration = 0f; // [g/L] produces smooth tumor
	//protected static float oxygenBulkConcentration = 1e-20f; // [g/L] produces fingers
	
	// Signal (S)
	//protected static float signalBulkConcentration = 1e-20f; // [g/L]
	protected static float signalBulkConcentration = 1e-6f;
	
	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	private static float oxygenDiffusivity = 5.76e6f; // [um^2/h]
	//private static float signalDiffusivity = 5.76e5f; // [um^2/h]
	private static float signalDiffusivity = 5.76e7f;
	
	//
	// Particulate species (biomass X)
	protected static float specificMassX = 3000f; // [g/L]

	// Yield coefficients
	// oxygen
	private static float YXO = 0.045f; // [gX/gS]
	// growth signal
	// private static float YXS = 0.045f; // [gX/gS]
	private static float YXS = 0.5f;

	// Processes
	// Growth (biomass production)
	protected static float uMax = 0.0000001f; // [1/h]
	// protected static float uMax = 0.0547f; //[1/h]
	
	// Growth signal production rate
	protected static float Rs = 10f; // [g/h]
	//protected static float Rs = 0f;
	
	// Growth signal production cost (this should be < 0 for public good production)
	protected static float c = 0f;

	// the "diffusivity" of the migrator biomass species
	protected static float migratorMotility = 0.01f; // [um^2/h]
	
	// Half-saturation constant for oxygen growth
	private static float KO = 3.5e-4f; // [g/L]
	
	// Half-saturation constant for growth signal growth
	private static float KS = 3.5e-5f; //[g/L]
	
	// half-max spread constant for migrators
	private static float migratorKspread = 0f; // [um^-1]
	//private static float migratorKspread = 100f;
	
	// half-max speed constant for migrators
	private static float migratorKspeed = 1e-6f; // [g/(L*um)]

	// Computation parameters
	protected static float systemSize = 200; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.005f;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	protected static float relativeBoundaryLayer = 0.05f;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]
	
	// r^2 detachment
	protected static float rdetach = 0f; // 0 = NO DETACHMENT

	// initial number of particles in the system (inoculum)
	protected static int initialMigratorNumber = 50;
	protected static int initialAttractorNumber = 20;
	
	// initial radius of inoculum
	protected static float inoculumRadius = 75f;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 * 
	 * Chemotaxis test model: only signal
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
//		// oxygen
//		SoluteSpecies oxygen = new SoluteSpecies("oxygen",
//				oxygenDiffusivity);
//		// set up the simplest type of bulk concentration: constant
//		oxygen.setBulkConcentration(new ConstantBulkConcentration(
//				oxygenBulkConcentration));
		// growth signal
		SoluteSpecies signal = new SoluteSpecies("chemoattractant", signalDiffusivity);
		signal.setBulkConcentration(new ConstantBulkConcentration(signalBulkConcentration));

		// 2a. Create the particulate species (solids)
		// Chemotactic cell active mass
		ParticulateSpecies activeMigrator = new ParticulateSpecies("activeMigrator",
				specificMassX, Color.red);

		// 2b. Create the particulate species (solids)
		// Attractant secretor cell active mass
		ParticulateSpecies activeAttractor = new ParticulateSpecies("activeAttractor",
				specificMassX, Color.blue);
		// array of fixed species that constitute migrator
		ParticulateSpecies[] spMigrator = { activeMigrator, activeAttractor };
		float[] fractionalVolumeCompositionMigrator = { 1, 0 };
		// array of fixed species that constitute attractor
		ParticulateSpecies[] spAttractor = { activeMigrator, activeAttractor };
		float[] fractionalVolumeCompositionAttractor = { 0, 1 };

		// 3a. Create the migrator biomass species
		ChemotacticBiomassSpecies speciesMigrator = new ChemotacticBiomassSpecies("migrator",
				spMigrator, fractionalVolumeCompositionMigrator, signal,
				migratorMotility, migratorKspeed, migratorKspread);
		//speciesMigrator.setShovingHierarchy(1);
		
		// 3b. Create hte attractor biomass species
		BiomassSpecies speciesAttractor = new BiomassSpecies("attractor", spAttractor, 
				fractionalVolumeCompositionAttractor);

		
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		//ProcessFactor mS = new Saturation(oxygen, KO);
		ProcessFactor signalS = new Saturation(signal, KS);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate

		// 5. Create the reactions
		// growth
		Reaction growthMigrator = new Reaction("growthMigrator", activeMigrator, uMax, 1);
//		Reaction growthX = new Reaction("growthX", activeX, uMax, 1);
//		//growthX.addFactor(mS);
//		growthX.addFactor(signalS);
//
//		Reaction growthY = new Reaction("growthY", activeY, uMax, 1);
//		//growthY.addFactor(mS);
//		growthY.addFactor(signalS);
			
		// growth signal production
		Reaction produceSignal = new Reaction("produceSignal", activeAttractor, Rs, 0);
		Reaction consumeSignal = new Reaction("consumeSignal", activeMigrator, 1, 0);

		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass
		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		
		NetReaction rsAttractorActive = new NetReaction(1);
		rsAttractorActive.addReaction(produceSignal, c);
		activeAttractor.setProcesses(rsAttractorActive);

		NetReaction rsMigratorActive = new NetReaction(1);
		rsMigratorActive.addReaction(growthMigrator, 1);
		activeMigrator.setProcesses(rsMigratorActive);
//		NetReaction rsYactive = new NetReaction(1);
//		rsYactive.addReaction(growthY, 1);
//		activeY.setProcesses(rsYactive);
		//
		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// oxygen
		//NetReaction rsOxygen = new NetReaction(2);
		//rsOxygen.addReaction(growthX, -(1 / YXO));
		//rsOxygen.addReaction(growthY, -(1 / YXO));
		//oxygen.setProcesses(rsOxygen);
		// This defines that oxygen consumption rate is -(1 / YXO)*(rX + rY)
		// growth signal
		NetReaction rsSignal = new NetReaction(2);
//		rsSignal.addReaction(growthX, -(1 / YXS));
//		rsSignal.addReaction(growthY, -(1 / YXS));
		rsSignal.addReaction(produceSignal, 1);
		rsSignal.addReaction(consumeSignal, -(1 / YXS));
		signal.setProcesses(rsSignal);
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesMigrator);
		addBiomassSpecies(speciesAttractor);
		//addSoluteSpecies(oxygen);
		addSoluteSpecies(signal);
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
		int[] nCells = { initialMigratorNumber, initialAttractorNumber };
		inoculateRandomlyInsideRadius(nCells, inoculumRadius);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here.
		DetachmentSpeedFunction df = new Radius2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// no command line arguments because this is supposed to be run locally
//		// read arguments from the command line
//		if (args.length < 6) {
//			throw new RuntimeException(
//					"input arguments missing: \n"
//							+ "1: output directory (CAUTION!!! directory will be erased \n"
//							+ "2: seed for random number generator \n"
//							+ "3: flag for running with graphics (1 on, 0 off) \n"
//							+ "4: mutation rate \n"
//							+ "5: signal production cost (this number should be NEGATIVE) \n" 
//							+ "6: 2D/3D (2:2D, 3:3D");
//		}
//		
//		outputDirectory = args[0];
//		int seed = Integer.parseInt(args[1]);
//		Model.model().setSeed(seed);
//		boolean runWithGraphics = (Integer.parseInt(args[2]) == 1);
//		mutationRate = Float.parseFloat(args[3]);
//		c = Float.parseFloat(args[4]);
//		geometry = Integer.parseInt(args[5]);
		
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new TestChemotaxis();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new RunLengthXSeries();
		// also record the Y run length
		VariableSeries width = new RunLengthYSeries();
		
		// The following code will be omitted if no vizuals are desired
		if (runWithGraphics) {
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// the biomass thickness visualizer
			app = new SeriesVizualizer(app, thickness);
		}
		
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// --- nothing to set in this case: constant bulk concentration
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new SoluteConcentrationWriter());
			app.addTimedStateWriter(new SolidsConcentrationWriter());
			app.addTimedStateWriter(new ParticlePositionWriter());
			// app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writer
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(width);
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getErodedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getSloughedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
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
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
		}
		try {
			// start iterating cycle
			// Model.model().setFinishIterationTime(simulationFinishTime);
			Model.model().setFinishIterationCount(simulationFinishIterationCount);
			Model.model().setMaxRunLength(simulationFinishRadius);
			Model.model().overrideTimeStep(fixedTimeStep);
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		try {
			System.out.println("Writing the output at end of simulation");
			app.forceWriteStateWithoutSerializing();
		} catch (ModelException e2) {
			System.err.println("Error trying to write output");			
		}
		System.out.println("Simulation finished.");
	}
}
