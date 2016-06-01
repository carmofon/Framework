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
import nl.tudelft.bt.model.multigrid.boundary_conditions.GranuleBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.NoBoundaryLayer;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.reaction.*;

/**
 * Tumor with two cell types
 * 
 * One macrophage placed at center, packed with cancer cells
 * Growth and migration rates are 0
 * 
 * Objective: to see how well CCs at different distances from the macrophage
 * align
 * 
 * @author Will Chang
 */
public class TestTwoCellTypesStatic extends GranuleModelHandler {
	// All the model parameters are defined here as static attributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program

	// default run with graphics
	private static boolean runWithGraphics = true;
	
	// output directory name
	protected static String outputDirectory = "/Users/wkc/xlab/results/tests";
	
	// How long to run the simulation?
	protected static float simulationFinishTime = Float.POSITIVE_INFINITY;
	
	// How many iterations to run the simulation?
	protected static int simulationFinishIterationCount = 100;
	
	// Alternatively, use max tumor radius as end criterium
	protected static float simulationFinishRadius = 180;
	
	// Max time step
	protected static float maximumTimeStep = 10f; // [h]

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
	
	//
	private static SoluteSpecies csf1;
	private static SoluteSpecies egf;
	private static SoluteSpecies oxygen;

	// Solute species
	// Oxygen (O)
	protected static float oxygenBulkConcentration = 1e-4f; // [g/L] produces smooth tumor
	//protected static float oxygenBulkConcentration = 1e-20f; // [g/L] produces fingers
	
	// CSF-1
//	protected static float CSF1BulkConcentration = 1e-6f; // [g/L] produces fingers with 1 cell type
	protected static float CSF1BulkConcentration = 0f; // [g/L] produces fingers with 1 cell type
	
	// EGF
	protected static float EGFBulkConcentration = 1e-4f; // [g/L]
//	protected static float EGFBulkConcentration = 0f; // [g/L]
	
	
	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	private static float oxygenDiffusivity = 5.76e6f; // [um^2/h]
	private static float CSF1Diffusivity = 1.78e6f; // [um^2/h]
	private static float EGFDiffusivity = 8.26e5f; // [um^2/h]
	
	// Nutrient flux parameters
	protected static float externalMassTransferCoefficient = 1e5f; // [h-1]

	//
	// Particulate species
	protected static float specificMassTumor = 3000f; // [g/L]
	protected static float specificMassMacrophage = 3000f; // [g/L]
	
	// motility diffusivity
	// mean induced speed ~ 24um/h (Goswami 2005)
	protected static float motility = 0f; // [um/h]
	
	// relative adhesive radius
	private static float RA = 1.5f; // [cell diameter]
	
	// weight of social influence on motility
	private static float tumorWS = 1f; // [dimensionless]

	// weight of individual factors
	private static float tumorWG = 1f; // [dimensionless]
	
	private static float macrophageWS = 0f;
	
	private static float macrophageWG = 1f;
	
	// weight of random polarization
	private static float WR = 0.001f;

	// Yield coefficients
	// oxygen
	private static float YXO = 0.045f; // [gX/gS]
	// yield of EGF from CSF1
	private static float EGFXCSF1 = 10f; // [gEGF/gCSF1]
	// yield of CSF1 from EGF
	private static float CSF1XEGF = 10f; // [gCSF1/gEGF]

	// Processes
	// Growth (biomass production)
	protected static float uMax = 0.0001f; // [1/h]
	
	// basal relative CSF-1 production rate
	protected static float basalRcsf1 = 0.001f;
	
	protected static float basalREGF = 0.001f;
	
	// stimulated relative CSF-1 production rate
	protected static float stimulatedRcsf1 = 0.1f;
	
	// CSF-1 production cost (this should be < 0)
	protected static float costCSF1 = 0f;
	
	// relative EGF production rate
	protected static float Regf = 0.01f;
	
	// EGF production cost
	// assume for now that macrophages can make EGF for free, since we are
	// not studying selection in macrophages (yet)
	protected static float costEGF = 0f;

	// the average mutation rate for tumor cells; between 0 and 1, preferably small
	protected static float mutationRate = 0f;
	
	// Half-saturation constant for oxygen growth
	private static float KO = 1e-4f; // [g/L]
	
	// Half-saturation constant for CSF-1 regulated behavior
	private static float KCSF1 = 1e-5f; //[g/L]
	
	// ...for EGF regulated behavior
	private static float KEGF = 1e-5f;
	
	// Gradient sensitivity constant (Tranquillo 1988)
	private static float KD = 0.02f; // [dimless]

	// Computation parameters
	protected static float systemSize = 200; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius1000 
	protected static float relativeMaximumRadius = 0.005f;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.1f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 0.05f;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]
	
	// r^2 detachment
	protected static float rdetach = 0; // 0 = NO DETACHMENT
	// r^4 detachment
	//protected static float rdetach = 1e-4f; // 0 = NO DETACHMENT

	// initial number of tumor cells in the system (inoculum)
	protected static int initialTumorCellNumber = 1000;
	
	// initial number of macrophages
	protected static int initialMacrophageNumber = 1;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// oxygen
		oxygen = new SoluteSpecies("oxygen",
				oxygenDiffusivity);
		// set up the simplest type of bulk concentration: constant
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		
		csf1 = new SoluteSpecies("csf1", CSF1Diffusivity);
		csf1.setBulkConcentration(new ConstantBulkConcentration(CSF1BulkConcentration));
		
		egf = new SoluteSpecies("egf", EGFDiffusivity);
		egf.setBulkConcentration(new ConstantBulkConcentration(EGFBulkConcentration));

		// 2a. Create the particulate species (solids)
		// wild-type tumor
		ParticulateSpecies tumorStationary = new ParticulateSpecies("stationary",
				specificMassTumor, Color.red);
		
		// mutant tumor
		ParticulateSpecies tumorMotile = new ParticulateSpecies("motile",
				specificMassTumor, Color.MAGENTA);

		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spTumor = { tumorStationary, tumorMotile };
		float[] fractionalVolumeCompositionTumor = { 1, 0 };
		
		// macrophage
		ParticulateSpecies macrophage = new ParticulateSpecies("macrophage", 
				specificMassMacrophage, Color.blue);
		ParticulateSpecies[] spMacrophage = {macrophage};
		float[] fractionalVolumeCompositionMacrophage = {1};

		// 3a. Create the biomass speciesTumor
		BiomassSpecies speciesTumor = new SwitchableVicsekBiomassSpecies("mutatorTumor",
				spTumor, fractionalVolumeCompositionTumor,
				tumorStationary, tumorMotile,
				egf, motility, RA, tumorWG, tumorWS, WR, KEGF, KD);
		speciesTumor.setShovingHierarchy(1);
		
		// assume macrophages do not mutate
		BiomassSpecies speciesMacrophage = new ZonalVicsekBiomassSpecies("macrophage", 
				spMacrophage, fractionalVolumeCompositionMacrophage, macrophage, 
				macrophage, 0, csf1, motility, 0, macrophageWG, macrophageWS, WR, KCSF1, KD);
		speciesMacrophage.setShovingHierarchy(2); // macrophage pushes tumor cells, but not v.v.
		
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor oxygenSaturation = new Saturation(oxygen, KO);
		ProcessFactor egfSaturation = new Saturation(egf, KEGF);
		ProcessFactor csf1Saturation = new Saturation(csf1, KCSF1);

		// 5. Create the reactions
		// growth
		Reaction growthStationary = new Reaction("growthStationary", tumorStationary, uMax, 1);
		growthStationary.addFactor(oxygenSaturation);
			
		// basal CSF-1 production by mutant tumor
//		Reaction basalProduceCSF1ByStationary = new Reaction(
//				"basalProduceCSF1ByStationary", tumorStationary, basalRcsf1, 0);
//		Reaction basalProduceCSF1ByMotile = new Reaction(
//				"basalProduceCSF1ByMotile", tumorMotile, basalRcsf1, 0);
		
		// basal EGF production by macrophage
		Reaction basalProduceEGF = new Reaction("basalProduceEGF", macrophage, basalREGF, 0);
				
		// stimulated CSF-1 production
		Reaction stimulatedProduceCSF1ByStationary = new Reaction(
				"stimulatedProduceCSF1ByStationary", tumorStationary, stimulatedRcsf1, 1);
		stimulatedProduceCSF1ByStationary.addFactor(egfSaturation);
		Reaction stimulatedProduceCSF1ByMotile = new Reaction(
				"stimulatedProduceCSF1ByMotile", tumorMotile, stimulatedRcsf1, 1);
		
		// EGF production
		Reaction produceEGF = new Reaction("produceEGF", macrophage, Regf, 1);
		produceEGF.addFactor(csf1Saturation);
		
		// fluxes
		Reaction fluxOxygen = new Flux("fluxOxygen", oxygen, oxygenBulkConcentration);
		Reaction fluxEGF = new Flux("fluxEGF", egf, EGFBulkConcentration);
		Reaction fluxCSF1 = new Flux("fluxCSF1", csf1, CSF1BulkConcentration);
		
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsTumorStationary = new NetReaction(2);
		rsTumorStationary.addReaction(growthStationary, 1);
//		rsTumorStationary.addReaction(basalProduceCSF1ByStationary, costCSF1);
		rsTumorStationary.addReaction(stimulatedProduceCSF1ByStationary, costCSF1);
		tumorStationary.setProcesses(rsTumorStationary);
		//
		NetReaction rsTumorMotile = new NetReaction(1);
//		rsTumorMotile.addReaction(basalProduceCSF1ByMotile, costCSF1);
		rsTumorMotile.addReaction(stimulatedProduceCSF1ByMotile, costCSF1);
		tumorMotile.setProcesses(rsTumorMotile);
		
		NetReaction rsMacrophage = new NetReaction(2);
		rsMacrophage.addReaction(basalProduceEGF, 0);
		rsMacrophage.addReaction(produceEGF, costEGF);
		macrophage.setProcesses(rsMacrophage);
		
		NetReaction rsOxygen = new NetReaction(2);
		rsOxygen.addReaction(growthStationary, -(1 / YXO));
		rsOxygen.addReaction(fluxOxygen, externalMassTransferCoefficient);
		oxygen.setProcesses(rsOxygen);
		
		NetReaction rsEGF = new NetReaction(3);
		rsEGF.addReaction(fluxEGF, externalMassTransferCoefficient);
		rsEGF.addReaction(basalProduceEGF, 1);
		rsEGF.addReaction(produceEGF, 1);
//		rsEGF.addReaction(stimulatedProduceCSF1ByStationary, -(1 / CSF1XEGF));
//		rsEGF.addReaction(stimulatedProduceCSF1ByMotile, -(1 / CSF1XEGF));
		egf.setProcesses(rsEGF);
		
		NetReaction rsCSF1 = new NetReaction(3);
//		rsCSF1.addReaction(basalProduceCSF1ByStationary, 1);
//		rsCSF1.addReaction(basalProduceCSF1ByMotile, 1);
		rsCSF1.addReaction(stimulatedProduceCSF1ByStationary, 1);
		rsCSF1.addReaction(stimulatedProduceCSF1ByMotile, 1);
//		rsCSF1.addReaction(produceEGF, -(1 / EGFXCSF1));
		rsCSF1.addReaction(fluxCSF1, externalMassTransferCoefficient);
		csf1.setProcesses(rsCSF1);
		
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesTumor);
		addBiomassSpecies(speciesMacrophage);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(egf);
		addSoluteSpecies(csf1);
	}
	
	/*
	 * (non-Javadoc)
	 * 
	 * Solve without boundary layer
	 * 
	 * @see nl.tudelft.bt.model.apps.ModelHandler#createBoundaryLayer(float)
	 */
	protected void createBoundaryLayer(float h)
		throws MultigridSystemNotSetException {
		_boundaryLayer = new NoBoundaryLayer();
		// create the boundary conditions
		MultigridVariable
		.setBoundaryConditions(new GranuleBoundaryConditions());
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
		int[] nMacrophage = {0, initialMacrophageNumber};
		int[] nTumor = {initialTumorCellNumber, 0};
		inoculateRandomly(nTumor);
		inoculateRandomly(nMacrophage);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Radius2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		Model.model().setSeed(0);
		// read arguments from the command line
		if (args.length > 0) {
			if (args.length < 5) {
				throw new RuntimeException("input arguments missing");
//						"input arguments missing: \n"
//						+ "1: output directory (CAUTION!!! directory will be erased \n"
//						+ "2: seed for random number generator \n"
//						+ "3: flag for running with graphics (1 on, 0 off) \n"
//						+ "4: mutation rate \n"
//						+ "5: signal production cost (this number should be NEGATIVE) \n" 
//						+ "6: 2D/3D (2:2D, 3:3D");
			}
			outputDirectory = args[0];
			int seed = Integer.parseInt(args[1]);
			Model.model().setSeed(seed);
			runWithGraphics = (Integer.parseInt(args[2]) == 1);
			geometry = Integer.parseInt(args[3]);
			tumorWS = Float.parseFloat(args[4]);
			System.out.println("seed = " + seed + "; tumor WS = " + tumorWS);
		}
		
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new TestTwoCellTypesStatic();
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
			app.initializeDiffusionReactionSystem(); // also innoculates
			app.initializeDetachmentFunction();
			app.intializeStateWriters(outputDirectory);
			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new ImageWriter(oxygen));
			app.addTimedStateWriter(new ImageWriter(egf));
			app.addTimedStateWriter(new ImageWriter(csf1));
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
			Model.model().setFinishIterationTime(simulationFinishTime);
			Model.model().setFinishIterationCount(simulationFinishIterationCount);
//			Model.model().setMaxRunLength(simulationFinishRadius);
//			Model.model().setMaximumTimeStep(maximumTimeStep);
//			Model.model().overrideTimeStep(maximumTimeStep);
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
