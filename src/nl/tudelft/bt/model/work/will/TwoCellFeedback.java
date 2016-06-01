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
 * Tumor with two cell types: chemoattractant producer and adhesive
 * chemotactic migrator
 * 
 * @author Will Chang
 */
public class TwoCellFeedback extends GranuleModelHandler {
	// All the model parameters are defined here as static attributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program
	
	// How long to run the simulation?
	protected static float simulationFinishTime = Float.POSITIVE_INFINITY;
//	protected static float simulationFinishTime = 1f;
	
	// How many iterations to run the simulation?
//	protected static int simulationFinishIterationCount = (int) Float.POSITIVE_INFINITY;
	protected static int simulationFinishIterationCount = 200;
	
	// Alternatively, use max tumor radius as end criterium
	protected static float simulationFinishRadius = Float.POSITIVE_INFINITY;

	// output directory name
	protected static String outputDirectory = "/Users/wkc/xlab/results/tests";
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
	private static float maxTimeStep = 0.01f;

	// Solute species
	// unified general bulk concentration
//	protected static float bulkConcentration = 1e-12f;
	protected static float bulkConcentration = 0f;
	
	// CSF-1 bulk concentration
	protected static float CSF1BulkConcentration = bulkConcentration; // [g/L]
	// in experiments macrophages induced with roughly 4-6e-5 g/L
//	protected static float CSF1BulkConcentration = 6e-5f; 
	
	// EGF bulk concentration
	protected static float EGFBulkConcentration = bulkConcentration; // [g/L]
	
	// diffusivity
	// CSF-1 diffusivity (calculated from molecular weight)
	private static float CSF1Diffusivity = 1.78e6f; // [um^2/h]

	// EGF diffusivity (calculated from molecular weight)
	private static float EGFDiffusivity = 8.26e5f; // [um^2/h]
	
	// Particulate species (biomass)	
	protected static float specificMassX = 3000f; // [g/L]

	// unified max chemoattractant production rate
	protected static float maxR = 0.005f; // [dimensionless] (fraction of biomass)
	
	// CSF-1 production rate
	protected static float Rc = maxR; // [1/h]

	// EGF production rate
	protected static float Re = maxR; // [1/h]	
	
	// relative minimum response (0-1)
	protected static float relativeMinResponse = 0.3f;
	
	// relative minimum CSF-1 production
	protected static float minimumCSF1Response = relativeMinResponse;
	
	// relative minimum EGF production
	protected static float minimumEGFResponse = relativeMinResponse;
	
	// motility diffusivity
	// mean induced speed ~ 24um/h (Goswami 2005)
	protected static float motility = 576; // [um^2/h]
	
	// macrophage motility
	protected static float macrophageMotility = motility; // [um/h]
	
	// tumor cell motility
	protected static float tumorMotility = motility; // [um/h]
	
	// symmetric half-saturation constant
	private static float Km = 1e-5f; // [g/L]
	
	// Half-saturation constant for EGF-induced processes
	private static float Ke = Km; // [g/L]
	
	// Half-saturation constant for CSF-1 induced processes
	private static float Kc = Km; // [g/L]
	
	// Gradient sensitivity constant (Tranquillo 1988)
	private static float KD = 0.002f; // [dimless]
	
	// half-saturation directedness constant for tumor cells
	private static float tumorKDirectedness = KD; // [dimensionless]
	
	// half-saturation directedness constant for macrophage
	private static float macrophageKDirectedness = KD; // [dimensionless]
	
	// weight of social influence on motility
	private static float WS = 0f; // [dimensionless]
	// set to 0 for fully individual migration
//	private static float WS = f;
	
	// tumor cell social factors weight
	private static float tumorWS = WS; // [dimless]
	
	// tumor cell individual factors weight
	private static float tumorWG = 1 - tumorWS; // [dimless]
	
	// macrophage social factors weight
	private static float macrophageWS = WS; // [dimensionless]
	
	// macrophage individual factors weight
	private static float macrophageWG = 1 - macrophageWS; // [dimensionless]
	
	// relative adhesive radius
	private static float RA = 1.5f; // [cell diameter]
	
	// relative adhesive radius of tumor cells
	private static float tumorRelativeAdhesiveRadius = RA; // [dimless]
	
	// relative adhesive radius of macrophages
	private static float macrophageRelativeAdhesiveRadius = RA; // [dimless]
	
	// Parameters for coloring cells by inductedness
	
	// tumor
	private static Color tumorUninducedColor = new Color(127, 0, 0);
	private static Color tumorInducedColor = new Color(255, 0, 0);
	
	// macrophage
	private static Color macrophageUninducedColor = new Color(0, 127, 0);
	private static Color macrophageInducedColor = new Color(0, 255, 0);
	
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
	// protected static float relativeBoundaryLayer = 0.05f;
	protected static float relativeBoundaryLayer = 0.005f;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]
	
	// r^2 detachment
	protected static float rdetach = 0f; // 0 = NO DETACHMENT

	// initial number of particles in the system (inoculum)
	protected static int tumorCellNumber = 1500;
//	protected static int tumorCellNumber = 0;
	protected static int macrophageNumber = 100;
	
	// initial radius of macrophage inoculum
	protected static float macrophageInoculateRadius = 45f;
	
	// initial radius of tumor inoculum
	protected static float tumorInoculateRadius = 40f;

	/**
	 * Define the solute and biomass species
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// TODO 1. Create the solutes
		// CSF-1
		SoluteSpecies CSF1 = new SoluteSpecies("CSF-1", CSF1Diffusivity);
		CSF1.setBulkConcentration(new ConstantBulkConcentration(CSF1BulkConcentration));
		
		// EGF
		SoluteSpecies EGF = new SoluteSpecies("EGF", EGFDiffusivity);
		EGF.setBulkConcentration(new ConstantBulkConcentration(EGFBulkConcentration));
		
		// TODO 2. Create the particulate species (solids)
		// Tumor cell active mass
		ParticulateSpecies activeTumor = new ParticulateSpecies("activeTumor",
				specificMassX, Color.red);

		// Macrophage active mass
		ParticulateSpecies activeMacrophage = new ParticulateSpecies("activeMacrophage",
				specificMassX, Color.blue);
		
		// array of fixed species that constitute tumor cell
		ParticulateSpecies[] spTumor = { activeTumor, activeMacrophage };
		float[] fractionalVolumeCompositionTumor = { 1, 0 };
		// array of fixed species that constitute macrophage
		ParticulateSpecies[] spMacrophage = { activeTumor, activeMacrophage };
		float[] fractionalVolumeCompositionMacrophage = { 0, 1 };

		// TODO 3a. Create the tumor biomass species
		BiomassSpecies speciesTumor = new AdhesiveBiomassSpecies("tumor",
				spTumor, fractionalVolumeCompositionTumor, EGF,
				tumorMotility, Ke, tumorKDirectedness, minimumEGFResponse,
				tumorWS, tumorWG, tumorRelativeAdhesiveRadius);
		speciesTumor.setShovingHierarchy(2);
		// color cells according to inductedness
		speciesTumor.setInducibleColor(EGF, Ke, tumorInducedColor, tumorUninducedColor);
		
		// 3b. Create the macrophage biomass species
		BiomassSpecies speciesMacrophage = new AdhesiveBiomassSpecies("macrophage", 
				spMacrophage, fractionalVolumeCompositionMacrophage, CSF1,
				macrophageMotility, Kc, macrophageKDirectedness, minimumCSF1Response,
				macrophageWS, macrophageWG, macrophageRelativeAdhesiveRadius);
		speciesMacrophage.setShovingHierarchy(1);
		speciesMacrophage.setInducibleColor(CSF1, Kc, macrophageInducedColor, macrophageUninducedColor);
		
		// TODO 4. Create the Reaction factors, Monod and inhibition coefficients
		// CSF-1
		ProcessFactor CSF1SaturationWithFloor = new SaturationWithFloor(CSF1, Kc, minimumCSF1Response);		
		
		// EGF
		ProcessFactor EGFSaturation = new Saturation(EGF, Ke);
		
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		ProcessFactor EGFSaturationWithFloor = new SaturationWithFloor(EGF, Ke, minimumEGFResponse);
		
		// TODO 5. Create the reactions

		// CSF-1 production
		Reaction produceCSF1 = new Reaction("produceCSF-1", activeTumor, Rc, 1);
		// add factor for induction by CSF-1
		produceCSF1.addFactor(EGFSaturationWithFloor);
		
		// EGF production
		Reaction produceEGF = new Reaction("produceEGF", activeMacrophage, Re, 1);
		// add factor for induction by EGF
		produceEGF.addFactor(CSF1SaturationWithFloor);
		
		// TODO 6. Assign reaction to the species through ReactionStoichiometries
		
		// active mass
		NetReaction rsTumorActive = new NetReaction(1);
		rsTumorActive.addReaction(produceCSF1, 0);
		activeTumor.setProcesses(rsTumorActive);

		NetReaction rsMacrophageActive = new NetReaction(1);
		rsMacrophageActive.addReaction(produceEGF, 0);
		activeMacrophage.setProcesses(rsMacrophageActive);
		// This defines that biomass growth rate is 1*rX
		
		// assign reaction stoichiometry to the solutes
		// CSF-1
		NetReaction rsCSF1 = new NetReaction(1);
		rsCSF1.addReaction(produceCSF1, 1);
		CSF1.setProcesses(rsCSF1);
		
		// EGF
		NetReaction rsEGF = new NetReaction(1);
		rsEGF.addReaction(produceEGF, 1);
		EGF.setProcesses(rsEGF);
		
		// TODO 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesTumor);
		addBiomassSpecies(speciesMacrophage);
		addSoluteSpecies(CSF1);
		addSoluteSpecies(EGF);
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
		int[] nTumor = { tumorCellNumber, 0 };
		int[] nMacrophage = { 0, macrophageNumber };
		inoculateRandomlyAtRadius(nMacrophage, macrophageInoculateRadius);
		inoculateRandomlyInsideRadius(nTumor, tumorInoculateRadius);
//		inoculateRandomly(nTumor);
//		inoculateRandomly(nMacrophage);
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
		Model.model().setSeed(1);
		// read arguments from the command line
		if (args.length > 0) {
			if (args.length < 5) {
				throw new RuntimeException(
						"input arguments missing: \n"
								+ "1: output directory (CAUTION!!! directory will be erased) \n"
								+ "2: flag for running with graphics (1 on, 0 off) \n"
								+ "3: 2D/3D (2:2D, 3:3D)"
								+ "4: seed for random number generator \n"
								+ "5: relative boundary layer \n"
				);
			}
			outputDirectory = args[0];
			runWithGraphics = (Integer.parseInt(args[1]) == 1);
			geometry = Integer.parseInt(args[2]);
			int seed = Integer.parseInt(args[3]);
			Model.model().setSeed(seed);
			relativeBoundaryLayer = Float.parseFloat(args[4]);
		}
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new TwoCellFeedback();
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
			Model.model().setMaximumTimeStep(maxTimeStep);
			Model.model().setFinishIterationTime(simulationFinishTime);
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
//