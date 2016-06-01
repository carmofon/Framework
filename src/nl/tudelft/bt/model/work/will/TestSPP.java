package nl.tudelft.bt.model.work.will;

import java.awt.Color;
import java.io.IOException;
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
import nl.tudelft.bt.model.work.benzoate.MotileBiomassSpecies;
import nl.tudelft.bt.model.work.vanni.WriteParsToFile;

/**
 * Tumor with two cell types: chemoattractant producer and adhesive
 * chemotactic migrator
 * 
 * @author Will Chang
 */
public class TestSPP extends GranuleModelHandler {
	// All the model parameters are defined here as static attributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program
	
	// How long to run the simulation?
	protected static float simulationFinishTime = Float.POSITIVE_INFINITY;

	// How many iterations to run the simulation?
//	protected static int simulationFinishIterationCount = (int) Float.POSITIVE_INFINITY;
	protected static int simulationFinishIterationCount = 3000;
	
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

	// Solute species
	// unified general bulk concentration
	protected static float bulkConcentration = 1e-4f;
	
	// diffusivity
	private static float diffusivity = 2e5f; // [um^2/h]
	
	// Particulate species (biomass)	
	protected static float specificMass = 3000f; // [g/L]
	
	// motility diffusivity
	// mean induced speed ~ 24um/h (Goswami 2005)
	protected static float motility = 576; // [um^2/h]

	// symmetric half-saturation constant
	private static float Km = 1e-5f; // [g/L]

	// concentration threshold for motility
	private static float rth = 0.25f;
	
	// Gradient sensitivity constant (Tranquillo 1988)
	private static float KD = 0.02f; // [dimless]
	
	// weight of social influence on motility
	private static float WS = 1f; // [dimensionless]

	// weight of individual factors
	private static float WG = 0f; // [dimensionless]
	
	// weight of random polarization
	private static float WR = 0f;

	// relative adhesive radius
	private static float RA = 1.5f; // [cell diameter]

	// Parameters for coloring cells by inductedness
	private static Color uninducedColor = new Color(127, 0, 0);
	private static Color inducedColor = new Color(255, 0, 0);
	
	// Computation parameters
	protected static float systemSize = 400; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.0025f;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// protected static float relativeBoundaryLayer = 0.05f;
	protected static float relativeBoundaryLayer = 0.005f;
	
	// Nutrient flux parameters
	protected static float externalMassTransferCoefficient = 0; // [h-1]

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]
	
	// r^2 detachment
	protected static float rdetach = 0f; // 0 = NO DETACHMENT

	// initial number of particles in the system (inoculum)
	protected static int cellNumber = 1000;
	
	// set a fixed time step
	private static float maxTimeStep = 0.1f;

	/**
	 * Define the solute and biomass species
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// TODO 1. Create the solutes
		// CSF-1
		SoluteSpecies solute = new SoluteSpecies("solute", diffusivity);
		solute.setBulkConcentration(new ConstantBulkConcentration(bulkConcentration));
		
		// TODO 2. Create the particulate species (solids)
		// Tumor cell active mass
		ParticulateSpecies activeBiomass = new ParticulateSpecies("activeMass",
				specificMass, Color.red);
		
		// array of fixed species that constitute tumor cell
		ParticulateSpecies[] sp = { activeBiomass };
		float[] fractionalVolumeComposition = { 1 };

		// TODO 3a. Create the tumor biomass species
		PolarizableBiomassSpecies species = new PolarizableBiomassSpecies("species", sp, 
				fractionalVolumeComposition, solute, motility, RA, WG, WS, WR, Km, KD, rth);
		// color cells according to inductedness
		species.
		setInducibleColor(solute, Km, inducedColor, uninducedColor);
		// if no graphics, set simulation to end on periodic boundary crossing
//		if (!runWithGraphics) species.setStopOnPeriodicBoundaryCrossing();
		
		Reaction produceSolute = new Reaction("produceSolute", activeBiomass, 0, 0);
		
		NetReaction rsBiomass = new NetReaction(1);
		rsBiomass.addReaction(produceSolute, 0);
		activeBiomass.setProcesses(rsBiomass);
		
		Reaction fluxSolute = new Flux("fluxSignal", solute, bulkConcentration);
		
		// assign reaction stoichiometry to the solutes
		NetReaction rsSolute = new NetReaction(1);
		rsSolute.addReaction(fluxSolute, externalMassTransferCoefficient);
		solute.setProcesses(rsSolute);
		
		// TODO 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(species);
		addSoluteSpecies(solute);
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
		int[] nCells = { cellNumber };
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here.
		DetachmentSpeedFunction df = new Radius2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}
	
	/*
	 * Parameter writing object c/o Vanni Bucci
	 */
	public static void writeParameters() throws IOException {
		String path = outputDirectory + "/parameters.txt";
		WriteParsToFile parameterWriter = new WriteParsToFile(path, true);
		parameterWriter.writeToFile("bulkConcentration\t" + bulkConcentration);
		parameterWriter.writeToFile("diffusivity\t" + diffusivity);
		parameterWriter.writeToFile("externalMassTransferCoefficient\t" + externalMassTransferCoefficient);
		parameterWriter.writeToFile("cellNumber\t" + cellNumber);
		parameterWriter.writeToFile("KD\t" + KD);
		parameterWriter.writeToFile("Km\t" + Km);
		parameterWriter.writeToFile("motility\t" + motility);
		parameterWriter.writeToFile("WS\t" + WS);
		parameterWriter.writeToFile("WG\t" + WG);
		parameterWriter.writeToFile("WR\t" + WR);
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		Model.model().setSeed(1);
		// read arguments from the command line
		if (args.length > 0) {
			if (args.length < 7) {
				throw new RuntimeException(
						"input arguments missing: \n"
								+ "1: output directory (CAUTION!!! directory will be erased) \n"
								+ "2: flag for running with graphics (1 on, 0 off) \n"
								+ "3: 2D/3D (2:2D, 3:3D)"
								+ "4: seed for random number generator \n"
								+ "5: WS \n"
								+ "6: WG \n"
								+ "7: WR \n"
				);
			}
			outputDirectory = args[0];
			runWithGraphics = (Integer.parseInt(args[1]) == 1);
			geometry = Integer.parseInt(args[2]);
			int seed = Integer.parseInt(args[3]);
			Model.model().setSeed(seed);
			WS = Float.parseFloat(args[4]);
			WG = Float.parseFloat(args[5]);
			WR = Float.parseFloat(args[6]);
		}
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new TestSPP();
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
			// write the parameters
//			writeParameters();
//			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new ImageWriter());
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