package nl.tudelft.bt.model.tests.granule.massbalance;

import java.awt.Color;
import java.util.Iterator;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Radius2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.reaction.*;

/**
 * Tests a system with a single biomass particle in granule geometry with the
 * dynamics of the bulk concentration of sybstrate derived from mass balances to
 * the the biomass produced
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class TestMassBalance1Cell extends GranuleModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	//output directory name
	protected static String outputDirectory;
	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	//Reaction parameters
	//Phospate accumulating organisms (XPAO)
	protected static float qSMaxPAO = 1f; //[gCOD-S/gCOD-XPAO/h]
	private static float KPAOS = 4e-3f; //[gCOD/L]
	private static float YPhbAerobicPAO = 1f; //[gCOD-PHB/gCOD-XPAO]

	// geometry (default is 3D) - change value to 2 for 2D
	protected static int geometry = 2;

	//Reactor properties
	protected static float reactorVolume = 3; // [L]
	protected static float feedFraction = 0.5f; // [dimensionless]

	//duration of the SBR cycle
	protected static float cycleTime = 30; //[h]
	//The SBR sycle
	private static SbrBulkConcentration.SbrCycle _cycle = new SbrBulkConcentration.SbrCycle(
			cycleTime);
	private static int _nCycles = 400; // write full data concerning cycle every
	// _nCycles

	//Feed composition
	protected static float substrateFeedConcentration = 396e-3f; //[gCOD/L]
	private static float substrateDiffusivity = 4e6f; //[um2/h]

	//Leave these values
	protected static float nConcentrationsPrecision = 0.3e-3f; //[gN/L]
	protected static float sConcentrationsPrecision = 4e-3f; //[gCOD/L]
	protected static float pConcentrationsPrecision = 0.15e-3f; //[gP/L]

	//
	//Particulate species (biomass H)
	protected static float specificMassBiomass = 150f; //[gCOD-H/L]
	protected static float specificMassPolymers = specificMassBiomass * 1e5f; //[gCOD-PHB/L]

	// Computation parameters
	//Size of computational volume (size of size of square)
	protected static float systemSize = 1600; // [um]

	//relativeMaximumRadius defines the maximum radius of the biomass particles
	//in relation to the system size
	//the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.007f;

	//Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 10 um
	protected static float relativeBoundaryLayer = 10 / systemSize;

	protected static float maximumGranuleRadius = 550; //[um]

	// other model parameters
	protected static int gridSide = 33; // multigrid grid side
	//Don't change this

	protected static float kShov = 1.0f; // shoving parameter[dim/less]
	//Don't change this

	//detachment rate
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]
	//leave at zero to form round granules

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumberXPAO = 1;

	//iteration finish time
	protected static float simulationFinishTime = cycleTime; //[h]

	//outpute (write results to file) every:
	protected static float outputEvery = 15.0f; //[h]

	//Computational volume multiplier
	protected static float nComp = 25.93f * 6.7e5f; //[dimensionless]
	//factor 25.93f makes substrate loading per biomass the same as in real
	// reactor

	///END OF PARAMETERS

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		//1. Create the solutes
		//substrate (S)
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		substrate.setBulkConcentration(new SbrBulkConcentration(
				substrateFeedConcentration, feedFraction,
				sConcentrationsPrecision, _cycle));
		//2. Create the particulate species (solids)
		//Pao (PAO)
		Color paoColor = new Color(204, 0, 204);
		ParticulateSpecies activePAO = new ParticulateSpecies("activePAO",
				specificMassBiomass, paoColor);
		// array of fixed species that constitute speciesH (in this case,
		ParticulateSpecies[] spPAO = {activePAO};
		float[] fractionalVolumeCompositionPAO = {1};
		//3. Create the biomass species
		BiomassSpecies speciesPAO = new BiomassSpecies("speciesPAO", spPAO,
				fractionalVolumeCompositionPAO);
		speciesPAO.setActiveMass(activePAO);
		//4. Create the Reaction factors, Monod and inhibition coefficients
		//for PAO
		ProcessFactor mPAOS = new Saturation(substrate, KPAOS);
		//for all organisms
		//PAO
		//substrate uptake
		Reaction sUptakePAO = new Reaction("sUptakePAO", activePAO, qSMaxPAO, 1);
		sUptakePAO.addFactor(mPAOS);
		//6. Assign reaction to the species through ReactionStoichiometries
		// active mass PAO
		NetReaction rsPaoActive = new NetReaction(1);
		rsPaoActive.addReaction(sUptakePAO, YPhbAerobicPAO);
		activePAO.setProcesses(rsPaoActive);
		//substrate (S)
		NetReaction rsSubstrate = new NetReaction(1);
		rsSubstrate.addReaction(sUptakePAO, -1);
		substrate.setProcesses(rsSubstrate);
		//7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesPAO);
		addSoluteSpecies(substrate);
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
		int[] nCells = {initialParticleNumberXPAO};
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Radius2MassDetachment(kdetach);
		setDetachmentHandler(df);
		// set the maximum granule radius
		try {
			Model.model().setMaximumBiofilmHeight(maximumGranuleRadius);
		} catch (InvalidValueException e) {
			System.out.println(e);
			System.exit(-1);
		}
	}

	/**
	 * Set the reactors dimensions
	 */
	private static void setTheReactorDimensions(ApplicationComponent app) {
		float rVIM = reactorVolume * 1e15f; //reactor volume in cubic
		// micrometer
		float carrierArea = nComp
				* Model.model().getComputationalVolumeCarrierArea();
		app.setReactorParameters(cycleTime, carrierArea, rVIM);
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line, unless
		//program is running localy (in my laptop)
		boolean runWithGraphics = false;
		if (args.length != 2) {
			throw new RuntimeException("Program arguments missing: "
					+ "2 program arguments should be supplied"
					+ " (1 - output directory,"
					+ " 2 - flag for running with graphics)");
		}
		// parse the input arguments
		outputDirectory = args[0];
		int arg1 = Integer.parseInt(args[1]);
		switch (arg1) {
			case 0 :
				runWithGraphics = false;
				break;
			case 1 :
				runWithGraphics = true;
				break;
			default :
				throw new RuntimeException("second program"
						+ " argument must be 0"
						+ " (for running with no graphics) "
						+ "or 1 (for running with graphics)");
		}
		//set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new TestMassBalance1Cell();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the biovolume series
		VaribleSeries biovolume = new BiovolumeSeries();
		VaribleSeries[] runLengthSeries = {new RunLengthXSeries(),
				new RunLengthYSeries(), new RunLengthZSeries()};
		// The following code will be omitted if no vizuals are desired
		if (runWithGraphics) {
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// bulk concentrations
			app = new BulkConcentrationVizualizer(app);
			// finally, the controller must be the last decorator to add
			app = new VizualModelControler(app);
		}
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// --- nothing to set in this case: constant bulk concentration
			//initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			//Pov witer is added twice
			app
					.addStateWriter(new TimedStateWriterDecorator(
							new PovRayWriter()));
			app
					.addStateWriter(new TimedStateWriterDecorator(
							new PovRayWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SoluteConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new ParticlePositionWriter()));
			//app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(biovolume);
			spw.addSeries(runLengthSeries[0]);
			spw.addSeries(runLengthSeries[1]);
			spw.addSeries(runLengthSeries[2]);
			spw.addSeries(Model.detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(Model.detachedBiomassContainer()
					.getErodedBiomassSeries());
			spw.addSeries(Model.detachedBiomassContainer()
					.getSloughedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			// intialize the reactor dimensions
			setTheReactorDimensions(app);
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
			// write the header in parameters.txt file
			spw.initializeParametersWriting();
		} catch (ModelException e) {
			System.out.println(e);
			e.printStackTrace();
			System.exit(-1);
		}
		try {
			// start iterating cycle
			//Model.model().setCompulsoryTimeStep(outputEvery);
			Model.model().setFinishIterationTime(simulationFinishTime);
			//app.waitForStartIteratingRequest();
			app.startIterating();
		} catch (Exception e1) {
			app.forceWriteState();
			e1.printStackTrace();
			System.out.println(e1);
		}
		System.out.println("Simulation finished.");
	}
}