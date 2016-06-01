package nl.tudelft.bt.model.tests.gfp;

import java.awt.Color;
import java.util.Iterator;
import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Height2VolumetricDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;

/**
 * GFP producing construct
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Gfp extends ModelHandler {
	// output directory name (
	protected static String outputDirectory = "C:/joao/results/GFP/earlyAddition";

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

	// //Microbial parameters
	// Maximum specific growth rate
	protected static float uMax = 0.5f; // [gX/gX/h]

	// Substrate saturation constant
	private static float KS = 4e-4f; // [gS/L]

	// Yield of EPS on produced biomass
	private static float YEPS = 0.5f; // [gEPS/gX]

	// Yield of consumed substrate per biomass produced
	private static float YS = 1.63f; // [gS/gX]

	// diffusivity of solutes
	private static float substrateDiffusivity = 4.2e6f; // [um^2/h]

	// Concentration of solutes
	protected static float substrateBulkConcentration = 1e-3f; // [gO/L]

	// Specific mass of particulates
	protected static float specificMassX = 200f; // [gX/L]

	protected static float specificMassEPS = specificMassX / 6f; // [gEPS/L]

	// GFP properties
	private static float YGFP = 50f; // [gGFP/gX]

	// very high density so that effect is not felt on volume of particles
	protected static float specificMassGfP = specificMassX * 1e10f; // [gGFP/L]

	// detachment
	protected static float rdetach = 1e-2f;

	// uncomment
	// inoculation
	protected static int initialParticleNumber = 200;

	// Numeric parameters
	protected static float systemSize = 200; // [um]

	// grid resolution
	protected static int gridSide = 33; // multigrid grid side

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 2 / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	protected static float relativeBoundaryLayer = 20 / systemSize;

	// Shoving parameter
	protected static float kShov = 1.0f;

	// outpute (write results to file) every:
	protected static float outputEvery = 2f; // [h]

	protected static float finishIteratinTime = 100f; // [h]

	// early addition
	 protected static float tInducer = 0f; // [h]
	 protected static float durationInducer = 3f; // [h]
	// late addition
//	protected static float tInducer = 70f; // [h]
//
//	protected static float durationInducer = finishIteratinTime - tInducer; // [h]

	// eps and eps_star must be atributes since they are also used
	// todefine the detachment function
	private ParticulateSpecies _activeX;

	private ParticulateSpecies _eps;

	private ParticulateSpecies _gfp;

	private BiomassSpecies _speciesX;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// //1. Create the solutes
		// substrate (S)
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));
		// 2. Create the particulate species (soliids)
		// X active mass
		_activeX = new ParticulateSpecies("activeX", specificMassX, new Color(
				0, 0, 0.5f));
		// EPS
		_eps = new ParticulateSpecies("EPS", specificMassEPS, Color.gray);
		// GFP
		_gfp = new ParticulateSpecies("GFP", specificMassGfP, Color.green);
		// array of fixed species that constitute speciesX (in this case
		// active mass, EPS and EPS*)
		ParticulateSpecies[] spX = { _activeX, _eps, _gfp };
		// volunmetric fractional components
		float fEpsV = 1 / (1 + specificMassEPS / (YEPS * specificMassX));
		float[] fractionalVolumeCompositionH1 = { 1 - fEpsV, fEpsV, 0 };
		// 3. Create the biomass species
		_speciesX = new BiomassSpecies("speciesX", spX,
				fractionalVolumeCompositionH1);
		_speciesX.setActiveMass(_activeX);
		ParticulateSpecies[] eps_components = { _eps };
		_speciesX.setEpsMass(eps_components);
		// speciesX.getColorFromGrowth();
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mS = new Saturation(substrate, KS);
		// 5. Create the reactions
		// growth
		Reaction growth = new Reaction("growth", _activeX, uMax, 1);
		growth.addFactor(mS);
		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactive = new NetReaction(1);
		rsXactive.addReaction(growth, 1);
		_activeX.setProcesses(rsXactive);
		// EPS
		NetReaction rsEps = new NetReaction(1);
		rsEps.addReaction(growth, YEPS);
		_eps.setProcesses(rsEps);
		// GFP
		NetReactionWithDynamicCoefficient rsGfp = new NetReactionWithDynamicCoefficient(
				1);
		rsGfp.addReaction(growth, YGFP, tInducer, durationInducer);
		_gfp.setProcesses(rsGfp);
		//
		// assign reaction stoichiometry to the solutes
		// substrate
		NetReaction rsSubstrate = new NetReaction(1);
		rsSubstrate.addReaction(growth, -YS);
		substrate.setProcesses(rsSubstrate);
		// 7. add the solute species and the biomass
		// species (which contain the
		// particulate species) to system
		addBiomassSpecies(_speciesX);
		addSoluteSpecies(substrate);
		// Implement the dynamic behaviour of the ceofficient
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
		int[] nCells = { initialParticleNumber };
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here.
		// DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		float rd = rdetach / (YEPS + 1)
				* (YEPS / specificMassEPS + 1 / specificMassX)*2;
		DetachmentSpeedFunction df = new Height2VolumetricDetachment(rd);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line, unless
		// program is running localy (in my computer)
		// Arguments should be:
		// args[0] - parameter file name
		// args[1] - epsDecayRate * pdpBulkConcentration
		// args[2] - pdpDiffusivity / pdpDecayRate
		boolean locationRemote = false;
		if (args.length > 0) {
			locationRemote = true;
			geometry = 3;
		}

		// set multigrid variables:
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new Gfp();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VaribleSeries thickness = new BiofilmMaximumThicknessSeries();
		if (!locationRemote) {
			// The following code will be omitted if no vizuals are desired
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// the biomass thickness visualizer
			app = new SeriesVizualizer(app, thickness);
			// add visualizer for solutes rates
			app = new SoluteRateSeriesVizualizer(app);
			// add visualizer for bulk concentrations
			app = new BulkConcentrationVizualizer(app);
			// add the controler
			app = new VizualModelControler(app);
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
			if (!locationRemote) {
				app.addStateWriter(new TimedStateWriterDecorator(
						new ParticlePositionWriter()));
				app.addStateWriter(new TimedStateWriterDecorator(
						new PovRayWriter()));
			}
			app.addStateWriter(new TimedStateWriterDecorator(
					new SoluteConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			// app.addTimedStateWriter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(Model.detachedBiomassContainer()
					.getErodedBiomassSeries());
			spw.addSeries(Model.detachedBiomassContainer()
					.getSloughedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			try {
				// initialize
				app.initializeDiffusionReactionSystem(); // also innoculates
			} catch (ModelException e) {
				e.printStackTrace();
				System.exit(-1);
			}
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
			// write the header in parameters.txt file
			spw.initializeParametersWriting();
			//
			Model.model().setMaximumTimeStep(outputEvery);
			Model.model().setFinishIterationTime(finishIteratinTime);
			// start iterating cycle
			// app.waitForStartIteratingRequest();
			app.startIterating();
		} catch (InterruptedException e1) {
			app.forceWriteState();
			e1.printStackTrace();
			System.out.println(e1);
		}
		System.out.println("Simulation finished.");
	}
}