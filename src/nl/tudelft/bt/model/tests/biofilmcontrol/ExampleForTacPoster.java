package nl.tudelft.bt.model.tests.biofilmcontrol;

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
import nl.tudelft.bt.model.detachment.levelset.functions.Height2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;

/**
 * Implements a single species biofilm system where an EPS degrading chemical is
 * used to affect the growth of an EPS producing biofilm
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ExampleForTacPoster extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/jxavier/results/lixo";

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
	// protected static float qSmax = 0.952f; //[gS/gX/h]
	protected static float qSmax = 2f; // [gS/gX/h]

	// Substrate saturation constant
	private static float KS = 4e-3f; // [gS/L]

	// Oxygen saturation constant
	private static float KO2 = 3.5e-4f; // [gO/L]

	// Yield of biomass on consumed substrate
	private static float YX_S = 0.2f; // [gX/gS]

	// Yield of EPS on produced substrate consumed
	private static float YEPS_S = 0.29f; // [gEPS/gS]

	// diffusivity of solutes
	private static float substrateDiffusivity = 4.17e6f; // [um^2/h]

	private static float oxygenDiffusivity = 8.33e6f; // [um^2/h]

	// Enzyme parameters
	private static float epsDecayRate = 10000f; // [gEPS/gEnz/h]

	private static float enzymeDecayRate = 1000f; // [gEnz/gEnz/h]

	private static float enzymeDiffusivity = 2.1e6f; // [um^2/h]

	private static float applicationTime = 1000f; // [h]

	private static float beginApplicationAt = 800f; // [h]

	// Concentration of solutes
	protected static float substrateBulkConcentration = 100e-3f; // [gS/L]

	// protected static float oxygenBulkConcentration = 0.8e-3f; //[gO/L]
	protected static float oxygenBulkConcentration = 1e-3f; // [gO/L]

	protected static float enzymeBulkConcentration = 1e-3f; // [gEnz/L]

	// Specific mass of particulates
	// protected static float specificMassX = 24f; //[gX/L]
	// protected static float specificMassEPS = 4f; //[gEPS/L]
	protected static float specificMassX = 200f; // [gX/L]

	protected static float specificMassEPS = specificMassX / 6f; // [gEPS/L]

	// detachment
	protected static float rdetach = 1e-3f;

	// protected static float rdetach = 0;

	// uncomment
	// inoculation
	protected static int initialParticleNumber = 100;

	// Numeric parameters
	protected static float systemSize = 2000; // [um]

	// grid resolution
	protected static int gridSide = 65; // multigrid grid side

	// comment
	// protected static int initialParticleNumber = 20;
	// protected static float systemSize = 500; // [um]
	// protected static int gridSide = 17; // multigrid grid side
	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 6 / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	protected static float relativeBoundaryLayer = 100 / systemSize;

	// Shoving parameter
	protected static float kShov = 1.0f;

	// outpute (write results to file) every:
	protected static float outputEvery = 15.0f; // [h]

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// //1. Create the solutes
		// substrate (S)
		// SoluteSpecies substrate = new SoluteSpecies("substrate",
		// substrateDiffusivity);
		// substrate.setBulkConcentration(new ConstantBulkConcentration(
		// substrateBulkConcentration));
		// substrate (O)
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		SoluteSpecies enzyme = new SoluteSpecies("enzyme", enzymeDiffusivity);
		enzyme.setBulkConcentration(new SinglePulseBulkConcentration(
				enzymeBulkConcentration, beginApplicationAt, applicationTime));
		// enzyme.setBulkConcentration(new ConstantBulkConcentration(0));
		// 2. Create the particulate species (soliids)
		// X active mass
		ParticulateSpecies activeX = new ParticulateSpecies("activeX",
				specificMassX, Color.red);
		// EPS
		ParticulateSpecies eps = new ParticulateSpecies("EPS", specificMassEPS,
				Color.yellow);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spX = { activeX, eps };
		float[] fractionalVolumeCompositionH1 = { 1, 0 };
		// 3. Create the biomass species
		BiomassSpecies speciesX = new BiomassSpecies("speciesX", spX,
				fractionalVolumeCompositionH1);
		speciesX.setActiveMass(activeX);
		speciesX.setEpsMass(eps);
		// speciesX.getColorFromGrowth();
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// ProcessFactor mS = new Saturation(substrate, KS);
		ProcessFactor mO2 = new Saturation(oxygen, KO2);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		//
		// 5. Create the reactions
		// growth
		Reaction growth = new Reaction("growth", activeX, qSmax, 1);
		// growth.addFactor(mS);
		growth.addFactor(mO2);
		// EPS decay
		Reaction epsDecay = new Reaction("epsDecay", enzyme, epsDecayRate, 0);
		// Enzyme decay
		Reaction enzymeDecay = new Reaction("enzymeDecay", enzyme,
				enzymeDecayRate, 0);
		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactive = new NetReaction(1);
		rsXactive.addReaction(growth, YX_S);
		activeX.setProcesses(rsXactive);
		// EPS
		NetReaction rsEps = new NetReaction(2);
		rsEps.addReaction(growth, YEPS_S);
		rsEps.addReaction(epsDecay, -1);
		eps.setProcesses(rsEps);
		//
		// assign reaction stoichiometry to the solutes
		// substrate
		// NetReaction rsSubstrate = new NetReaction(1);
		// rsSubstrate.addReaction(growth, -(1 / YX));
		// substrate.setProcesses(rsSubstrate);
		// oxygen
		NetReaction rsOxygen = new NetReaction(1);
		rsOxygen.addReaction(growth, -(1 - YX_S - YEPS_S));
		oxygen.setProcesses(rsOxygen);
		// enzyme
		NetReaction rsEnzyme = new NetReaction(1);
		rsEnzyme.addReaction(enzymeDecay, -1);
		enzyme.setProcesses(rsEnzyme);
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesX);
		// addSoluteSpecies(substrate);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(enzyme);
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
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		// DetachmentSpeedFunction df = new Height2EpsDetachment(rdetach,
		// specificMassEPS);
		DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line, unless
		// program is running localy (in my laptop)
		boolean locationRemote = false;
		if (args.length > 0) {
			rdetach = Float.parseFloat(args[0]);
			outputDirectory = "results/rdetach_" + args[0];
			locationRemote = true;
			System.out.println("simulation started with rdetach = " + args[0]);
		}
		// set multigrid variables:
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new ExampleForTacPoster();
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
			app
					.addStateWriter(new TimedStateWriterDecorator(
							new PovRayWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SoluteConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new ParticlePositionWriter()));
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
			// write the header in parameters.txt file
			spw.initializeParametersWriting();
			//
			Model.model().setCompulsoryTimeStep(outputEvery);
			// start iterating cycle
			app.startIterating();
		} catch (Exception e1) {
			app.forceWriteState();
			e1.printStackTrace();
			System.out.println(e1);
		}
		System.out.println("Simulation finished.");
	}
}