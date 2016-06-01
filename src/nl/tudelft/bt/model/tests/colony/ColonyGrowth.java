package nl.tudelft.bt.model.tests.colony;

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
 * Simulates the growth of a colony to demonstrate individual based modelling in
 * presentations
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ColonyGrowth extends GranuleModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory = "C:joao/results/lixo";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	// //Reaction parameters
	// //ammonium Oxydizers (XNH)
	protected static float uMaxXNH = 0.0854f; // [gCOD-XNH/gCOD-XNH/h]

	// private static float K_XNH_O2 = 0.6e-3f; //[gO2/L]
	// private static float K_XNH_NH4 = 2.4e-3f; //[gN/L]

	private static float Y_XNH_NH4 = 0.150f; // [gCOD-XNH/gN]

	// //Nitrite oxydizers (XNO)
	// protected static float uMaxNO = 0.0604f; //[gCOD-XNO/gCOD-XNO/h]
	// private static float K_XNO_O2 = 2.2e-3f; //[gO2/L]
	// private static float K_XNO_NO2 = 5.5e-3f; //[gN/L]
	// private static float Y_XNO_NO2 = 0.041f; // [gCOD-XNO/gN]
	// //Heterotrophs (XH)
	// protected static float qSMaxXH = 0.952f; //[gCOD-XH/gCOD-XH/h]
	// private static float K_XH_O2 = 3.5e-4f; //[gO2/l]
	// private static float K_XH_S = 4e-3f; //[gCOD/l]
	// private static float Y_H_S = 0.495f; //[gCOD-XH/gCOD-S]
	// private static float fPHBmaxXH = 0.33f; //[gCOD-PHB/gCOD-XH]
	// private static float Y_PHB_S_XH = 0.667f; //[gCOD-PHB/gCOD-S]
	// protected static float kPhbXH = 0.15f; //[gCOH-PHB/gCOD-PHB/h]
	// private static float Y_XH_PHB = 0.668f; //[gCOD-H/gCOD-PHB]
	// //Phospate accumulating organisms (XPAO)
	// protected static float qSMaxPAO = 0.4029f; //[gCOD-S/gCOD-XPAO/h]
	// private static float KPAOO2 = 7e-4f; //[gO2/L]
	// private static float KPAOS = 4e-3f; //[gCOD/L]
	// private static float KPAOPP = 0.01e-3f; //[gP/L]
	// private static float KPAOGLY = 0.01e-3f; //[gCOD-Gly/L]
	// private static float Y_PHB_S_of_PAO = 1.5f; //[gCOD-PHB/gCOD-S]
	// private static float Y_PO4_S = 0.3f; //[gP/gCOD-S]
	// private static float fPHBmaxPAO = 1.0f; //[gCOD-PHB/gCOD-XPAO]
	// //protected static float kPhbPAO = 0.3f; //[gCOH-PHB/gCOD-PHB/h]
	// protected static float kPhbPAO = 0.15f; //[gCOH-PHB/gCOD-PHB/h]
	// private static float YPhbAerobicPAO = 1.39f; //[gCOD-PHB/gCOD-XPAO]
	// //private static float YPhbAerobicPAO = 1.5f; //[gCOD-PHB/gCOD-XPAO]
	// private static float YPhbAnoxicPAO = 1.7f; //[gCOD-PHB/gCOD-XPAO]
	// private static float KPAOPO4 = 3.1e-3f; //[gP/L]
	// protected static float kPPPAO = 0.0187f; //[gP/gCOD-PHB/h]
	// private static float YPpAerobicPAO = 4.37f; //[gP/gCOD-PAO]
	// private static float YPpAnoxicPAO = 2.92f; //[gP/gCOD]
	// private static float fPPmaxPAO = 0.5f; //[gP/gCOD-XPAO]
	// private static float Y_GLY_PHB = 0.73f; //[gCOD-Gly/gCOD-PHB]
	// protected static float kGlyPAO = 0.0454f; //[gCOD-Gly/gCOD-PHB/h]
	// private static float YGlyAerobicPAO = 1.12f; //[gCOD-Gly/gCOD-XPAO]
	// private static float YGlyAnoxicPAO = 1.15f; //[gCOD-Gly/gCOD]
	// private static float fGLYmaxPAO = 0.45f; //[gCOD-Gly/gCOD-PAO]
	// //Common to all species
	// protected static float kDecay = 0.0033f; //[gCOD/gCOD/h]
	// private static float fIDecay = 0.4f; //[gCOD-I/gCOD-X]
	// private static float iNBM = 0.07f; //[gN/gCOD-X]
	// private static float iNXI = 0.02f; //[gN/gCOD-I]
	// private static float KNH4 = 0.05e-3f; //[gN/L]
	// private static float KNO3 = 1e-3f; //[gN/L]
	// private static float eta = 0.2f;

	// geometry (default is 3D) - change value to 2 for 2D
	protected static int geometry = 2;

	// //Reactor properties
	// protected static float reactorVolume = 3; // [L]
	// protected static float feedFraction = 0.5f; // [dimensionless]
	//
	// //duration of the SBR cycle
	// protected static float cycleTime = 3; //[h]
	// //The SBR sycle
	// private static SbrBulkConcentration.SbrCycle _cycle = new
	// SbrBulkConcentration.SbrCycle(
	// cycleTime);
	// private static int _nCycles = 400; // write full data concerning cycle
	// every
	// // _nCycles
	//
	// //Feed composition
	protected static float oxygenBulkConcentration = 10e-3f; // [gO2/L] DO

	// 100
	// //protected static float oxygenBulkConcentration = 2e-3f; //[gO2/L] DO 20
	//
	// protected static float ammoniumFeedConcentration = 71e-3f; //[gN/L]
	//
	// protected static float nitriteFeedConcentration = 0; //[gN/L]
	//
	// protected static float nitrateFeedConcentration = 0; //[gN/L]
	//
	// protected static float substrateFeedConcentration = 396e-3f; //[gCOD/L]
	//
	// protected static float phosphateBulkConcentration = 20e-3f; //[gP/L]
	//

	private static float oxygenDiffusivity = 8.3e6f; // [um^2/h]

	//
	// private static float substrateDiffusivity = 4e6f; //[um2/h]
	//
	// private static float phosphateDiffusivity = 4e6f; //[um2/h] ASSUMED!
	//
	// private static float ammoniumDiffusivity = 7.083e6f; //[um^2/h]
	//
	// private static float nitrateDiffusivity = 6.6667e6f; //[um^2/h]
	//
	// private static float nitriteDiffusivity = 6.6667e6f; //[um^2/h]
	//
	// //Leave these values
	// protected static float nConcentrationsPrecision = 0.3e-3f; //[gN/L]
	//
	// protected static float sConcentrationsPrecision = 4e-3f; //[gCOD/L]
	//
	// protected static float pConcentrationsPrecision = 0.15e-3f; //[gP/L]

	//
	// Particulate species (biomass H)
	protected static float specificMassBiomass = 150f; // [gCOD-H/L]

	// protected static float specificMassPolymers = specificMassBiomass * 1e5f;
	// // [gCOD-PHB/L]

	// Computation parameters
	// Size of computational volume (size of size of square)
	protected static float systemSize = 1600; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.007f;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 10 um
	protected static float relativeBoundaryLayer = 10 / systemSize;

	protected static float maximumGranuleRadius = 550; // [um]

	// other model parameters
	protected static int gridSide = 33; // multigrid grid side

	// Don't change this

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// Don't change this

	// detachment rate
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]

	// leave at zero to form round granules

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumberXNH = 1;

	// protected static int initialParticleNumberXNO = 20;
	//
	// protected static int initialParticleNumberXH = 20;
	//
	// protected static int initialParticleNumberXPAO = 20;

	// iteration finish time
	protected static float simulationFinishTime = 80f; // [h]

	// outpute (write results to file) every:
	protected static float outputEvery = 0.3f; // [h]

	// Computational volume multiplier
	// protected static float nComp = 25.93f * 6.7e5f; // [dimensionless]

	// factor 25.93f makes substrate loading per biomass the same as in real
	// reactor

	// /END OF PARAMETERS

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		// 2. Create the particulate species (solids)
		// NH
		ParticulateSpecies activeNH = new ParticulateSpecies("activeNH",
				specificMassBiomass, Color.yellow);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		ParticulateSpecies[] spNH = { activeNH };
		float[] fractionalVolumeCompositionNH = { 1.0f };
		// 3. Create the biomass species
		BiomassSpecies speciesNH = new BiomassSpecies("speciesNH", spNH,
				fractionalVolumeCompositionNH);
		speciesNH.setActiveMass(activeNH);
		// 5. Create the reactions
		// NH
		// growth NH
		Reaction growthNH = new Reaction("growthNH", activeNH, uMaxXNH, 0);
		//
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass NH
		NetReaction rsNHactive = new NetReaction(1);
		rsNHactive.addReaction(growthNH, 1);
		activeNH.setProcesses(rsNHactive);
		// assign reaction stoichiometry to the solutes
		// oxygen
		NetReaction rsOxygen = new NetReaction(1);
		rsOxygen.addReaction(growthNH, 1 - 3.43f / Y_XNH_NH4);
		oxygen.setProcesses(rsOxygen);
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesNH);
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
		int[] nCells = { initialParticleNumberXNH };
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
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line, unless
		// program is running localy (in my laptop)
		boolean locationRemote = false;
		if (args.length > 0) {
			outputDirectory = args[0];
			locationRemote = true;
		}
		// set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new ColonyGrowth();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the biovolume series
		VaribleSeries biovolume = new BiovolumeSeries();
		VaribleSeries[] runLengthSeries = { new RunLengthXSeries(),
				new RunLengthYSeries(), new RunLengthZSeries() };
		// The following code will be omitted if no vizuals are desired
		if (!locationRemote) {
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// add vizualizer for solutes rates
			// app = new SoluteRateSeriesVizualizer(app);
			// the biovolume visualizer
			// app = new SeriesVizualizer(app, biovolume);
			// runlength
			// app = new SeriesVizualizer(app, runLengthSeries, "Run Length");
			// detached biomass
			// app = new DetachedBiomassVizualizer(app);
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
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			// Pov witer is added twice
			app
					.addStateWriter(new TimedStateWriterDecorator(
							new PovRayWriter()));
			// app
			// .addStateWriter(new TimedStateWriterDecorator(
			// new PovRayWriter()));
			// app.addStateWriter(new TimedStateWriterDecorator(
			// new SoluteConcentrationWriter()));
			// app.addStateWriter(new TimedStateWriterDecorator(
			// new SolidsConcentrationWriter()));
			// app.addStateWriter(new TimedStateWriterDecorator(
			// new SolidsConcentrationWriter()));
			// app.addStateWriter(new TimedStateWriterDecorator(
			// new ParticlePositionWriter()));
			// app.addStateWritter(new DetachmentLevelSetWriter());
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
			Model.model().setCompulsoryTimeStep(outputEvery);
			Model.model().setFinishIterationTime(simulationFinishTime);
			// app.waitForStartIteratingRequest();
			app.startIterating();
		} catch (Exception e1) {
			app.forceWriteState();
			e1.printStackTrace();
			System.out.println(e1);
		}
		System.out.println("Simulation finished.");
	}
}