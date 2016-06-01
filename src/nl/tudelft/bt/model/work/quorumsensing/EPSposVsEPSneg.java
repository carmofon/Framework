package nl.tudelft.bt.model.work.quorumsensing;

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
import nl.tudelft.bt.model.multigrid.boundary_conditions.BiofilmBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.FixedOnTopBoundaryLayer;
import nl.tudelft.bt.model.reaction.*;

/**
 * Competition between EPS+ strain and EPS- strain. Both strais will be
 * producing EPS, but at different investments. This model is used for ESS
 * determination of optimal EPS investment
 */
public class EPSposVsEPSneg extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	// protected static String outputDirectory =
	// "/users/careynadell/BiofilmResults/practice";
	protected static String outputDirectory; // =

	// "/Volumes/lsdivfs2.cgr.harvard.edu/users/cnadell/results/practice";

	// protected static String outputDirectory =
	// "C:/results/quorumsensing/test";

	protected static boolean graphics;

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

	// Solute species
	// Substrate (S)
	protected static float substrateBulkConcentration = 5e-4f; // [g/L]

	private static float substrateDiffusivity = 5.76e6f; // [um^2/h]

	// Autoinducer (A)
	protected static float autoinducerBulkConcentration = 0.0f;

	protected static float autoinducerDiffusivity = 1.0e4f;

	//
	// Particulate species (biomass X)
	protected static float specificMassX = 350f; // [g/L]

	// EPS
	protected static float specificMassEPS = specificMassX / 6; // [g/L]

	// Yield coefficients
	// Biomass on substrate
	private static float YXS = 0.5f; // [gX/gS]

	// Processes
	// Max growth rate (Mourino-Perez et al., 2003)
	protected static float uMax = 1.0f; // [1/h]

	// proportion of resources invested in EPS production
	protected static float f_EPSpos = 0.5f; // [gEPS/gX]

	protected static float f_EPSneg = f_EPSpos; // [gEPS/gX]

	// substrate saturation constant
	private static float KS = 3.5e-4f * 0.1f; // [g/L]

	// Computation parameters
	protected static float systemSize = 250; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.003f; // 0.003

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0002f; // 0.0001

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 0.4f;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// protected static float rdetach = 3e-5f;
	// protected static float rdetach = 1e-3f; // 1e-6
	protected static float rdetach = 0f; // Joao

	protected static float detachmentImpact = 0f;

	// Initial particle number split into separate values for quorum sensing and
	// non-quorum sensing strains
	protected static int initialParticleNumber_EPSpos = 50;

	protected static int initialParticleNumber_EPSneg = 50;

	static float outputEvery = 350.0f;

	static float simulationFinishTime = 340.0f;

	/**
	 * Define QS and EPS-only species; specify solutes, particulates, and
	 * reactions
	 */
	private void defineSpeciesAndReactions() throws ModelException {

		// 1. Create the solutes
		// ***********************************************************************

		// Substrate
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		// set up the simplest type of bulk concentration: constant
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));

		// 2. Create the particulate species (solids)
		// ***********************************************************************

		// Active mass of EPS-positive strain
		ParticulateSpecies activeX_EPSpos = new ParticulateSpecies(
				"activeX_EPSpos", specificMassX, Color.red);

		// EPS
		ParticulateSpecies eps = new ParticulateSpecies("EPS", specificMassEPS,
				Color.yellow);

		// array of particulate species that constitute sp_EPSpos
		ParticulateSpecies[] sp_EPSpos = { activeX_EPSpos, eps };
		// set initial composition of innoculating agents
		float[] fractionalVolumeCompositionH1 = { 1f, 0f };

		// Active Mass of EPS-negative strain
		ParticulateSpecies activeX_EPSneg = new ParticulateSpecies(
				"activeX_EPSneg", specificMassX, Color.blue);

		// array of fixed species that constitute species_QSneg
		ParticulateSpecies[] sp_EPSneg = { activeX_EPSneg, eps };
		float[] fractionalVolumeCompositionH2 = { 1f, 0f };

		// 3. Create the QS-positive and -negative biomass species
		// ***********************************************************************

		// QS-positive species, overrides activeX_QSpos color to reflect quorum
		// state
		BiomassSpecies species_EPSpos = new BiomassSpecies("species_EPSpos",
				sp_EPSpos, fractionalVolumeCompositionH1);
		species_EPSpos.setActiveMass(activeX_EPSpos);
		species_EPSpos.setEpsMass(eps);

		// QS-negative species
		BiomassSpecies species_EPSneg = new BiomassSpecies("species_EPSpos",
				sp_EPSneg, fractionalVolumeCompositionH2);
		species_EPSneg.setActiveMass(activeX_EPSneg);
		species_EPSneg.setEpsMass(eps);

		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// ***********************************************************************

		// Saturation factor for bacterial growth
		ProcessFactor mS = new Saturation(substrate, KS);

		// 5. Create the reactions
		// ***********************************************************************

		// growth of QS-positive strain, saturating with substrate availability
		Reaction growth_EPSpos = new Reaction("growth_EPSpos", activeX_EPSpos,
				uMax, 1);
		growth_EPSpos.addFactor(mS);
		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass

		// growth of QS-negative strain, saturating with substrate availability
		Reaction growth_EPSneg = new Reaction("growth_EPSneg", activeX_EPSneg,
				uMax, 1);
		growth_EPSneg.addFactor(mS);

		// 6. Assign reaction to the species through ReactionStoichiometries
		// ***********************************************************************

		// active mass of EPSpos strain
		NetReaction rsXactive_EPSpos = new NetReaction(1);
		rsXactive_EPSpos.addReaction(growth_EPSpos, 1-f_EPSpos);
		activeX_EPSpos.setProcesses(rsXactive_EPSpos);

		// active mass of EPSneg strain
		NetReaction rsXactive_QSneg = new NetReaction(1);
		rsXactive_QSneg.addReaction(growth_EPSneg, 1 - f_EPSneg);
		activeX_EPSneg.setProcesses(rsXactive_QSneg);

		// EPS
		NetReaction rsEPS = new NetReaction(2);
		rsEPS.addReaction(growth_EPSpos, f_EPSpos);
		rsEPS.addReaction(growth_EPSneg, f_EPSneg);
		eps.setProcesses(rsEPS);

		// EPS, if we want QS-neg to have its own EPS type
		// NetReaction rsEPS2 = new NetReaction(1);
		// rsEPS2.addReaction(growth_QSneg, f_QSneg);
		// eps_QSneg.setProcesses(rsEPS2);

		// 7. Create Solute Reaction Stoichiometries
		// ***********************************************************************

		// Substrate - consumed at equal rates by both strains
		NetReaction rsSubstrate = new NetReaction(2);
		rsSubstrate.addReaction(growth_EPSpos, -(1 / YXS));
		rsSubstrate.addReaction(growth_EPSneg, -(1 / YXS));
		substrate.setProcesses(rsSubstrate);
		// This defines that substrate consumption rate is -(1 / YXS)*rX

		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(species_EPSpos);
		addBiomassSpecies(species_EPSneg);
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
		int[] nCells = { initialParticleNumber_EPSpos,
				initialParticleNumber_EPSneg };
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rate of detachment = 0
		DetachmentSpeedFunction df = new Height2MassDetachment(0);
		setDetachmentHandler(df);
	}

	protected void createBoundaryLayer(float h)
			throws MultigridSystemNotSetException {
		// uncomment the type of boundary layer wanted
		// _boundaryLayer = new FrontBoundaryLayer(h);
		_boundaryLayer = new FixedOnTopBoundaryLayer();
		// create the boundary conditions
		MultigridVariable.setBoundaryConditions(new BiofilmBoundaryConditions());
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void run() {
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new EPSposVsEPSneg();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		// detached biomass
		VariableSeries detached = Model.model().detachedBiomassContainer()
				.getTotalDetachedBiomassSeries();
		// The following code will be omitted if no vizuals are desired
		// start decorationg the application
		// /*
		if (graphics) {
			app = new BiomassVizualizer(app);
			// the biomass thickness visualizer
			app = new SeriesVizualizer(app, thickness);
			app = new SeriesVizualizer(app, detached); // */
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
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
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
				spw.addSeries(new DetachedParticulateSpecies(s));
			}
			// start iterating cycle
			Model.model().setCompulsoryTimeStep(outputEvery);
			Model.model().setFinishIterationTime(simulationFinishTime);
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
		}
		try {
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		System.out.println("Simulation finished.");
	}

	public static void main(String[] args) {
		if (args.length < 3)
			throw new RuntimeException(
					"First argument: output directory; Second argument: graphics on/off "
							+ "Third argument: seed");
		if (!args[0].contains("results"))
			throw new RuntimeException(
					"Output directory doesn't contain the word 'results'");
		outputDirectory = args[0];
		graphics = (Integer.parseInt(args[1]) == 0 ? false : true);
		// seting the random number generator seed
		int seed = Integer.parseInt(args[2]);
		Model.model().setSeed(seed);
		outputEvery = 300;
		run();
	}

}