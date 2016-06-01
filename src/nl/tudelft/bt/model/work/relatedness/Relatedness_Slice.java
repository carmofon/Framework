package nl.tudelft.bt.model.work.relatedness;

import java.awt.Color;
import java.util.Iterator;
import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BiofilmBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.FixedOnTopBoundaryLayer;
import nl.tudelft.bt.model.reaction.*;

/**
 * Competition between quorum-sensing EPS+ strain and non-QS EPS+ strain
 */
public class Relatedness_Slice extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	protected static String outputDirectory;

	// this variable switches between downregulation of EPS (if true)
	// and upregulation (if false)
	//protected static boolean downRegulation = true;

	// protected static String outputDirectory =
	// "C:/results/quorumsensing/test";

	protected static boolean graphics;
	protected static boolean gradientsOn;

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory containing important.
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
	//protected static float autoinducerBulkConcentration = 0.0f;
	//protected static float autoinducerDiffusivity = 1.0e4f;

	//
	// Particulate species (biomass X)
	protected static float specificMassX = 350f; // [g/L]

	// EPS
	//protected static float specificMassEPS = specificMassX / 6; // [g/L]

	// Yield coefficients
	// Biomass on substrate
	private static float Yxs = 0.5f; // [gX/gS]

	// Autoinducer on biomass
	//private static float YAX = 20f / 2.65e-3f; // [A/gX]

	// Processes
	// Max growth rate (Mourino-Perez et al., 2003)
	protected static float beta = 1.0f; // [1/h]

	// proportion of resources invested in EPS production
	protected static float f_QSpos = 0.5f; // [gEPS/gX]

	protected static float f_QSneg = f_QSpos; // [gEPS/gX]

	protected static float initialBiomassQSpos = 1;

	protected static float initialBiomassQSneg = 1;

	// substrate saturation constant
	//private static float KS = 3.5e-4f * 0.1f; // [g/L]

	// Autoinducer production constant; "sigma" changed to "aiProductionRate"
	//protected static float aiProductionRate = 0.2f; // [1/h]

	// Autoinducer quorum concentration threshold
	// "alpha" changed to "aiQuorumThreshold"
	//protected static float QSthreshold = 1f; // [1/L^2]

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
	protected static int initialParticleNumber_QSpos = 50;

	protected static int initialParticleNumber_QSneg = 50;

	// declare the eps particulate species out of defineSpeciesAndReactions
	// scope
	// so that we can access it from the detachment definition
	private ParticulateSpecies eps;

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

		// Autoinducer
		SoluteSpecies autoinducer = new SoluteSpecies("auotinducer",
				autoinducerDiffusivity);
		// autinducer also has constant bulk concentration
		autoinducer.setBulkConcentration(new ConstantBulkConcentration(
				autoinducerBulkConcentration));

		// 2. Create the particulate species (solids)
		// ***********************************************************************

		// Active mass of QS-positive strain; color will be overridden below
		ParticulateSpecies activeX_QSpos = new ParticulateSpecies(
				"activeX_QSpos", specificMassX, Color.yellow);

		// EPS (remember: eps has been declared outside of the current method)
		eps = new ParticulateSpecies("EPS", specificMassEPS, Color.yellow);

		// array of particulate species that constitute sp_QSpos
		ParticulateSpecies[] sp_QSpos = { activeX_QSpos, eps };
		// set initial composition of innoculating agents
		// float[] fractionalVolumeCompositionH1 = { initialBiomassQSpos,
		// 1-initialBiomassQSpos };
		float[] fractionalVolumeCompositionH1 = { 1f, 0f };

		// Active Mass of QS-negative strain
		ParticulateSpecies activeX_QSneg = new ParticulateSpecies(
				"activeX_QSneg", specificMassX, Color.blue);

		// If we want the QS-neg strain to have a unique EPS type . . .
		//
		// ParticulateSpecies eps_QSneg = new ParticulateSpecies("eps_QSneg",
		// specificMassEPS,
		// Color.yellow);

		// array of fixed species that constitute species_QSneg
		ParticulateSpecies[] ep_QSneg = { activeX_QSneg, eps };
		// float[] fractionalVolumeCompositionH2 = { initialBiomassQSneg,
		// 1-initialBiomassQSneg };
		float[] fractionalVolumeCompositionH2 = { 1f, 0f };

		// 3. Create the QS-positive and -negative biomass species
		// ***********************************************************************

		// QS-positive species, overrides activeX_QSpos color to reflect quorum
		// state
		BiomassSpecies species_QSpos = new QSBiomassSpecies("speciesQS",
				sp_QSpos, fractionalVolumeCompositionH1, autoinducer,
				QSthreshold, Color.green, Color.red);
		species_QSpos.setActiveMass(activeX_QSpos);
		species_QSpos.setEpsMass(eps);

		// QS-negative species
		BiomassSpecies species_QSneg = new BiomassSpecies("species_QSneg",
				ep_QSneg, fractionalVolumeCompositionH2);
		species_QSneg.setActiveMass(activeX_QSneg);
		species_QSneg.setEpsMass(eps);

		// 4. Create the Reaction factors, Monod and inhibition coefficients
		// ***********************************************************************

		// Saturation factor for bacterial growth
		ProcessFactor mS = new Saturation(substrate, KS);

		// Step function factor to simulate quorum sensing
		ProcessFactor quorumSensing;
		if (downRegulation)
			quorumSensing = new DownRegulationStep(autoinducer, QSthreshold);
		else
			quorumSensing = new UpRegulationStep(autoinducer, QSthreshold);

		// Constant factor to represent autoinducer production
		// ProcessFactor aiConst = new Constant(sigma);

		// 5. Create the reactions
		// ***********************************************************************

		// growth of QS-positive strain, saturating with substrate availability
		Reaction growth_QSpos = new Reaction("growth_QSpos", activeX_QSpos,
				uMax, 1);
		growth_QSpos.addFactor(mS);
		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass

		// EPS production for QS-positive strain, either full blast (below
		// quorum)
		// or off (above quorum)
		Reaction epsProduction = new Reaction("epsProduction", activeX_QSpos,
				uMax, 2);
		epsProduction.addFactor(mS);
		epsProduction.addFactor(quorumSensing);

		// Autoinducer production by both strains, constant
		// Reaction aiProduction = new Reaction("aiProduction", activeX, 1, 1);
		// aiProduction.addFactor(aiConst);
		Reaction aiProduction_QSpos = new Reaction("aiProduction_QSpos",
				activeX_QSpos, aiProductionRate, 0);
		Reaction aiProduction_QSneg = new Reaction("aiproduction_QSneg",
				activeX_QSneg, aiProductionRate, 0);

		// growth of QS-negative strain, saturating with substrate availability
		Reaction growth_QSneg = new Reaction("growth_QSneg", activeX_QSneg,
				uMax, 1);
		growth_QSneg.addFactor(mS);

		// EPS production for QS-negative strain, constant
		Reaction epsProduction2 = new Reaction("epsProduction2", activeX_QSneg,
				uMax, 1);
		epsProduction2.addFactor(mS);

		// 6. Assign reaction to the species through ReactionStoichiometries
		// ***********************************************************************

		// active mass of QS-positive strain
		NetReaction rsXactive_QSpos = new NetReaction(3);
		rsXactive_QSpos.addReaction(growth_QSpos, 1);
		rsXactive_QSpos.addReaction(epsProduction, -f_QSpos);
		rsXactive_QSpos.addReaction(aiProduction_QSpos, -(1 / YAX));
		activeX_QSpos.setProcesses(rsXactive_QSpos);

		// active mass of QS-negative strain
		NetReaction rsXactive_QSneg = new NetReaction(2);
		rsXactive_QSneg.addReaction(growth_QSneg, 1 - f_QSneg);
		rsXactive_QSneg.addReaction(aiProduction_QSneg, -(1 / YAX));
		activeX_QSneg.setProcesses(rsXactive_QSneg);

		// EPS
		NetReaction rsEPS = new NetReaction(2);
		rsEPS.addReaction(epsProduction, f_QSpos);
		rsEPS.addReaction(epsProduction2, f_QSneg);
		eps.setProcesses(rsEPS);

		// EPS, if we want QS-neg to have its own EPS type
		// NetReaction rsEPS2 = new NetReaction(1);
		// rsEPS2.addReaction(growth_QSneg, f_QSneg);
		// eps_QSneg.setProcesses(rsEPS2);

		// 7. Create Solute Reaction Stoichiometries
		// ***********************************************************************

		// Substrate - consumed at equal rates by both strains
		NetReaction rsSubstrate = new NetReaction(2);
		rsSubstrate.addReaction(growth_QSpos, -(1 / Yxs));
		rsSubstrate.addReaction(growth_QSneg, -(1 / Yxs));
		substrate.setProcesses(rsSubstrate);
		// This defines that substrate consumption rate is -(1 / YXS)*rX

		// Autoinducer - produced by both strains
		NetReaction rsAutoinducer = new NetReaction(2);
		rsAutoinducer.addReaction(aiProduction_QSpos, 1);
		rsAutoinducer.addReaction(aiProduction_QSneg, 1);
		autoinducer.setProcesses(rsAutoinducer);

		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(species_QSpos);
		addBiomassSpecies(species_QSneg);
		addSoluteSpecies(substrate);
		addSoluteSpecies(autoinducer);
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
		int[] nCells = { initialCooperator, initialCheater };
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new EPSonlyFractionDetachment(rdetach,
				eps, detachmentImpact);
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
		ApplicationComponent app = new Relatedness_Slice();
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
		if (args.length < 4)
			throw new RuntimeException(
					"First argument: output directory; Second argument: graphics on/off "
							+ "Third argument: AI production rate [0,1]; Fourth argument: seed");
		if (!args[0].contains("results"))
			throw new RuntimeException(
					"Output directory doesn't contain the word 'results'");
		outputDirectory = args[0];
		graphics = (Integer.parseInt(args[1]) == 0 ? false : true);
		float alpha_tilde = Float.parseFloat(args[2]);
		float L2 = 100 * 100;
		QSthreshold = alpha_tilde * L2 * aiProductionRate * specificMassX
				/ autoinducerDiffusivity;
		// seting the random number generator seed
		int seed = Integer.parseInt(args[3]);
		Model.model().setSeed(seed);
		outputEvery = 300;
		run();

	}

}