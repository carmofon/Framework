package nl.tudelft.bt.model.work.carlos;

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
import nl.tudelft.bt.model.multigrid.boundary_conditions.BiofilmBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.FixedOnTopBoundaryLayer;
import nl.tudelft.bt.model.multigrid.boundary_layers.SphericalDilationBoundaryLayer;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.particlebased.tumor.TumorModelHandler;
import nl.tudelft.bt.model.reaction.*;

/**
 * @author Carlos Carmona Fontaine (carmonac@mskcc.org)
 */
public class ModelWithAngiogenesis extends ModelHandler {
	// public class Example1KWithTumorHandler extends TumorModelHandler {
	// All the model parameters are defined here as static attributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/CarmonaFontaine/Documents/Xlab/Models/Angiogenesis/Test";

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

	/**
	 * Constants
	 */

	// /////////// Computation parameters

	protected static float systemSize = 500; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 5f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.7f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 0f;

	// other model parameters
	protected static int gridSide = 33; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less] NEED
											// THIS?

	protected static float rdetach = 1e-10f; // NEED THIS?

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 100;

	// /////////// Biochemical parameters //

	protected static SoluteSpecies oxygen;
	protected static SoluteSpecies lactate;
	protected static SoluteSpecies glucose;
	protected static SoluteSpecies vegf;

	protected static float oxygenBulkConcentration = 0.3f; // [g/L]Joao's paper 
															// Should be 16-22 ml O2/100ml arterial blood
															// and 13-16 ml in venous blood ()
	protected static float lactateBulkConcentration = 0; // [g/L]
	protected static float glucoseBulkConcentration = 1f; //1g/L blood plasma
	protected static float vegfBulkConcentration = 0f; 


	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	// Glucose diffusivity (Casciari et al 1988) 1.1e-6 cm2/s = 3.96e5um2/h
	// mouse; 5.5e-7-2.3e-7 cm2/s = 1.98e4um2/h to 8.2e3 um2/h
	//private static float oxygenDiffusivity = 7.2e7f; // 10x Diffusivity
	private static float oxygenDiffusivity = 7.2e6f;//[um^2/h]Thomlinson: 2e-5
														// cm2/s
	private static float lactateDiffusivity = 1.04e6f; // [um^2/h] 2.9e-6 cm2/s.
														// Measurement of
														// Effective
														// Diffusivities of
														// Lactose and Lactic
														// Acid in 3% Agarose
														// Gel Membrane
														// Amarjeet S. Bassi,
														// Sohrab Rohani and D.
														// G. Macdonaid
	private static float glucoseDiffusivity = 1.04e6f; // CHECK
	private static float vegfDiffusivity = 1.04e6f; // CHECK

	// protected static SoluteSpecies Hydrogen;

	/**
	 * Variables
	 */

	// /////////// Biochemical parameters

	// Yield coefficients
	//private static float gYo2 = 0.8f; // [gX/gS] Glucose Yield Aerobic
	private static float gYana = 0.5f; // [gX/gS] Lactate Yield 
	private static float oY = 0.5f; // [gX/gS] Oxygen Yield
	private static float lY = 1f - gYana; // [gX/gS] Lactate production Yield
	private static float clY = lY * 0.01f; // [gX/gS] Lactate consumption Yield
//	private static float switchValue = lY * 0.000000000001f;
	private static float switchValue = oY * 0.25f;

	
	// Satutration constants

	private static float oKST = oxygenBulkConcentration / 2f; // [g/L]
	private static float oKSM = oxygenBulkConcentration / 20f; // [g/L]

	private static float lKST = 1e-2f; // [g/L] saturation constant lactate
									// effect on tumor cells
	private static float lKSM = 1e-2f; // [g/L] saturation constant lactate
										// effect on macrophages

	// /////////// Cell parameters

	protected static int numberofCellGroups = 2;
	//
	// Particulate species (biomass X)
	protected static float specificMassT = 100f; // [g/L]
	protected static float specificMassM = 100f; // [g/L]
	protected static float specificMassE = 100f; // [g/L]

	// Max growth rates
	protected static float uMaxT = 0.01f; // [1/h]
	protected static float uMaxM = uMaxT * 1f; // [1/h]

	// Death rates
	protected static float kdT = uMaxT * 4f;// lactate effect strength on tumor cells (relative to growth)
	protected static float kdM = kdT * 3f;// lactate effect strength on macrophages (relative to cancer sensitivity)

	
	/**
	 * Definitions and Processes
	 */

	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// 1. Create the solutes
		// Oxygen
		oxygen = new SoluteSpecies("Oxygen", oxygenDiffusivity);
		// set up the simplest type of bulk concentration: constant
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		// Glucose
		glucose = new SoluteSpecies("Glucose", glucoseDiffusivity);
		// set up the simplest type of bulk concentration: constant
		glucose.setBulkConcentration(new ConstantBulkConcentration(
				glucoseBulkConcentration));
		
		// Lactate
		lactate = new SoluteSpecies("Lactate", lactateDiffusivity);
		// set up the simplest type of bulk concentration: constant
		lactate.setBulkConcentration(new ConstantBulkConcentration(
				lactateBulkConcentration));
		// VEGF
		vegf = new SoluteSpecies("VEGF", vegfDiffusivity);
		// set up the simplest type of bulk concentration: constant
		vegf.setBulkConcentration(new ConstantBulkConcentration(
				vegfBulkConcentration));

		// 2. Create the particulate species (solids)

		// X active mass
		ParticulateSpecies activeT = new ParticulateSpecies("activeT",
				specificMassT, Color.blue);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spT = { activeT };
		float[] fractionalVolumeCompositionT = { 1.0f };

		ParticulateSpecies activeM = new ParticulateSpecies("activeM",
				specificMassM, Color.green);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spM = { activeM };
		float[] fractionalVolumeCompositionM = { 1.0f };
		
		ParticulateSpecies activeE = new ParticulateSpecies("activeE",
				specificMassE, Color.red);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spE = { activeE };
		float[] fractionalVolumeCompositionE = { 1.0f };

		// 3. Create the biomass species
		BiomassSpecies speciesT = new BiomassSpecies("speciesT", spT,
				fractionalVolumeCompositionT);
		speciesT.setActiveMass(activeT);
		// speciesCC.getColorFromGrowth();

		BiomassSpecies speciesM = new BiomassSpecies("speciesM", spM,
				fractionalVolumeCompositionM);
		speciesM.setActiveMass(activeM);
		speciesM.setInducibleColor(oxygen, switchValue, Color.green,
				new Color(1f, 0, 1f));//HOW MAKE DEPENDENT OF LACTATE + OXYGEN
		// speciesTAM.getColorFromGrowth();
		
		BiomassSpecies speciesE = new BiomassSpecies("speciesE", spE,
				fractionalVolumeCompositionE);
		speciesE.setActiveMass(activeE);

		// 4. Create the Reaction factors, Monod and inhibition coefficients

		// Growth factor production


		// ProcessFactor mSo = new Saturation(glucose, oKS);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		// growth
		// aerobic growth weight
		ProcessFactor TumorAeroWeight = new Saturation(oxygen, oKST);
		// anaerobic growth weight
		ProcessFactor TumorAnaWeight = new Inhibition(oxygen, oKST);

		// aerobic growth weight
		ProcessFactor MacroAeroWeight = new Saturation(oxygen, oKSM);
		// anaerobic growth weight
		ProcessFactor MacroAnaWeight = new Inhibition(oxygen, oKSM);
		// lactate growth weight
		ProcessFactor MacroLactWeight = new Saturation(lactate, lKSM);

		// ProcessFactor Endogenous = new Saturation(lactate, lKST);
		ProcessFactor LactateEffT = new Saturation(lactate, lKST);
		ProcessFactor LactateEffM = new Saturation(lactate, lKSM);

		// 5. Create the reactions
		// Aerobic growth
		Reaction aeroGrowthT = new Reaction("aeroGrowthT", activeT, uMaxT, 1);
		aeroGrowthT.addFactor(TumorAeroWeight);

		Reaction aeroGrowthM = new Reaction("aeroGrowthM", activeM, uMaxM, 1);
		aeroGrowthM.addFactor(MacroAeroWeight);

		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass

		// Anaerobic growth
		Reaction anaGrowthT = new Reaction("anaGrowthT", activeT, uMaxT, 1);
		anaGrowthT.addFactor(TumorAnaWeight);

		Reaction anaGrowthM = new Reaction("anaGrowthM", activeM, uMaxM, 1);
		anaGrowthM.addFactor(MacroAnaWeight);

		// Endogenous Decay and death

		Reaction LactateKillingT = new Reaction("LactateKillingT", activeT,
				kdT, 1);
		LactateKillingT.addFactor(LactateEffT);

		Reaction LactateKillingM = new Reaction("LactateKillingM", activeM,
				kdM, 1);
		LactateKillingM.addFactor(LactateEffM);

		// Growth due to lactate consumption		
		
		Reaction lactGrowthM = new Reaction("lactGrowthM", activeM, uMaxM, 1);
		lactGrowthM.addFactor(MacroLactWeight);
		
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactiveT = new NetReaction(3);
		rsXactiveT.addReaction(aeroGrowthT, 1);
		rsXactiveT.addReaction(anaGrowthT, 1);
		rsXactiveT.addReaction(LactateKillingT, -1);
		activeT.setProcesses(rsXactiveT);

		NetReaction rsXactiveM = new NetReaction(4);
		rsXactiveM.addReaction(aeroGrowthM, 1);
		rsXactiveM.addReaction(anaGrowthM, 1);
		rsXactiveM.addReaction(lactGrowthM, 1);
		rsXactiveM.addReaction(LactateKillingM, -1);
		activeM.setProcesses(rsXactiveM);
		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// substrate

		NetReaction rsOxygen = new NetReaction(2);
		rsOxygen.addReaction(aeroGrowthT, -(1 / oY));
		rsOxygen.addReaction(aeroGrowthM, -(1 / oY));
		oxygen.setProcesses(rsOxygen);

		NetReaction rsLactate = new NetReaction(3);
		rsLactate.addReaction(anaGrowthT, (lY));
		rsLactate.addReaction(anaGrowthM, (lY));
		rsLactate.addReaction(aeroGrowthM,-(1 / clY));
		lactate.setProcesses(rsLactate);
		
//		NetReaction rsIschemia = new NetReaction(2);
//		rsIschemia.addReaction(rsLactate, (lY));
//		rsIschemia.addReaction(anaGrowthM, (lY));
//		ischemia.setProcesses(rsIschemia);

		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesT);
		addBiomassSpecies(speciesM);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(lactate);
	//	addSoluteSpecies(ischemia);

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
		int[] nCells = { initialParticleNumber, initialParticleNumber };
		// inoculateRandomlyInsideRadius(nCells, initialColonyRadius);
	//	inoculateRandomlyEverywhere(nCells);
		inoculateRandomly(nCells);

	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Height2MassDetachment(0);
		setDetachmentHandler(df);
	}

	protected void createBoundaryLayer(float h)
			throws MultigridSystemNotSetException {
		_boundaryLayer = new FixedOnTopBoundaryLayer();
		// create the boundary conditions
		MultigridVariable
				.setBoundaryConditions(new BiofilmBoundaryConditions());
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new ModelWithAngiogenesis();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// The following code will be omitted if no vizuals are desired
		// start decorationg the application
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		//app = new SeriesVizualizer(app, biomass);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// --- nothing to set in this case: constant bulk concentration
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			app.initializeDiffusionReactionSystem(); // also innoculates

			// grafics
			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new ImageWriter());
			app.addTimedStateWriter(new ImageWriter(lactate));
			app.addTimedStateWriter(new ImageWriter(oxygen));
			// app.addTimedStateWriter(new ImageWriter(egf));
			// app.addTimedStateWriter(new ImageWriter(csf1));
			// app.addTimedStateWriter(new ImageWriter(IL4));
			app.addTimedStateWriter(new SoluteConcentrationWriter());
			app.addTimedStateWriter(new SolidsConcentrationWriter());
			app.addTimedStateWriter(new ParticlePositionWriter());
			// app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getErodedBiomassSeries());
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getSloughedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
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
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		System.out.println("Simulation finished.");
	}
}