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
public class Angiogenesis_May2016 extends TumorModelHandler {
//	public class Angiogenesis_May2016 extends ModelHandler {

	// public class Example1KWithTumorHandler extends TumorModelHandler {
	// All the model parameters are defined here as static attributes
	// at the beginning of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/CarmonaFontaine/Documents/Xlab/Models/160503_Angiogenesis/Test";

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

	protected static float systemSize = 1000; // [um]

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

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 0f; // 

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 100;

	///////////// PARTICULATES //

	protected static int numberofCellGroups = 2;

	protected static float specificMassC = 100f; // [g/L]
	protected static float specificMassB = 100f; // [g/L]
	
	
	///////////// SOLUTES //

	protected static SoluteSpecies oxygen;
	protected static SoluteSpecies lactate;
	protected static SoluteSpecies vegf;

	
	
	///////////// PARAMETERS //

	protected static float oxygenBulkConcentration = 0f; // [g/L]Joao's paper 
															// Should be 16-22 ml O2/100ml arterial blood
															// and 13-16 ml in venous blood ()
	protected static float lactateBulkConcentration = 0f; // [g/L]
	protected static float vegfBulkConcentration = 0f; 
	

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	// Glucose diffusivity (Casciari et al 1988) 1.1e-6 cm2/s = 3.96e5um2/h
	// mouse; 5.5e-7-2.3e-7 cm2/s = 1.98e4um2/h to 8.2e3 um2/h
	//private static float oxygenDiffusivity = 7.2e7f; // 10x Diffusivity
	private static float oxygenDiffusivity = 7.2e6f;//[um^2/h]Thomlinson: 2e-5
														// cm2/s
	private static float lactateDiffusivity = 1.04e6f; //1.04e6f [um^2/h] 2.9e-6 cm2/s.
														// Measurement of
														// Effective
														// Diffusivities of
														// Lactose and Lactic
														// Acid in 3% Agarose
														// Gel Membrane
														// Amarjeet S. Bassi,
														// Sohrab Rohani and D.
														// G. Macdonaid
	private static float vegfDiffusivity = 3.7e5f; // 3.7e5f 113um^2/s = 374400m^2/h
	//from: Computational Model of Vascular Endothelial Growth Factor Spatial Distribution in Muscle and Pro-Angiogenic Cell Therapy

	// protected static SoluteSpecies Hydrogen;

	/**
	 * Variables
	 */

	// /////////// Biochemical parameters

	// Yield coefficients

	private static float alpha = 0.9f; // [gX/gS] Oxygen Yield
	private static float beta = 0.1f; // [gX/gS] Lactate Production Yield 
	private static float gamma = 0.1f; // [gX/gS] Acidosis Yield
	private static float delta = 0.1f; // [gX/gS] Angiogenesis cost

	
	// Saturation constants

	private static float ko2 = 30f / 2f; // [g/L]
	private static float ko2v = ko2;
	private static float kLac = 1e-2f; // [g/L] saturation constant lactate killing effect
	private static float kLacv = kLac *10f; // [g/L] saturation constant lactate effect on vegf secretion



	// Rates
	protected static float uM1 = 0.1f; // [1/h] max aerobic growth rate
	protected static float uM2 = uM1 * 0.001f; // [1/h] max anaerobic growth rate
	protected static float uM3 = uM1 * 100f; // [1/h] max killing rate by acidosis
	protected static float uM4 = uM1 * 10f; // [1/h] max oxygen perfusion rate
	protected static float uM5 = uM1 * 10e7f; // [1/h] max lactate perfusion rate
	protected static float uM6 = uM1 * 10f; // [1/h] max vegf secretion rate
	private static float lambda = 5f; // Angiogenic factor decay
	


	
	/**
	 * Definitions and Processes
	 */

	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// Oxygen
		oxygen = new SoluteSpecies("Oxygen", oxygenDiffusivity);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(oxygenBulkConcentration));
			
		// Lactate
		lactate = new SoluteSpecies("Lactate", lactateDiffusivity);
		lactate.setBulkConcentration(new ConstantBulkConcentration(lactateBulkConcentration));
	
		// VEGF
		vegf = new SoluteSpecies("VEGF", vegfDiffusivity);
		vegf.setBulkConcentration(new ConstantBulkConcentration(vegfBulkConcentration));

		// 2. Create the particulate species (solids)

		// X active mass
		ParticulateSpecies activeC = new ParticulateSpecies("activeC",
				specificMassC, Color.blue);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spC = { activeC };
		float[] fractionalVolumeCompositionC = { 1.0f };

		ParticulateSpecies activeB = new ParticulateSpecies("activeB",
				specificMassB, Color.red);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spB = { activeB };
		float[] fractionalVolumeCompositionB = { 1.0f };
		
	

		// 3. Create the biomass species
		BiomassSpecies speciesC = new BiomassSpecies("speciesC", spC,
				fractionalVolumeCompositionC);
		speciesC.setActiveMass(activeC);
		//speciesC.getColorFromGrowth();

		BiomassSpecies speciesB = new BiomassSpecies("speciesB", spB,
				fractionalVolumeCompositionB);
		speciesB.setActiveMass(activeB);
		//speciesM.setInducibleColor(oxygen, switchValue, Color.green,
			//	new Color(1f, 0, 1f));//HOW MAKE DEPENDENT OF LACTATE + OXYGEN
		// speciesTAM.getColorFromGrowth();


		// 4. Create the Reaction factors, Monod and inhibition coefficients

		ProcessFactor monodAeroGrowth = new Saturation(oxygen, ko2);
		ProcessFactor monodAnaGrowth = new Inhibition(oxygen, 0.1f *ko2);
		ProcessFactor monodAcidosis = new Saturation(lactate, kLac); //
		ProcessFactor linearVegfPerfusion = new Linear(vegf);
		ProcessFactor linearLactatePerfusion = new Linear(lactate);
		ProcessFactor monodVegfSecretionLac = new Saturation(lactate,kLacv);
		ProcessFactor inhibitionVegfSecretionOxy = new Inhibition(oxygen,ko2v);
//		ProcessFactor linearVegfCost = new Saturation(vegf, kLac); 

		ProcessFactor monodBloodGrowth = new Saturation(oxygen, ko2); 

		// 5. Create the reactions
		Reaction aeroGrowthC = new Reaction("aeroGrowthC", activeC, uM1, 1);
		aeroGrowthC.addFactor(monodAeroGrowth);
		
		Reaction anaGrowthC = new Reaction("anaGrowthC", activeC, uM2, 1);
		anaGrowthC.addFactor(monodAnaGrowth);
		
		Reaction LactateKilling = new Reaction("LactateKilling", activeC, uM3, 1);
		LactateKilling.addFactor(monodAcidosis);

		Reaction bloodOxygenPerfusion = new Reaction("bloodOxygenPerfusion", activeB, uM4, 1);
		bloodOxygenPerfusion.addFactor(linearVegfPerfusion);
		
		Reaction bloodLactatePerfusion = new Reaction("bloodLactatePerfusion", activeB, uM5, 2);
		bloodLactatePerfusion.addFactor(linearLactatePerfusion);
		bloodLactatePerfusion.addFactor(linearVegfPerfusion);
		
		Reaction vegfSecretion = new Reaction("vegfSecretion", activeC, uM6,2);
		vegfSecretion.addFactor(monodVegfSecretionLac);
		vegfSecretion.addFactor(inhibitionVegfSecretionOxy);
		
		Reaction bloodGrowthB = new Reaction("bloodGrowthB", activeB, 0f, 1);
		bloodGrowthB.addFactor(monodBloodGrowth);
		
		Reaction vegfDecay = new Flux("vegf decay", vegf, 0f);
		
		
		// 6. Assign reaction to the species through ReactionStoichiometries
		NetReaction rsXactiveC = new NetReaction(4);
		rsXactiveC.addReaction(aeroGrowthC, 1);
		rsXactiveC.addReaction(anaGrowthC, 1);
		rsXactiveC.addReaction(LactateKilling, -(gamma));// *** ASKK is it ok to put stoichiometry here??
		rsXactiveC.addReaction(vegfSecretion, -(delta));
		activeC.setProcesses(rsXactiveC);

		NetReaction rsXactiveB = new NetReaction(1);
		rsXactiveB.addReaction(bloodGrowthB, 1);
		activeB.setProcesses(rsXactiveB);

		NetReaction rsOxygen = new NetReaction(2);
		rsOxygen.addReaction(aeroGrowthC, -(alpha));
		rsOxygen.addReaction(bloodOxygenPerfusion, 1);
		oxygen.setProcesses(rsOxygen);

		NetReaction rsLactate = new NetReaction(2);
		rsLactate.addReaction(anaGrowthC, (beta));	
		rsLactate.addReaction(bloodLactatePerfusion, -1);
		lactate.setProcesses(rsLactate);
		
		NetReaction rsVegf = new NetReaction(2);
		rsVegf.addReaction(anaGrowthC, (delta)); // OR SHOULD BE ONE?		
		rsVegf.addReaction(vegfDecay, lambda);
		vegf.setProcesses(rsVegf);
		
		
		
		addBiomassSpecies(speciesB);
		addBiomassSpecies(speciesC);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(lactate);
		addSoluteSpecies(vegf);

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
		int[] nCells = { initialParticleNumber, initialParticleNumber};
		// inoculateRandomlyInsideRadius(nCells, initialColonyRadius);
		//inoculateRandomlyEverywhere(nCells);
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

//	protected void createBoundaryLayer(float h)
//			throws MultigridSystemNotSetException {
//		_boundaryLayer = new FixedOnTopBoundaryLayer();
//		// create the boundary conditions
//		MultigridVariable
//				.setBoundaryConditions(new BiofilmBoundaryConditions());
//	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new Angiogenesis_May2016();
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
			app.addTimedStateWriter(new ImageWriter(vegf));

			app.addTimedStateWriter(new SoluteConcentrationWriter());
			app.addTimedStateWriter(new SolidsConcentrationWriter());
			app.addTimedStateWriter(new ParticlePositionWriter());

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