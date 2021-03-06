package nl.tudelft.bt.model.work.phobia;

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
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.reaction.*;

/**
 * Example 1: Version with a single organism phototrophic)
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class Phobia3 extends ModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/xavierj/results/test/";

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
	// Carbon Dioxide(CO2)
	protected static float CO2BulkConcentration = 1e-5f; // [g/L]
	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	private static float CO2Diffusivity = 5.76e6f; // [um^2/h]
	// Oxygen (O2)
	protected static float O2BulkConcentration = 0.0f; // [g/L]
	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	private static float O2Diffusivity = 5.76e6f; // [um^2/h]
	//

	//
	// Particulate species (biomass X)
	protected static float specificMassBio = 170f; // [g/L]
	protected static float specificMassEPS = 170f/6f; // [g/L]
	// Yield coefficients
	private static float YXS = 0.045f; // [gX/gS]

	// Processes
	// Growth (biomass production)
	protected static float uMaxPhot = 0.1f; // [1/h]
	protected static float uMaxHet = 0.5f; // [1/h]
	// protected static float uMax = 0.0547f; //[1/h]

	private static float KS = 3.5e-4f; // [g/L]
	private static float KSO2 = 6.25e-7f; // [g/L]

	// Computation parameters
	protected static float systemSize = 100; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.005f;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 0.05f;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	protected static float rdetach = 0; // NO DETACHMENT PRESENT IN THIS CASE

	// initial number of particles in the system (inoculum)
	protected static int[] initialParticleNumber = {10,10};

	protected static SoluteSpecies CO2;
	
	protected static SoluteSpecies O2;
	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes

		// CO2
		CO2 = new SoluteSpecies("CO2",
				CO2Diffusivity);
		// set up the simplest type of bulk concentration: constant
		CO2.setBulkConcentration(new ConstantBulkConcentration(
				CO2BulkConcentration));
		// O2
		O2 = new SoluteSpecies("O2",
						O2Diffusivity);
		// set up the simplest type of bulk concentration: constant
		O2.setBulkConcentration(new ConstantBulkConcentration(
						O2BulkConcentration));
		
		// 2. Create the particulate species (solids)
		// Phototrophic biomass
		ParticulateSpecies PhotoMass = new ParticulateSpecies("PhotoMass",
				specificMassBio, Color.green);
		ParticulateSpecies EPS = new ParticulateSpecies("EPS",
				specificMassEPS, Color.yellow);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spPhot = { PhotoMass, EPS };
		float[] fractionalVolumeCompositionP = { 1.0f, 0.0f };

		// 3. Create the biomass species
		// Phototrophs
		BiomassSpecies Photos = new BiomassSpecies("Photos", spPhot,
				fractionalVolumeCompositionP);
		Photos.setActiveMass(PhotoMass);
		Photos.setEpsMass(EPS);
		Photos.getColorFromGrowth();
		
		// Heterotrophs
		ParticulateSpecies HetMass = new ParticulateSpecies("HetMass",
				specificMassBio, Color.red);
		
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spHet = { HetMass };
		float[] fractionalVolumeCompositionH = { 1.0f };

		// 3. Create the biomass species
		BiomassSpecies Hets = new BiomassSpecies("Hets", spHet,
				fractionalVolumeCompositionH);
		Hets.setActiveMass(HetMass);

		

		// 4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mCO2 = new Saturation(CO2, KS);
		ProcessFactor mO2 = new Saturation(O2, KSO2);
		ProcessFactor heightDependent = new HeightDependent(10);
		// Phototrophs
		// Phototrophs
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate

		// 5. Create the reactions
		// growth phototrophs
		Reaction growthPhoto = new Reaction("growthPhoto", PhotoMass, uMaxPhot, 2);
		growthPhoto.addFactor(mCO2);
		growthPhoto.addFactor(heightDependent);
		// growth heterotrophs
		Reaction growthHet = new Reaction("growthHet", HetMass, uMaxHet, 1);
		growthHet.addFactor(mO2);
		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass
		//
		
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactive = new NetReaction(1);
		rsXactive.addReaction(growthPhoto, 1);
		PhotoMass.setProcesses(rsXactive);
		//
		NetReaction rsEPS = new NetReaction(1);
		rsEPS.addReaction(growthPhoto, 0.5f);
		EPS.setProcesses(rsEPS);
		//
		NetReaction rsHets = new NetReaction(1);
		rsHets.addReaction(growthHet, 1f);
		HetMass.setProcesses(rsHets);
		// This defines that biomass growth rate is 1*rX
		//
		// assign reaction stoichiometry to the solutes
		// substrate
		NetReaction rsCO2 = new NetReaction(1);
		rsCO2.addReaction(growthPhoto, -(1 / YXS));
		CO2.setProcesses(rsCO2);
		//
		NetReaction rsO2 = new NetReaction(2);
		rsO2.addReaction(growthPhoto, 1);
		rsO2.addReaction(growthHet, -0.01f);
		O2.setProcesses(rsO2);
		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(Photos);
		addBiomassSpecies(Hets);
		addSoluteSpecies(CO2);
		addSoluteSpecies(O2);
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
		inoculateRandomly(initialParticleNumber);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		boolean graphicsOn = true;
		if (args.length >= 2) {
			initialParticleNumber[0] = Integer.parseInt(args[0]);
			initialParticleNumber[0] = Integer.parseInt(args[0]);
			graphicsOn = (Integer.parseInt(args[1]) == 1) ;
		}
		MultigridVariable.setSteps(10, 500);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new Phobia3();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		if (graphicsOn) {
			// The following code will be omitted if no vizuals are desired
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
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			app.intializeStateWriters(outputDirectory);
			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new ImageWriter(O2));
			app.addTimedStateWriter(new SoluteConcentrationWriter());
			app.addTimedStateWriter(new SolidsConcentrationWriter());
			app.addTimedStateWriter(new ParticlePositionWriter());
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