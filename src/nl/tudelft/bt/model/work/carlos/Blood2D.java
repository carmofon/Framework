// CHANGE NAME OF PACKAGE TO MATCH THE LOCAITON IN THE FRAMEWORK
package nl.tudelft.bt.model.work.carlos;

import java.awt.Color;
import java.io.IOException;
import java.util.Iterator;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Radius2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.multigrid.boundary_conditions.GranuleBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.NoBoundaryLayer;
import nl.tudelft.bt.model.particlebased.granule.GranuleBiomassSpecies;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.reaction.*;

/**
 * Tests whether cells can go anywhere in the computational space
 * 
 * @author Joao Xavier (xavierj@mskcc.org)
 */
public class Blood2D extends GranuleModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory = "/Users/CarmonaFontaine/Documents/Xlab/Models/Blood2D";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.

	// basal production
	private static float alpha = 0.1f;
	private static float alpha1 = alpha;
	private static float alpha2 = alpha;
	private static float gamma = 1.9f;
	private static float gamma1 = gamma;
	private static float gamma2 = gamma;
	private static float switchValue = 4e-4f;

	// Nutrient flux parameters
	protected static float externalMassTransferCoefficient = 1e4f; // [h-1]
	// protected static float externalMassTransferCoefficient = 1e-4f; // [h-1]

	// Reaction parameters
	// Ammonium oxydizers (XNH)
	protected static float uMax = 1f; // [h-1]

	private static float KS = 0.3e-3f; // [gO2/L] Merle

	private static float YS = 1.5f; // [gCOD-PHB/gCOD-S]

	// private static float fMaintenance = 0.5f; // fraction of uMax that goes
	// to
	private static float fMaintenance = 1e-3f; // fraction of uMax that goes to

	// maintenance

	// geometry Must be 2 for 2D colony growth
	protected static int geometry = 2;

	// The SBR sycle
	protected static float nutrientConcentration = 3e-3f; // [gO2/L] DO 20

	private static float nutrientDiffusivity = 4e5f; // [um2/h]
	
	
	private static float egfDiffusivity = nutrientDiffusivity; // [um2/h]
	private static float csf1Diffusivity = nutrientDiffusivity; // [um2/h]

	//
	// Particulate species (biomass H)
	protected static float specificMassBiomass = 150f; // [gCOD-H/L]

	// Computation parameters
	// Size of computational volume (size of size of square)
	protected static float systemSize = 500; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 2f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 10 um
	protected static float relativeBoundaryLayer = 0 / systemSize;

	// This turns-off cell shedding from colony
	protected static float maximumColonyRadius = Float.POSITIVE_INFINITY; // [um]

	protected static float initialColonyRadius = 200f; // [um]

	// other model parameters
	protected static int gridSide = 33; // multigrid grid side

	// Don't change this
	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// Don't change this, detachment rate, leave at zero to form round granules
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 100;

	// iteration finish time
	protected static float simulationFinishTime = 5f; // [h]

	// outpute (write results to file) every:
	protected static float outputEvery = 0.2f; // [h]
	

	// /END OF PARAMETERS

	private static SoluteSpecies egf;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	protected void defineSpeciesAndReactions() throws ModelException {
		// nutrient (S)
		SoluteSpecies nutrient = new SoluteSpecies("nutrient",
				nutrientDiffusivity);
		nutrient.setBulkConcentration(new ConstantBulkConcentration(
				nutrientConcentration));
		// EGF - secreted by macrophages
		egf = new SoluteSpecies("egf", egfDiffusivity);
		egf.setBulkConcentration(new ConstantBulkConcentration(0));
		// CSF-1 - secreted by cancer cells
		SoluteSpecies csf1 = new SoluteSpecies("csf1", csf1Diffusivity);
		csf1.setBulkConcentration(new ConstantBulkConcentration(0));
		// 2. Create the particulate species (solids)
		// The active mass of green strain
		ParticulateSpecies activeM = new ParticulateSpecies("activeM",
				specificMassBiomass, Color.green);
		// The active mass of red strain
		ParticulateSpecies activeC = new ParticulateSpecies("activeC",
				specificMassBiomass, Color.red);
		// array of fixed species that constitute speciesH (in this case,
		// speciesH is entirely constituted by active mass)
		// green strain
		ParticulateSpecies[] spM = { activeM };
		float[] fractionalVolumeCompositionM = { 1.0f };
		// 3. Create the biomass species
		BiomassSpecies speciesMacrophage = new ImotileBiomassSpecies(
				"speciesMacrophage", spM, fractionalVolumeCompositionM);
		speciesMacrophage.setActiveMass(activeM);
		speciesMacrophage.setInducibleColor(csf1, switchValue, Color.green,
				new Color(0, 0.5f, 0));
		// red strain
		ParticulateSpecies[] spC = { activeC };
		float[] fractionalVolumeCompositionRed = { 1.0f };
		// 3. Create the biomass species
		BiomassSpecies speciesCancer = new GranuleBiomassSpecies("speciesCancer", spC,
				fractionalVolumeCompositionRed);
		speciesCancer.setActiveMass(activeC);
		speciesCancer.setInducibleColor(egf, switchValue, Color.red,
				new Color(0.5f, 0, 0));
		// 4. Create the Reaction factors
		ProcessFactor mS = new Saturation(nutrient, KS);
		// 5. Create the reactions
		// growth Macrophage
		Reaction aerobicGrowthGreen = new Reaction("growthGreen", activeM,
				uMax, 1);
		aerobicGrowthGreen.addFactor(mS);
		// growth Cancer
		Reaction aerobicGrowthRed = new Reaction("growthRed", activeC, uMax, 1);
		aerobicGrowthRed.addFactor(mS);
		// EGF production
		Reaction egfBasalProduction = new Reaction("egfBasal", activeM, alpha1,
				0);
		Reaction egfUpregulatedProduction = new Reaction("egfUp", activeM,
				gamma1, 1);
		egfUpregulatedProduction.addFactor(new Step(csf1, switchValue));
		// CSF1 production
		Reaction csf1BasalProduction = new Reaction("csf1Basal", activeC,
				alpha2, 0);
		Reaction csf1UpregulatedProduction = new Reaction("egfUp", activeC,
				gamma2, 1);
		csf1UpregulatedProduction.addFactor(new Step(egf, switchValue));
		// fluxes
		Reaction fluxNutrient = new Flux("flux", nutrient,
				nutrientConcentration);
		Reaction fluxEgf = new Flux("fluxEgf", egf, 0);
		Reaction fluxCsf1 = new Flux("fluxCsf1", csf1, 0);

		// 6. Assign reaction to the species through ReactionStoichiometries
		// active Green
		NetReaction rsActiveGreen = new NetReaction(1);
		rsActiveGreen.addReaction(aerobicGrowthGreen, 0f);
		activeM.setProcesses(rsActiveGreen);
		// active Red
		NetReaction rsActiveRed = new NetReaction(1);
		rsActiveRed.addReaction(aerobicGrowthRed, 1);
		activeC.setProcesses(rsActiveRed);
		// assign reaction stoichiometry to the solutes
		// substrate (S)
		NetReaction rsSubstrate = new NetReaction(3);
		rsSubstrate.addReaction(aerobicGrowthGreen, -1 / YS);
		rsSubstrate.addReaction(aerobicGrowthRed, -1 / YS);
		rsSubstrate.addReaction(fluxNutrient, externalMassTransferCoefficient);
		nutrient.setProcesses(rsSubstrate);
		// egf
		NetReaction rsEgf = new NetReaction(3);
		rsEgf.addReaction(egfBasalProduction, 1);
		rsEgf.addReaction(egfUpregulatedProduction, 1);
		rsEgf.addReaction(fluxEgf, externalMassTransferCoefficient);
		egf.setProcesses(rsEgf);
		// csf1
		NetReaction rsCsf1 = new NetReaction(3);
		rsCsf1.addReaction(csf1BasalProduction, 1);
		rsCsf1.addReaction(csf1UpregulatedProduction, 1);
		rsCsf1.addReaction(fluxCsf1, externalMassTransferCoefficient);
		csf1.setProcesses(rsCsf1);
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesMacrophage);
		addBiomassSpecies(speciesCancer);
		addSoluteSpecies(nutrient);
		addSoluteSpecies(egf);
		addSoluteSpecies(csf1);
	}

	
	
	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	
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
	protected void inoculate() {
		int[] nCells = { initialParticleNumber, initialParticleNumber };
		inoculateRandomlyEverywhere(nCells);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		DetachmentSpeedFunction df = new Radius2MassDetachment(kdetach);
		setDetachmentHandler(df);
	}

	/**
	 * @param app
	 */
	protected static void setSystemParametersAndInitializeSystemSpace(
			ApplicationComponent app) {
		// create the space
		app.setSystemSpaceParameters(geometry, systemSize,
				relativeMaximumRadius, relativeMinimumRadius,
				relativeBoundaryLayer, gridSide, kShov);
		// initialize
		app.initializeSystemSpace();
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line
		/*if (args.length < 3) {
			throw new RuntimeException("input arguments missing:\n"
					+ "1: output directory (CAUTION!!!"
					+ " directory will be erased\n"
					+ "2: seed for random number generator\n"
					+ "3: flag for running with graphics (1 on, 0 off)\n"
					+ "4: externalMassTransferCoefficient");
		}
		// set external mass transfer, the 4th, optional, input
		if (args.length > 3) {
			externalMassTransferCoefficient = Float.parseFloat(args[3]);
		}*/

		// parse inputs
		//outputDirectory = args[0] + externalMassTransferCoefficient;
		int seed = 2;
		Model.model().setSeed(seed);
		boolean runWithGraphics = true;
		// set numerics for multigrid
		MultigridVariable.setSteps(2, 20);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new Blood2D();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the biovolume series
		VariableSeries biovolume = new BiovolumeSeries();
		VariableSeries[] runLengthSeries = { new RunLengthXSeries(),
				new RunLengthYSeries(), new RunLengthZSeries() };
		// The following code will be omitted if no vizuals are desired
		if (runWithGraphics) {
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// bulk concentrations
			// app = new BulkConcentrationVizualizer(app);
			// finally, the controller must be the last decorator to add
			// app = new VizualModelControler(app);
		}
		try {
			// create the space
			setSystemParametersAndInitializeSystemSpace(app);
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			// initialize
			app.intializeStateWriters(outputDirectory);
			// Pov witer is added twice
			app
					.addStateWriter(new TimedStateWriterDecorator(
							new PovRayWriter()));
            app.addStateWriter(new ImageWriter());
			app.addStateWriter(new TimedStateWriterDecorator(
					new SoluteConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new SolidsConcentrationWriter()));
			app.addStateWriter(new TimedStateWriterDecorator(
					new ParticlePositionWriter()));
			// app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(biovolume);
			spw.addSeries(runLengthSeries[0]);
			spw.addSeries(runLengthSeries[1]);
			spw.addSeries(runLengthSeries[2]);
			spw.addSeries(Model.model().detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			// add the time constraints writer
			app.addStateWriter(new TimeConstraintsWriter());
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
			e.printStackTrace();
			System.exit(-1);
		}
		try {
			// start iterating cycle
			// Model.model().setCompulsoryTimeStep(outputEvery);
			Model.model().setFinishIterationTime(simulationFinishTime);
			// start the iteration
			app.writeState(); // write iteration 0
			app.startIterating();
		} catch (Exception e1) {
			try {
				app.forceWriteState();
			} catch (IOException e2) {
				System.err.println("Error serializing state:");
				System.err.println("");
				e2.printStackTrace();
			}
			System.err.println("");
			System.err.println("Program failed due to :");
			System.err.println("");
			e1.printStackTrace();
			System.out.println(e1);
		}
		System.out.println("Simulation finished.");
		System.exit(0);
	}
}