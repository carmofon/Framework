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
public class Angiogenesis_OrthogonalFlow_June2016 extends GranuleModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name
	protected static String outputDirectory = "/Users/CarmonaFontaine/Documents/Xlab/Models/AngioOrthogonal/Targeted";

	// WARNING: the contents of the outputdirectory will be deleted!!
	// Be sure not to choose a directory were you have important information
	// stored.
	
	// Computation parameters
	// Size of computational volume (size of size of square)
	protected static float systemSize = 1000; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 5f / systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.2f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 10 um
	protected static float relativeBoundaryLayer = 0 / systemSize;

	protected static float maximumColonyRadius = Float.POSITIVE_INFINITY; // [um]
	protected static float initialColonyRadius = 200f; // [um]

	// other model parameters
	protected static int gridSide = 33; // multigrid grid side

	// Don't change this
	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// Don't change this, detachment rate, leave at zero to form round granules
	protected static float kdetach = 0f; // [1e-15 gCOD-H/um^4/h]

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 50;

	// iteration finish time
	protected static float simulationFinishTime = 500f; // [h]

	// outpute (write results to file) every:
	protected static float outputEvery = 0.2f; // [h]
	
	
	
	
	

	// Yield coefficients

	private static float alpha = 1f; // [gX/gS] 1/Oxygen Yield
	private static float delta = 4f; // [gX/gS] Angiogenesis cost
	
	// Saturation constants

	private static float ko2 = 1f; // [g/L]
	
	////************* SWITCH ko2v to try different models
	private static float ko2v = 0.25f*ko2;// [g/L] saturation constant on vegf secretion
	//private static float ko2v = 1e20f;// [g/L] saturation constant on vegf secretion


	// Rates
	protected static float uM1 = 1f; // [1/h] max aerobic growth rate
	protected static float uM4 = uM1 * 1e6f; // [1/h] max oxygen perfusion rate
	protected static float uM6 = uM1 * 0.2f; // [1/h] max vegf secretion rate
	private static float lambda = 0.1f; // Angiogenic factor decay

	// Nutrient flux parameters
	protected static float externalMassTransferCoefficient = 1e2f; // [h-1]
	// protected static float externalMassTransferCoefficient = 1e-4f; // [h-1]


	
	private static float fMaintenance = 0f; // fraction of uMax that goes to maintenance

	// geometry Must be 2 for 2D colony growth
	protected static int geometry = 2;

	// Nutrient properties
	//protected static float nutrientConcentration = 3e-3f; // [gO2/L] DO 20
	protected static float nutrientBulkConcentration = 0f; // [gO2/L] DO 20
	private static float nutrientDiffusivity = 4e5f; // [um2/h]
	
	protected static float vegfBulkConcentration = 0f; 
	private static float vegfDiffusivity = 3.7e4f; // [um2/h]
	
	

	//
	// Particulate species (biomass H)
	protected static float specificMassBiomass = 150f; // [gCOD-H/L]


	

	///////////// PARTICULATES //

	protected static int numberofCellGroups = 2;
	protected static float specificMassC = 100f; // [g/L]
	protected static float specificMassB = 100f; // [g/L]
	
	
	///////////// SOLUTES //

	protected static SoluteSpecies nutrient;
	protected static SoluteSpecies vegf;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	protected void defineSpeciesAndReactions() throws ModelException {
		// nutrient (S)
		nutrient = new SoluteSpecies("nutrient", nutrientDiffusivity);
		nutrient.setBulkConcentration(new ConstantBulkConcentration(nutrientBulkConcentration));
	
		
		// VEGF - or any angiogenic molecule
		vegf = new SoluteSpecies("vegf", vegfDiffusivity);
		vegf.setBulkConcentration(new ConstantBulkConcentration(vegfBulkConcentration));
		
		
		
		// 2. Create the particulate species (solids)
		// The active mass of green strain
		
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
		speciesC.setInducibleColor(vegf, uM6*(0.25f*ko2), Color.yellow,
				new Color(0, 0, 1f));
		//speciesC.getColorFromGrowth();

		BiomassSpecies speciesB = new ImotileBiomassSpecies("speciesB", spB,
				fractionalVolumeCompositionB);
		speciesB.setActiveMass(activeB);
		//speciesM.setInducibleColor(oxygen, switchValue, Color.green,
			//	new Color(1f, 0, 1f));//HOW MAKE DEPENDENT OF LACTATE + OXYGEN
		// speciesTAM.getColorFromGrowth();
		
		
		// 4. Create the Reaction factors, Monod and inhibition coefficients

		ProcessFactor monodAeroGrowth = new Saturation(nutrient, ko2);
		ProcessFactor linearVegfPerfusion = new Linear(vegf);
		ProcessFactor inhibitionVegfSecretionOxy = new Inhibition(nutrient,ko2v);
		ProcessFactor monodBloodGrowth = new Saturation(nutrient, ko2); 
		
	
		
			
		// 5. Create the reactions
		Reaction aeroGrowthC = new Reaction("aeroGrowthC", activeC, uM1, 1);
		aeroGrowthC.addFactor(monodAeroGrowth);
		
		Reaction bloodOxygenPerfusion = new Reaction("bloodOxygenPerfusion", activeB, uM4, 1);
		bloodOxygenPerfusion.addFactor(linearVegfPerfusion);
		
		Reaction vegfSecretion = new Reaction("vegfSecretion", activeC, uM6, 1);
		vegfSecretion.addFactor(inhibitionVegfSecretionOxy);
		
		Reaction bloodGrowthB = new Reaction("bloodGrowthB", activeB, 0, 1);
		bloodGrowthB.addFactor(monodBloodGrowth);

		// fluxes
		//Reaction vegfDecay = new Flux("vegf decay", vegf, lambda);
		Reaction vegfDecay = new Flux("vegf decay", vegf, 0);
		Reaction fluxNutrient = new Flux("flux", nutrient, nutrientBulkConcentration);

		
		
		// 6. Assign reaction to the species through ReactionStoichiometries
		NetReaction rsXactiveC = new NetReaction(2);
		rsXactiveC.addReaction(aeroGrowthC, 1);
		rsXactiveC.addReaction(vegfSecretion, -(delta));
		activeC.setProcesses(rsXactiveC);
		
		NetReaction rsXactiveB = new NetReaction(1);
		rsXactiveB.addReaction(bloodGrowthB, 1);
		activeB.setProcesses(rsXactiveB);

		NetReaction rsNutrient = new NetReaction(3);
		rsNutrient.addReaction(aeroGrowthC, -(alpha));
		rsNutrient.addReaction(bloodOxygenPerfusion, 1);
		rsNutrient.addReaction(fluxNutrient, externalMassTransferCoefficient);
		nutrient.setProcesses(rsNutrient);
		
		NetReaction rsVegf = new NetReaction(2);
		rsVegf.addReaction(vegfSecretion, (1));		
	//	rsVegf.addReaction(vegfDecay, (-1));
		rsVegf.addReaction(vegfDecay, externalMassTransferCoefficient);
		vegf.setProcesses(rsVegf);
		
		addBiomassSpecies(speciesC);
		addBiomassSpecies(speciesB);				
		addSoluteSpecies(nutrient);
		addSoluteSpecies(vegf);

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
		int[] nCells = { initialParticleNumber, 3*initialParticleNumber };
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
		ApplicationComponent app = new Angiogenesis_OrthogonalFlow_June2016();
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
			app.addStateWriter(new TimedStateWriterDecorator(
							new PovRayWriter()));
            app.addStateWriter(new ImageWriter());
			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new ImageWriter());
			app.addTimedStateWriter(new ImageWriter(nutrient));
			app.addTimedStateWriter(new ImageWriter(vegf));

				
            
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