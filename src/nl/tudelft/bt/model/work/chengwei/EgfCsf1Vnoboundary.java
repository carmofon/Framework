package nl.tudelft.bt.model.work.chengwei;

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
import nl.tudelft.bt.model.multigrid.boundary_conditions.*;
import nl.tudelft.bt.model.multigrid.boundary_layers.NoBoundaryLayer;
import nl.tudelft.bt.model.particlebased.QSBiomassSpecies;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.particlebased.tumor.TumorModelHandler;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.util.ExtraMath;

/**
 * fourth version of EGF/CSF-1 paracrine loop
 * 
 * @author Chengwei Peng
 */
public class EgfCsf1Vnoboundary extends TumorModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program
	
	// output directory name (
	

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
	protected static String outputDirectory = "/Users/xlabguest/results/v4/test125";
	// Solute species
	// Substrate (S) - the only solute species used here
	protected static float glucoseBulkConcentration = 7.2e-15f; // [g/L]

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	private static float glucoseDiffusivity = 1.8e6f; // [um^2/h]
	// private static float glucoseDiffusivity = 4e10f;
	//protected static float p =0.03f;
	protected static float p =0.15f;
	
	protected static float qEgf =  0.02f; // [nmol/gcell/h]
	protected static float deltaqEgf =  60f; // [nmol/gcell/h]
	protected static float n = 1f;
	protected static float csf1T = 2.5e-15f; // [nM]
	// protected static float decayRate = .75f;

	protected static float qCsf1 =  0.02f;
	protected static float deltaqCsf1 =  60f;
	protected static float m = 1f;
	protected static float egfT = 2.5e-15f;

	protected static float alphaegf = 0.3f / qCsf1; // [nmol/gcell/h]
	protected static float betaegf = 30f;

	protected static float alphaCsf1 = 0.3f / qEgf;
	protected static float betaCsf1 = 30f;
	//
	// Particulate species (biomass X)
	protected static float specificMassX = 70e-15f; // [g/L]

	// Yield coefficients
	private static float YXS = 4.5e-15f; // [gX/gS]
	private static float c = 1e-15f;

	// Processes
	// Growth (biomass production)
	// protected static float uMax = 0.1f; // [1/h]
	protected static float uMax = 3f; // [1/h]

	private static float KS = 3.5e-19f; // [g/L]
	private static float ES = 13.02e-15f; // [nmol/L]
	private static float CS = 13.02e-15f; // [nmol/L]

	// Computation parameters
	protected static float systemSize = 1000; // [um]

	protected static float maxCellRadius = 5; // [um]

	protected static float thicknessOfBoundaryLayer = 0; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = maxCellRadius / systemSize;

	
	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	// protected static float relativeBoundaryLayer = thicknessOfBoundaryLayer
	// systemSize;

	protected static float relativeBoundaryLayer = 0 / systemSize;

	// other model parameters
	//protected static int gridSide = 65; // multigrid grid side
	protected static int gridSide =65; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	// protected static float rdetach = 2e-7f; // NO DETACHMENT PRESENT IN THIS
	// CASE
	protected static float rdetach = 0f; // NO DETACHMENT PRESENT IN THIS CASE

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 10;

	protected static SoluteSpecies glucose;
	protected static SoluteSpecies egf;
	protected static SoluteSpecies csf1;

	protected static BiomassSpecies macrophageBiomass;
	protected static BiomassSpecies cancer_cheaterBiomass;
	protected static BiomassSpecies cancer_producerBiomass;
	protected static BiomassSpecies cancer_cheaterBiomass_s;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {

		glucose = new SoluteSpecies("glucose", glucoseDiffusivity);
		// set up the simplest type of bulk concentration: constant
		glucose.setBulkConcentration(new ConstantBulkConcentration(
				glucoseBulkConcentration));

		egf = new SoluteSpecies("egf", 1.8e5f);
		// set up the simplest type of bulk concentration: constant
		egf.setBulkConcentration(new ConstantBulkConcentration(0));

		csf1 = new SoluteSpecies("csf1", 1.8e5f);
		// set up the simplest type of bulk concentration: constant
		csf1.setBulkConcentration(new ConstantBulkConcentration(0));

		// cancer producer active mass
		ParticulateSpecies cancer_producer = new ParticulateSpecies(
				"cancer_producer", specificMassX, Color.red);
		ParticulateSpecies[] sp_p = { cancer_producer };
		float[] fractionalVolumeCompositionProducers = { 1.0f };

		cancer_producerBiomass = new BiomassSpecies("cancer producer Biomass",
				sp_p, fractionalVolumeCompositionProducers);
		cancer_producerBiomass.setInducibleColor(egf, egfT, Color.red,
				new Color(0.3f, 0.0f, 0.0f));
		cancer_producerBiomass.setActiveMass(cancer_producer);

		ParticulateSpecies cancer_cheater = new ParticulateSpecies(
				"cancer_cheater", specificMassX, Color.blue);
		ParticulateSpecies[] sp_c = { cancer_cheater };
		float[] fractionalVolumeCompositionCheaters = { 1.0f };

		// cancer cheater Biomass
		cancer_cheaterBiomass = new BiomassSpecies("cancer cheater Biomass",
				sp_c, fractionalVolumeCompositionCheaters);
		cancer_cheaterBiomass.setInducibleColor(egf, egfT, Color.blue,
				new Color(0.0f, 0.0f, 0.3f));
		cancer_cheaterBiomass.setActiveMass(cancer_cheater);
		
		ParticulateSpecies cancer_cheater_s = new ParticulateSpecies(
				"cancer_cheater_s", specificMassX, Color.pink);
		ParticulateSpecies[] sp_c_s = { cancer_cheater_s };
		float[] fractionalVolumeCompositionCheaters_s = { 1.0f };

		// cancer cheater Biomass
		cancer_cheaterBiomass_s = new BiomassSpecies("cancer cheater_s Biomass",
				sp_c_s, fractionalVolumeCompositionCheaters_s);
		cancer_cheaterBiomass_s.setInducibleColor(egf, egfT, Color.magenta,
				Color.pink);
		cancer_cheaterBiomass_s.setActiveMass(cancer_cheater_s);

		ParticulateSpecies macrophage = new ParticulateSpecies("macrophage",
				specificMassX, Color.green);

		ParticulateSpecies[] sp_m = { macrophage };
		float[] fractionalVolumeCompositionY = { 1.0f };

		// macrophage biomass
		macrophageBiomass = new BiomassSpecies("macrophageBiomass", sp_m,
				fractionalVolumeCompositionY);
		macrophageBiomass.setInducibleColor(csf1, csf1T, Color.green,
				new Color(0.0f, 0.3f, 0.0f));
		macrophageBiomass.setActiveMass(macrophage);

		Reaction influxGlucose = new Flux("Inflx of glucose", glucose,
				glucoseBulkConcentration, 3);
		influxGlucose.addFactor(new Inhibition(cancer_producer, 10e-20f));
		influxGlucose.addFactor(new Inhibition(cancer_cheater, 10e-20f));
		influxGlucose.addFactor(new Inhibition(macrophage, 10e-20f));

		Reaction decayCsf1f = new Flux("decay of CSF", csf1, 0);
		Reaction decayEgf = new Flux("decay of EGF", egf, 0);

		ProcessFactor mS = new Saturation(glucose, KS);
		ProcessFactor eS = new Saturation(egf, ES);
		ProcessFactor cS = new Saturation(csf1, CS);

		Reaction growth_p = new Reaction("growth_p", cancer_producer, uMax, 2);
		growth_p.addFactor(new Hill(egf, egfT,n));
		growth_p.addFactor(mS);

		Reaction growth_c = new Reaction("growth_c", cancer_cheater, uMax, 2);
		growth_c.addFactor(new Hill(egf, egfT,n));
		growth_c.addFactor(mS);
		
		Reaction growth_c_s = new Reaction("growth_c", cancer_cheater_s, uMax, 2);
		growth_c_s.addFactor(new Hill(egf, egfT,n));
		growth_c_s.addFactor(mS);

		Reaction egf_comsumption_p = new Reaction("comsumption_egf_p",
				cancer_producer, betaegf, 1);
		egf_comsumption_p.addFactor(eS);

		Reaction egf_comsumption_c = new Reaction("comsumption_egf_c",
				cancer_cheater, betaegf, 1);
		egf_comsumption_c.addFactor(eS);
		
		Reaction egf_comsumption_c_s = new Reaction("comsumption_egf_c",
				cancer_cheater_s, betaegf, 1);
		egf_comsumption_c_s.addFactor(eS);

		Reaction csf_comsumption = new Reaction("comsumption_csf_p",
				macrophage, betaCsf1, 1);
		csf_comsumption.addFactor(cS);


		Reaction basalEgf = new Reaction("basal EGF secretion", macrophage,
				qEgf, 0);
		// basalEgf.addFactor(mS);
		Reaction upregEgf = new Reaction("upregulated EGF secretion",
				macrophage, deltaqEgf, 1);
		// upregEgf.addFactor(mS);
		upregEgf.addFactor(new Hill(csf1, csf1T,m));

		Reaction basalCsf1_p = new Reaction("basal CSF-1 secretion producer",
				cancer_producer, qCsf1, 0);
		// basalCsf1_p.addFactor(mS);
		Reaction upregCsf1_p = new Reaction("upregulated CSF secretion",
				cancer_producer, deltaqCsf1, 1);
		upregCsf1_p.addFactor(new Hill(csf1, csf1T,m));
		// upregCsf1_p.addFactor(mS);

		Reaction basalCsf1_c = new Reaction("basal CSF-1 secretion cheater",
				cancer_cheater, qCsf1, 0);
		// basalCsf1_c.addFactor(mS);
		Reaction basalCsf1_c_s = new Reaction("basal CSF-1 secretion cheater",
				cancer_cheater_s, qCsf1, 0);

		// 6. Assign reaction to the species through ReactionStoichiometries
		// cancer producers
		NetReaction rscancer_producer = new NetReaction(2);
		rscancer_producer.addReaction(growth_p, 1);
		//rscancer_producer.addReaction(basalCsf1_p, -c);
		rscancer_producer.addReaction(upregCsf1_p, -c);
		cancer_producer.setProcesses(rscancer_producer);

		NetReaction rscancer_cheater = new NetReaction(1);
		rscancer_cheater.addReaction(growth_c, 1);
		//rscancer_cheater.addReaction(basalCsf1_c, -c);
		cancer_cheater.setProcesses(rscancer_cheater);
		
		NetReaction rscancer_cheater_s = new NetReaction(1);
		rscancer_cheater_s.addReaction(growth_c_s, 1);
		//rscancer_cheater.addReaction(basalCsf1_c, -c);
		cancer_cheater_s.setProcesses(rscancer_cheater_s);

		NetReaction netReactionMacrophage = new NetReaction(0);
		macrophage.setProcesses(netReactionMacrophage);

		NetReaction netReactionCsf1 = new NetReaction(6);
		netReactionCsf1.addReaction(basalCsf1_p, 1);
		netReactionCsf1.addReaction(basalCsf1_c, 1);
		netReactionCsf1.addReaction(basalCsf1_c_s, 1);
		netReactionCsf1.addReaction(upregCsf1_p, 1);
		netReactionCsf1.addReaction(csf_comsumption, -1);
		//netReactionCsf1.addReaction(egf_comsumption_c, -1);
		// netReactionCsf1.addReaction(basalEgf, -alphaCsf1);
		// netReactionCsf1.addReaction(upregEgf, -betaCsf1);
		netReactionCsf1.addReaction(decayCsf1f, 0.0036f);
		csf1.setProcesses(netReactionCsf1);

		NetReaction netReactionEgf = new NetReaction(6);
		netReactionEgf.addReaction(basalEgf, 1);
		netReactionEgf.addReaction(upregEgf, 1);
		netReactionEgf.addReaction(egf_comsumption_p, -1);
		netReactionEgf.addReaction(egf_comsumption_c, -1);
		netReactionEgf.addReaction(egf_comsumption_c_s, -1);
		// netReactionEgf.addReaction(basalCsf1_p, -alphaegf);
		// netReactionEgf.addReaction(basalCsf1_c,-alphaegf);
		// netReactionEgf.addReaction(upregCsf1_p, -betaegf);
		netReactionEgf.addReaction(decayEgf, 0.0036f);
		egf.setProcesses(netReactionEgf);

		NetReaction netReactionGlucose = new NetReaction(4);
		netReactionGlucose.addReaction(growth_p, -(1 / YXS));
		netReactionGlucose.addReaction(growth_c, -(1 / YXS));
		netReactionGlucose.addReaction(growth_c_s, -(1 / YXS));
		netReactionGlucose.addReaction(influxGlucose, 1e-15f);
		glucose.setProcesses(netReactionGlucose);

		// This defines that substrate consumption rate is -(1 / YXS)*rX

		// add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(cancer_producerBiomass);
		addBiomassSpecies(cancer_cheaterBiomass);
		addBiomassSpecies(cancer_cheaterBiomass_s);
		addBiomassSpecies(macrophageBiomass);

		addSoluteSpecies(glucose);
		addSoluteSpecies(csf1);
		addSoluteSpecies(egf);

	}

	/*
	 * (non-Javadoc)
	 */

	public void initializeDiffusionReactionSystem() throws ModelException {
		defineSpeciesAndReactions();
		super.initializeDiffusionReactionSystem();
	}

	protected void createBoundaryLayer(float h)
			throws MultigridSystemNotSetException {
		_boundaryLayer = new NoBoundaryLayer();
		// create the boundary conditions
		MultigridVariable
		.setBoundaryConditions(new TissueBoundaryConditions(p));
		//MultigridVariable
		//.setBoundaryConditions(new SurroundedByBulkBoundaryConditions());
	}

	/*
	 * (non-Javadoc)
	 */
	protected void inoculate() {
		int[] nMacrophages = {75};
		BiomassSpecies[] macrophages = { macrophageBiomass };
		//placeCellsRandomly(nMacrophages, macrophages);
		//placeCellsRandomlybetweenRadius(nMacrophages, macrophages, 50f, 120f);
//		placeCellsAtRadiusandAngle(nMacrophages, macrophages, 80f, 0.7f);
//		placeCellsAtRadiusandAngle(nMacrophages, macrophages, 80f, 0.3f);
//		placeCellsAtRadiusandAngle(nMacrophages, macrophages, 80f, 0.1f);
//		placeCellsAtRadiusandAngle(nMacrophages, macrophages, 80f, 0.9f);
//		placeCellsAtRadiusandAngle(nMacrophages, macrophages, 80f, 0.49f);
		
//		placeCellsAtRadiusandAngle(nMacrophages, macrophages, 120f, 0.7f);
		//placeCellsAtRadiusandAngle(nMacrophages, macrophages, 80f, 0.3f);
		//placeCellsAtRadiusandAngle(nMacrophages, macrophages, 80f, 0.1f);
//		placeCellsAtRadiusandAngle(nMacrophages, macrophages, 120f, 0.9f);
//		placeCellsAtRadiusandAngle(nMacrophages, macrophages, 120f, 0.49f);
		//placeCellsAtRadius(nMacrophages, macrophages, 200f);
		//placeCellsAtRadius(nMacrophages, macrophages, 50f);
		//placeCellsAtRadius(nMacrophages,macrophages, 100f);
		int[] nCancers_m = { 10,10,10};
		BiomassSpecies[] cancers = { cancer_producerBiomass,
				cancer_cheaterBiomass,macrophageBiomass};
		placeCellsAtCenter(nCancers_m, cancers);
		//switchOneRandomnlyPickedCell(cancer_cheaterBiomass, cancer_cheaterBiomass_s);
		//switchOneRandomnlyPickedCell(cancer_cheaterBiomass, cancer_cheaterBiomass_s);
		// int[] nCells = {20,20,20};
		// inoculateRandomly(nCells);
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
		int seed =33;
		int graphicsOn = 1;
		if (args.length == 3) {
			outputDirectory = args[0];
			graphicsOn = Integer.parseInt(args[1]);
			seed = Integer.parseInt(args[2]);
			//c = Float.parseFloat(args[3]);
			//YXS = Float.parseFloat(args[4]);
			System.out.println("Starting simulation with seed " + seed);
			System.out.println("cost:" + c);
			System.out.println("betaegf:" + betaegf);
			System.out.println("qegf:"+ qEgf);
			System.out.println("deltaqEgf:" + deltaqEgf);
			System.out.println("qCsf1:" + qCsf1);
			System.out.println("p"+ p);
			System.out.println("glucosebulk"+ glucoseBulkConcentration);
			System.out.println("glucose consumption"+ YXS);
			
		} else if (args.length > 0) {
			System.out.println("User must supply 3 input arguments");
			System.out.println("1: the output directory");
			System.out.println("2: switch for graphics (1 for on, 0 for off)");
			System.out.println("3: the seed for the random number generator");
			System.out.println("4: betacsf1 values");
			System.exit(-1);
		}
		//TODO change steps for cluster simulations
		MultigridVariable.setSteps(50, 500);
		Model.model().setSeed(seed);
		// create a handle for the application, which will be decorated
		ApplicationComponent app = new EgfCsf1Vnoboundary();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		if (graphicsOn == 1) {
			// The following code will be omitted if no vizuals are desired
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// the biomass thickness visualizer
			app = new SeriesVizualizer(app, biomass);
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
			app.addTimedStateWriter(new ImageWriter(glucose));
			app.addTimedStateWriter(new ImageWriter(egf));
			app.addTimedStateWriter(new ImageWriter(csf1));
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
			//Model.model().setMaxRunLength(600f);
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
