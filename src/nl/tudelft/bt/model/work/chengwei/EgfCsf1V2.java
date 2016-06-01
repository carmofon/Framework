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
import nl.tudelft.bt.model.particlebased.QSBiomassSpecies;
import nl.tudelft.bt.model.particlebased.granule.GranuleModelHandler;
import nl.tudelft.bt.model.particlebased.tumor.TumorModelHandler;
import nl.tudelft.bt.model.reaction.*;

/**
 * First version of EGF/CSF-1 paracrine loop
 * 
 * @author Chengwei Peng
 */
public class EgfCsf1V2 extends TumorModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	// output directory name (
	protected static String outputDirectory = "/Users/xlabguest/results/v2_preliminary/test4";

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

	//UNITS: mass:g length:um time:h
	// Solute species
	// Substrate (S) - the only solute species used here
	protected static float glucoseBulkConcentration = 0.72f; // [g/L]

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	private static float glucoseDiffusivity = 1.8e6f; // [um^2/h]

	//

	//
	// Particulate species (biomass X)
	protected static float specificMassX = 70f; // [g/L]

	protected static float qEgf = 1.3886e-4f;   //[nmol/g/h]
	protected static float deltaqEgf = 1.3886f;	//[nmol/g/h]
	protected static float n = 2f;
	protected static float csf1T = 50f;         //[nM]

	protected static float qCsf1 = 1.3886e-4f;
	protected static float deltaqCsf1 = 1.3886f;
	protected static float m = 1.7f;
	protected static float egfT = 50f;

	// Yield coefficients
	private static float YXS = 0.045f; // [gX/gS]

	// Cost of CSF-1
	private static float c = 1e-12f; // [gX/gS]

	
	// Processes
	// Growth (biomass production)
	protected static float uMax = 0.1f; // [1/h]
	// protected static float uMax = 0.0547f; //[1/h]

	private static float KS = 3.5e-4f; // [g/L]

	// Computation parameters
	protected static float systemSize = 1000; // [um]
	
	protected static float maxCellRadius = 5; //[um]

	protected static float thicknessOfBoundaryLayer = 10; //[um]
	

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = maxCellRadius/systemSize;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = thicknessOfBoundaryLayer/systemSize;

	// other model parameters
	protected static int gridSide = 65; // multigrid grid side
	// protected static int gridSide = 33; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	//protected static float rdetach = 2e-7f; // NO DETACHMENT PRESENT IN THIS CASE
	protected static float rdetach = 0f; // NO DETACHMENT PRESENT IN THIS CASE

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 10;

	protected static SoluteSpecies glucose;
	protected static SoluteSpecies egf;
	protected static SoluteSpecies csf1;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		
		// glucose
		glucose = new SoluteSpecies("glucose", glucoseDiffusivity);
		// set up the simplest type of bulk concentration: constant
		glucose.setBulkConcentration(new ConstantBulkConcentration(
				glucoseBulkConcentration));
		// EGF
		egf = new SoluteSpecies("egf", 1.8e5f);
		// set up the simplest type of bulk concentration: constant
		egf.setBulkConcentration(new ConstantBulkConcentration(0.005f));
		// CSF-1
		csf1 = new SoluteSpecies("csf1", 1.8e5f);
		// set up the simplest type of bulk concentration: constant
		csf1.setBulkConcentration(new ConstantBulkConcentration(0));

		// 2. Create the particulate species (solids)
		
		// cancer producer active mass
		ParticulateSpecies cancer_producer = new ParticulateSpecies("cancer_producer",
				specificMassX, Color.red);
		ParticulateSpecies[] sp_p = { cancer_producer };
		float[] fractionalVolumeCompositionProducers = { 1.0f };

		// cancer producer Biomass
		BiomassSpecies cancer_producerBiomass = new QSBiomassSpecies(
				"cancer producer Biomass", sp_p, fractionalVolumeCompositionProducers, egf, egfT,
				new Color(0.5f, 0.0f, 0.0f), Color.red);
		cancer_producerBiomass.setActiveMass(cancer_producer);
		
		//cancer cheater active mass
		ParticulateSpecies cancer_cheater = new ParticulateSpecies("cancer_cheater",
				specificMassX, Color.blue);
		ParticulateSpecies[] sp_c = { cancer_cheater };
		float[] fractionalVolumeCompositionCheaters = { 1.0f };

		// cancer cheater Biomass
		BiomassSpecies cancer_cheaterBiomass = new QSBiomassSpecies(
				"cancer cheater Biomass", sp_c, fractionalVolumeCompositionCheaters, egf, egfT,
				new Color(0.0f, 0.0f, 0.5f), Color.blue);
		cancer_cheaterBiomass.setActiveMass(cancer_cheater);
		

		// Create macrophages
		ParticulateSpecies macrophage = new ParticulateSpecies("macrophage",
				specificMassX, Color.green);
		
		ParticulateSpecies[] sp_m = { macrophage};
		float[] fractionalVolumeCompositionY = { 1.0f};

		// macrophage biomass
		BiomassSpecies macrophageBiomass = new QSBiomassSpecies(
				"macrophageBiomass", sp_m, fractionalVolumeCompositionY, csf1, csf1T,
				new Color(0.0f, 0.5f, 0.0f), Color.green);
		macrophageBiomass.setActiveMass(macrophage);
		

		// 4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mS = new Saturation(glucose, KS);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate

		// 5. Create the reactions
		// growth of cancer producer cells
		Reaction growth_p = new Reaction("growth_p", cancer_producer, uMax, 2);
		growth_p.addFactor(mS);
		growth_p.addFactor(new Hill(egf, egfT, m));
		
		// growth of cancer cheater cells
		Reaction growth_c = new Reaction("growth_c", cancer_cheater, uMax, 2);
		growth_c.addFactor(mS);
		growth_c.addFactor(new Hill(egf, egfT, m));
		
		// secretion of EGF
		Reaction basalEgf = new Reaction("basal EGF secretion", macrophage,
				qEgf, 0);
		Reaction upregEgf = new Reaction("upregulated EGF secretion",
				macrophage, deltaqEgf, 1);
		upregEgf.addFactor(new Hill(csf1, csf1T, n));
		
		
		// secretion of CSF-1 by producers
		Reaction basalCsf1_p = new Reaction("basal CSF-1 secretion producer", cancer_producer,
				qCsf1, 0);
		Reaction upregCsf1_p = new Reaction("upregulated CSF secretion", cancer_producer,
				deltaqCsf1, 1);
		upregCsf1_p.addFactor(new Hill(egf, egfT,m));
		
		//secretion of CSF1 by cheaters (basal only)
		Reaction basalCsf1_c = new Reaction("basal CSF-1 secretion cheater", cancer_producer,
				qCsf1, 0);
		

		// 6. Assign reaction to the species through ReactionStoichiometries
		// cancer producers
		NetReaction rscancer_producer = new NetReaction(3);
		rscancer_producer.addReaction(growth_p, 1);
		rscancer_producer.addReaction(basalCsf1_p, -c);
		rscancer_producer.addReaction(upregCsf1_p, -c);		
		cancer_producer.setProcesses(rscancer_producer);
		
		//cancer cheaters
		NetReaction rscancer_cheater = new NetReaction(2);
		rscancer_cheater.addReaction(growth_c, 1);
		rscancer_cheater.addReaction(basalCsf1_c, -c);	
		cancer_cheater.setProcesses(rscancer_cheater);

		//macrophages
		NetReaction netReactionMacrophage = new NetReaction(0);
		macrophage.setProcesses(netReactionMacrophage);
		
	
		// Netreactions for signals
		NetReaction netReactionEgf = new NetReaction(2);
		netReactionEgf.addReaction(basalEgf, 1);
		netReactionEgf.addReaction(upregEgf, 1);
		egf.setProcesses(netReactionEgf);

		NetReaction netReactionCsf1 = new NetReaction(3);
		netReactionCsf1.addReaction(basalCsf1_p, 1);
		netReactionCsf1.addReaction(upregCsf1_p, 1);
		netReactionCsf1.addReaction(basalCsf1_c, 1);
		csf1.setProcesses(netReactionCsf1);

		
		NetReaction netReactionGlucose = new NetReaction(2);
		netReactionGlucose.addReaction(growth_p, -(1 / YXS));
		netReactionGlucose.addReaction(growth_c, -(1 / YXS));
		glucose.setProcesses(netReactionGlucose);
		// This defines that substrate consumption rate is -(1 / YXS)*rX
		
		// add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(cancer_producerBiomass);
		addBiomassSpecies(cancer_cheaterBiomass);
		addBiomassSpecies(macrophageBiomass);
		addSoluteSpecies(glucose);
		addSoluteSpecies(egf);
		addSoluteSpecies(csf1);
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
		int[] nCells = {20,20,20};
		inoculateRandomly(nCells);
		//inoculateRandomlyInGranuleAtRadius(nCells, 200f);
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
		
		MultigridVariable.setSteps(100, 1000);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new EgfCsf1V2();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		// The following code will be omitted if no vizuals are desired
		// start decorationg the application
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		app = new SeriesVizualizer(app, biomass);
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
			app.addTimedStateWriter(new ImageWriter(csf1));
			app.addTimedStateWriter(new ImageWriter(egf));
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