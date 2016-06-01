package nl.tudelft.bt.model.tests.detachmentstudy;
import java.awt.Color;
import java.text.DecimalFormat;
import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.Height2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;
/**
 * Simple 2D monospecies system, can be used to study the effect of diffusion
 * limitation on biofilm structure (change growth rate and density, for
 * instance) With EPS production
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DRGDetachmentEpsDecayApp extends ModelHandler {
	//output directory name
	private final static String outputDirectory = "D:\\results\\drg-detachment-eps-decay\\";
	// geometry (2D or 3D)
	private static int geometry = 2;
	// heterptroph biomass
	private static float uMax = 0.547f; // default val 5.47e-2f
	private static final float YH = -0.045f; // [g O2 / g H]
	private static float densityH = 70.0f; // [g/l = 10^-15g/um^3]
	private static float inactivationConstantH = 1e-3f; // [h^-1]
	// EPS
	private static final float YEPS = 0.5f; // [g EPS / g h]
	private static float densityEPS = 11.67f; // [3.5f g/l = 10^-15g/um^3]
	private static float decayRateEPS = 1.4e-3f; // [1.4e-3 h-1]
	// Oxygen
	private static final float diffusivityO2 = 5.76e6f; //[um2/h]
	private static float constantBulkConcentrationO2 = 4e-3f; //[2.4f g/l]
	private static final float KO2 = 3.5e-4f; //[gO2/l]
	private static final float mO2 = -1.08e-1f; //[KgO2/kgH/h]
	// Operational parameters
	private static float rdetach = 1.0e-3f; // detachment constant[g/l/h]
	private static float relativeBoundaryLayer = 0.1f;
	//simulation parameters
	private static float systemSize = 500; // [um]
	//Dimensionless numbers
	private static float relativeMaximumRadius = 0.008f;
	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;
	private static final float minimumMassRatio = 0.001f;
	private static int gridSide = 33; // multigrid grid side
	private static float kShov = 1.0f; // shoving parameter[dim/less]
	private static int initialCellNumber = 60; // inoculation
	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// create the chemicals
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", diffusivityO2);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				constantBulkConcentrationO2));
		// create the species
		// active mass
		ParticulateSpecies activeMassH = new ParticulateSpecies(
				"heterotrophActiveMass", densityH, Color.gray);
		// EPS
		ParticulateSpecies eps = new ParticulateSpecies("eps", densityEPS, Color.gray);
		// inert material
		ParticulateSpecies inert = new ParticulateSpecies("inert", densityH, Color.gray);
		// array of fixed species that constitute the denitrifier
		ParticulateSpecies[] spa = {activeMassH, eps, inert};
		float[] fractionalVolumeComposition = {1.0f, 0, 0};
		// create the biomass species
		BiomassSpecies heterotroph = new BiomassSpecies("heterotroph", spa,
				fractionalVolumeComposition);
		heterotroph.setActiveMass(activeMassH);
		heterotroph.setEpsMass(eps);
		heterotroph.setInertMass(inert);
		// create the reactions
		//growth
		Reaction growth = new Reaction("growth", activeMassH,
				uMax, 1);
		growth.addFactor(new Saturation(oxygen, KO2));
		//eps decay
		Reaction epsDecay = new Reaction("epsDecay", eps,
				decayRateEPS, 1);
		epsDecay.addFactor(new Inhibition(oxygen, KO2));
		//create the maintenance reaction
		Reaction maintenance = new Reaction("maintenance",
				activeMassH, mO2, 1);
		maintenance.addFactor(new Saturation(oxygen, KO2));
		// inactivation
		Reaction inactivation = new Reaction("inactivation",
				activeMassH, inactivationConstantH, 1);
		inactivation.addFactor(new Inhibition(oxygen, KO2));
		// assign reactions to the biomasses through ReactionStoichiometry
		// active mass
		NetReaction rsb = new NetReaction(3);
		rsb.addReaction(growth, 1.0f);
		rsb.addReaction(maintenance, YH);
		rsb.addReaction(inactivation, -1.0f);
		activeMassH.setProcesses(rsb);
		// eps mass
		NetReaction rseps = new NetReaction(2);
		rseps.addReaction(growth, YEPS);
		rseps.addReaction(epsDecay, -1.0f);
		eps.setProcesses(rseps);
		// inert formation
		NetReaction rinert = new NetReaction(1);
		rinert.addReaction(inactivation, 1.0f);
		inert.setProcesses(rinert);
		// assign reaction stoichiometry to the chemical
		NetReaction rsc = new NetReaction(2);
		rsc.addReaction(growth, 1 / YH);
		rsc.addReaction(maintenance, 1.0f);
		oxygen.setProcesses(rsc);
		// add the species to system
		addBiomassSpecies(heterotroph);
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
		try {
			int[] nCells = {initialCellNumber};
			inoculateRandomly(nCells);
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
			return;
		}
	}
	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.ApplicationComponent#initializeDetachmentFunction()
	 */
	public void initializeDetachmentFunction() {
		DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}
	/**
	 * Simulation of 1000 iterative steps, storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		MultigridVariable.setSteps(10, 100);
		ApplicationComponent app = new DRGDetachmentEpsDecayApp();
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		VaribleSeries thickness = new BiofilmMaximumThicknessSeries();
		app = new SeriesVizualizer(app, thickness);
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		//app = new SeriesVizualizer(app, prod); // 
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		app = new SeriesVizualizer(app, biomass);
		// add vizualizer for solutes rates
		app = new SoluteRateSeriesVizualizer(app);
		// detached biomass
		app = new DetachedBiomassVizualizer(app);
		// finally, the controller must be the last decorator to add
		app = new VizualModelControler(app);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// set reactor dimensions
			// set the global mass balance parameters
			// --- nothing to set in this case: constant bulk concentration
			//initialize
			app.initializeSystemSpace();
			// construct directory name
			DecimalFormat df = new DecimalFormat("0.00E00");
			String dir = outputDirectory + "umax" + df.format(uMax) + "kdet"
					+ df.format(rdetach) + "kd_eps" + df.format(decayRateEPS)
					+ "ki" + df.format(inactivationConstantH) + "\\";
			app.intializeStateWriters(dir);
			//app.addStateWritter(new PovRayWriter());
			app.addStateWriter(new SoluteConcentrationWriter());
			app.addStateWriter(new SolidsConcentrationWriter());
			app.addStateWriter(new ParticlePositionWriter());
			//app.addStateWritter(new DetachmentLevelSetWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.detachedBiomassContainer().getTotalDetachedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			// write the header in parameters.txt file
			spw.initializeParametersWriting();
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
		}
		try {
			//wait for user to press start iteration
			app.waitForStartIteratingRequest();
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		//Program exits
		System.out.println("done.");
	}
}
