package nl.tudelft.bt.model.tests.monospecies;
import java.awt.Color;
import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;
import nl.tudelft.bt.model.util.ExtraMath;
/**
 * Simple 2D monospecies system, can be used to study the effect of diffusion
 * limitation on biofilm structure (change growth rate and density, for
 * instance)
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class MonoSpeciesDRGConstantDepth extends ModelHandler {
	//output directory name
	private final static String outputDirectory = "D:\\results\\1Dcomparison"
			+ "\\monospecies\\C_02_5e-3\\";
	// geometry (2D or 3D)
	private static int geometry = 2;
	// Diffusion to reaction:
	private static final float diffusivity = 8333333f; //[um2/h]
	private static float uMax = 0.2499f; // default val 5.47e-2f
	// change to 5.47f to observe
	// finger formation [1/h]
	private static float systemSize = 2000; // [um]
	// carbon source to biomass convertion:
	private static float density = 10.0f; // [g/l = 10^-15g/um^3]
	// change to 150.0f to observe
	// finger formation
	private static final float YS = 0.63f;
	//private static float inputConcentration = 1e-1f; //[2.4f g/l]
	private static float inputConcentration = 5e-3f; //[2.4f g/l]
	private static final float KS = 2e-4f; //[gCOD/l]
	//Dimensionless numbers
	private static float relativeMaximumRadius = 0.008f;
	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;
	private static final float minimumMassRatio = 0.001f;
	private static float relativeBoundaryLayer = 0f;
	// other model parameters
	private static int gridSide = 33; // multigrid grid side
	private static float kShov = 1.0f; // shoving parameter[dim/less]
	private static float maximumThickness = 540; // [um]
	private static int initialCellNumber = 24;
	//reactor parameters
	private static float reactorVolume = 2e15f; //[um^3]
	private static float carrierArea = 9.59e13f; //[9.59e10f um^2]
	private static float residenceTime = 0.6696f; //[um^2]
	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// create the chemicals
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", diffusivity);
		oxygen.setBulkConcentration(new DynamicBulkConcentrationExplicit(
				inputConcentration));
		// create the species
		ParticulateSpecies sp = new ParticulateSpecies("heterotrophActiveMass", density, Color.gray);
		// array of fixed species that constitute the denitrifier
		ParticulateSpecies[] spa = {sp};
		float[] fractionalVolumeComposition = {1.0f};
		// create the biomass species
		BiomassSpecies heterotroph = new BiomassSpecies("heterotroph", spa,
				fractionalVolumeComposition);
		// create the reaction
		Reaction growth = new Reaction("growth", sp, uMax, 1);
		growth.addFactor(new Saturation(oxygen, KS));
		// assign reaction to the species through a ReactionStoichiometry
		NetReaction rsb = new NetReaction(1);
		rsb.addReaction(growth, 1.0f);
		sp.setProcesses(rsb);
		// assign reaction stoichiometry to the chemical
		NetReaction rsc = new NetReaction(1);
		rsc.addReaction(growth, -(1-YS)/YS);
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
		// To use a detachment function set here
		//DetachmentSpeedFunction df =
		//	new HeightMassDetachment(rdetach);
		//setDetachmentFuntion(df);
		try {
			Model.model().setMaximumBiofilmHeight(maximumThickness);
		} catch (InvalidValueException e) {
			throw new ModelRuntimeException(e.toString());
		}
	}
	/**
	 * Simulation of 1000 iterative steps, storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		MultigridVariable.setSteps(5, 50);
		ApplicationComponent app = new MonoSpeciesDRGConstantDepth();
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		VaribleSeries thickness = new BiofilmMaximumThicknessSeries();
		app = new SeriesVizualizer(app, thickness);
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		app = new SeriesVizualizer(app, biomass);
		// add vizualizer for solutes rates
		app = new SoluteRateSeriesVizualizer(app);
		// detached biomass
		app = new DetachedBiomassVizualizer(app);
		// finally, the controller must be the last decorator to add
		app = new VizualModelControler(app);
		float dX = systemSize * (1 - relativeBoundaryLayer);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// set reactor dimensions
			// set the global mass balnance parameters
			float computationalVolumeRatio = carrierArea
					/ ExtraMath.sq(systemSize);
			app.setReactorParameters(residenceTime, computationalVolumeRatio,
					reactorVolume);
			//initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			//app.addStateWritter(new PovRayWriter());
			app.addStateWriter(new SoluteConcentrationWriter());
			app.addStateWriter(new SolidsConcentrationWriter());
			app.addStateWriter(new ParticlePositionWriter());
			//app.addStateWritter(new SimulationParametersWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.detachedBiomassContainer().getTotalDetachedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			//
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
		//TODO remove this line
		System.out.println("done.");
	}
}
