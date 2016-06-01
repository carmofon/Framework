package nl.tudelft.bt.model.tests.monospecies;

import java.awt.Color;

import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.HeightMassDetachment;
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
public class MonoSpeciesDRGDetachment extends ModelHandler {
	// output directory name
	private final static String outputDirectory = "/Users/jxavier/lixo/";

	// geometry (2D or 3D)
	private static int geometry = 2;

	// Diffusion to reaction:
	private static final float diffusivity = 2e6f; // [um2/h]

	private static float uMax = 4e-1f; // default val 5.47e-2f

	// change to 5.47f to observe
	// finger formation [1/h]
	private static float systemSize = 10000; // [um]

	// carbon source to biomass convertion:
	private static float density = 70.0f; // [g/l = 10^-15g/um^3]

	// change to 150.0f to observe
	// finger formation
	private static final float YS = -0.5f;

	private static float inputConcentration = 0.1f; // [2.4f g/l]

	private static final float KS = 3.5e-4f; // [gCOD/l]

	// Dimensionless numbers
	private static float relativeMaximumRadius = 0.008f;

	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	private static final float minimumMassRatio = 0.001f;

	private static float relativeBoundaryLayer = 0.1f;

	// other model parameters
	private static int gridSide = 33; // multigrid grid side

	private static float kShov = 1.0f; // shoving parameter[dim/less]

	private static float rdetach = 50f / 700; // detachment rate
												// [10^-15g/um^2/h]

	private static int initialCellNumber = 24;

	// reactor parameters
	private static float reactorVolume = 2e15f; // [um^3]

	private static float carrierArea = 9.59e13f; // [9.59e10f um^2]

	private static float residenceTime = 0.6696f; // [um^2]

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// create the chemicals
		SoluteSpecies glucose = new SoluteSpecies("glucose", diffusivity);
		glucose.setBulkConcentration(new DynamicBulkConcentrationExplicit(
				inputConcentration));
		// create the species
		ParticulateSpecies sp = new ParticulateSpecies("heterotrophActiveMass",
				density, Color.gray);
		// array of fixed species that constitute the denitrifier
		ParticulateSpecies[] spa = { sp };
		float[] fractionalVolumeComposition = { 1.0f };
		// create the biomass species
		BiomassSpecies heterotroph = new BiomassSpecies("heterotroph", spa,
				fractionalVolumeComposition);
		// create the reaction
		Reaction growth = new Reaction("growth", sp, uMax, 1);
		growth.addFactor(new Saturation(glucose, KS));
		// assign reaction to the species through a ReactionStoichiometry
		NetReaction rsb = new NetReaction(1);
		rsb.addReaction(growth, 1.0f);
		sp.setProcesses(rsb);
		// assign reaction stoichiometry to the chemical
		NetReaction rsc = new NetReaction(1);
		rsc.addReaction(growth, YS);
		glucose.setProcesses(rsc);
		// add the species to system
		addBiomassSpecies(heterotroph);
		addSoluteSpecies(glucose);
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
		int[] nCells = { initialCellNumber };
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.ApplicationComponent#initializeDetachmentFunction()
	 */
	public void initializeDetachmentFunction() {
		DetachmentSpeedFunction df = new HeightMassDetachment(rdetach);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation of 1000 iterative steps, storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		MultigridVariable.setSteps(5, 50);
		ApplicationComponent app = new MonoSpeciesDRGDetachment();
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
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			app.addStateWriter(new PovRayWriter());
			app.addStateWriter(new SoluteConcentrationWriter());
			app.addStateWriter(new SolidsConcentrationWriter());
			app.addStateWriter(new ParticlePositionWriter());
			app.addStateWriter(new SimulationResultsWriter());
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
		}
		try {
			// wait for user to press start iteration
			app.waitForStartIteratingRequest();
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		// TODO remove this line
		System.out.println("done.");
	}

}
