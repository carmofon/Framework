package nl.tudelft.bt.model.tests.detachmentstudy;

import java.awt.Color;
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
 * instance)
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class DiffusionReactionGrowthApp extends ModelHandler {
	// output directory name
	private final static String outputDirectory = "/Users/jxavier/lixo/";

	// geometry (2D or 3D)
	private static int geometry = 2;

	// Diffusion to reaction:
	private static final float diffusivity = 5.76e6f; // [um2/h]

	private static float uMax = 5.47e-1f; // default val 5.47e-2f

	// change to 5.47f to observe
	// finger formation [1/h]
	private static float systemSize = 500; // [um]

	// carbon source to biomass convertion:
	private static float density = 70.0f; // [g/l = 10^-15g/um^3]

	// change to 150.0f to observe
	// finger formation
	private static final float YS = -0.045f;

	private static float constantBulkConcentration = 4e-3f; // [2.4f g/l]

	private static final float KS = 3.5e-4f; // [gCOD/l]

	// maintenance coefficient
	private static final float mS = -1.08e-1f; // [KgS/kgX/h]

	// Dimensionless numbers
	private static float relativeMaximumRadius = 0.008f;

	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	private static final float minimumMassRatio = 0.001f;

	private static float relativeBoundaryLayer = 0.02f;

	// other model parameters
	private static int gridSide = 33; // multigrid grid side

	private static float kShov = 1.0f; // shoving parameter[dim/less]

	private static float rdetach = 1.0e-3f;

	// detachment constant[g/l/h]
	private static int initialCellNumber = 15;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// create the chemicals
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", diffusivity);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				constantBulkConcentration));
		// create the species
		ParticulateSpecies activeH = new ParticulateSpecies(
				"heterotrophActiveMass", density, Color.gray);
		// array of fixed species that constitute the denitrifier
		ParticulateSpecies[] spa = { activeH };
		float[] fractionalVolumeComposition = { 1.0f };
		// create the biomass species
		BiomassSpecies heterotroph = new BiomassSpecies("heterotroph", spa,
				fractionalVolumeComposition);
		// create the reaction
		Reaction growth = new Reaction("growth", activeH, uMax, 1);
		growth.addFactor(new Saturation(oxygen, KS));
		// create the maintenance reaction
		Reaction maintenance = new Reaction("maintenance", activeH, mS, 1);
		maintenance.addFactor(new Saturation(oxygen, KS));
		// assign reaction to the species through a ReactionStoichiometry
		NetReaction rsb = new NetReaction(2);
		rsb.addReaction(growth, 1.0f);
		rsb.addReaction(maintenance, YS);
		activeH.setProcesses(rsb);
		// assign reaction stoichiometry to the chemical
		// TODO put factor back to 2 to consider maintenance
		NetReaction rsc = new NetReaction(2);
		rsc.addReaction(growth, 1 / YS);
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
		int[] nCells = { initialCellNumber };
		inoculateRandomly(nCells);
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
		ApplicationComponent app = new DiffusionReactionGrowthApp();
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		VaribleSeries thickness = new BiofilmMaximumThicknessSeries();
		app = new SeriesVizualizer(app, thickness);
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// app = new SeriesVizualizer(app, prod);
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// app = new SeriesVizualizer(app, biomass);
		// add vizualizer for solutes rates
		app = new SoluteRateSeriesVizualizer(app);
		// detached biomass
		// app = new DetachedBiomassVizualizer(app);
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
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			// app.addStateWritter(new PovRayWriter());
			app.addStateWriter(new SoluteConcentrationWriter());
			app.addStateWriter(new SolidsConcentrationWriter());
			app.addStateWriter(new ParticlePositionWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
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