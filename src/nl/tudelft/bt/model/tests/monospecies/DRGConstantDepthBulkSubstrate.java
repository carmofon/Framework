package nl.tudelft.bt.model.tests.monospecies;
import java.awt.Color;
import java.util.Iterator;
import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.ModelHandler;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
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
public class DRGConstantDepthBulkSubstrate extends ModelHandler {
	//output directory name
	private final static String outputDirectory = "D:\\results\\1Dcomparison"
			+ "\\monospecies\\bulk2Dimplicit\\";
	// geometry (2D or 3D)
	private static int geometry = 2;
	//
	// Solute species
	// O2
	private static final float diffusivityO2 = 8333333f; //[um2/h]
	private static final float KO2 = 2e-4f; //[gO2/l]
	private static float bulkConcentrationO2 = 5e-3f; //[2.4f g/l]
	// S
	private static final float diffusivityS = 4166666f; //[um2/h]
	private static final float KS = 4e-3f; //[gCOD/l]
	private static float inputConcentrationS = 0.03f; //[gCOD/l]
	//
	// Particulate species
	private static float uMax = 0.2499f; // default val 5.47e-2f
	private static float density = 10.0f; // [g/l = 10^-15g/um^3]
	private static final float YS = 0.63f;
	// change to 5.47f to observe
	// finger formation [1/h]
	//Computational parameters
	private static float systemSize = 2000; // [um]
	private static float relativeMaximumRadius = 0.008f;
	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;
	private static final float minimumMassRatio = 0.001f;
	private static int gridSide = 33; // multigrid grid side
	private static float kShov = 1.0f; // shoving parameter[dim/less]
	// reactor operation parameters
	private static float relativeBoundaryLayer = 0f;
	private static float maximumThickness = 545.4f; // [um]
	private static int initialCellNumber = 24;
	private static float reactorVolume = 1.25e15f; //[um^3]
	private static float carrierArea = 1E+11f; //[um^2]
	private static float residenceTime = 1.5f; //[h]
	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// create the solutes
		// oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", diffusivityO2);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				bulkConcentrationO2));
		// carbon source
		SoluteSpecies carbonSource = new SoluteSpecies("carbonSource",
				diffusivityS);
		carbonSource.setBulkConcentration(new DynamicBulkConcentrationImplicit(
				inputConcentrationS));
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
		growth.addFactor(new MonodTwoSubstrates(oxygen, KO2, carbonSource, KS));
		// assign reaction to the species through a ReactionStoichiometry
		NetReaction rsH = new NetReaction(1);
		rsH.addReaction(growth, 1.0f);
		sp.setProcesses(rsH);
		// assign reaction stoichiometry to the solutes
		// oxygen
		NetReaction rsO2 = new NetReaction(1);
		rsO2.addReaction(growth, -(1 - YS) / YS);
		oxygen.setProcesses(rsO2);
		// carbon source
		NetReaction rsS = new NetReaction(1);
		rsS.addReaction(growth, -1 / YS);
		carbonSource.setProcesses(rsS);
		// add the species to system
		addBiomassSpecies(heterotroph);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(carbonSource);
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
		MultigridVariable.setSteps(100, 1000);
		ApplicationComponent app = new DRGConstantDepthBulkSubstrate();
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
		// bulk concentrations
		app = new BulkConcentrationVizualizer(app);
		// finally, the controller must be the last decorator to add
		app = new VizualModelControler(app);
		float dX = systemSize * (1 - relativeBoundaryLayer);
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// set reactor dimensions
			// set the global mass balance parameters
			app.setReactorParameters(residenceTime, carrierArea, reactorVolume);
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
			// add bulk concentrations of all solutes as variable series
			//
			app.addStateWriter(spw);
			//
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
			app.initializeDetachmentFunction();
			for (Iterator iter = Model.model().getSoluteSpecies().iterator(); iter
					.hasNext();) {
				SoluteSpecies s = (SoluteSpecies) iter.next();
				spw.addSeries(s.getBulkConcentrationSeries());
			}
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
