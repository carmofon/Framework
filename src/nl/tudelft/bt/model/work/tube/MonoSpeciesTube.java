package nl.tudelft.bt.model.work.tube;

import java.awt.Color;
import java.util.Iterator;
import nl.tudelft.bt.model.*;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.components.*;
import nl.tudelft.bt.model.apps.output.*;
import nl.tudelft.bt.model.bulkconcentrations.*;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.detachment.levelset.functions.RadialProfileMassDetachment;
import nl.tudelft.bt.model.detachment.levelset.functions.TubeRadius2MassDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.particlebased.tube.TubeModelHandler;
import nl.tudelft.bt.model.particlebased.tube.detachment.DetachmentInducedByShear;
import nl.tudelft.bt.model.particlebased.tube.detachment.ParticleBasedDetachmentInducedByShear;
import nl.tudelft.bt.model.reaction.*;

/**
 * Growth of a monospecies biofilm with a single substrate species. No
 * detachment forces present and growth is represented in 2D space. The biofilm
 * developed shows a rough morphology (heterogeneous shape) as a consequance of
 * high diffusion limitation. The effect of diffusion limitation is represented
 * by the dimesionless G number, defined as
 * 
 * Note: um representes micro-meter
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class MonoSpeciesTube extends TubeModelHandler {
	// All the model parameters are defined here as static attrributes
	// at the begining of the Class. This way, they can be easily changed
	// without changing the remaining program

	public static String outputDirectory;

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
	// Substrate (S) - the only solute species used here
	protected static float substrateBulkConcentration = 10e-3f; // [g/L]

	// in the paper Ds = 1.6e-9 m2/s = 5.76e6 um2/h
	private static float substrateDiffusivity = 5.76e6f; // [um^2/h]

	// 

	//
	// Particulate species (biomass X)
	protected static float specificMassX = 70f; // [g/L]

	// Yield coefficients
	private static float YXS = 0.45f; // [gX/gS]

	// Processes
	// Growth (biomass production)
	protected static float uMax = 0.547f; // [1/h]

	protected static float maintenanceCoefficient = uMax * 0f; // [1/h]

	private static float KS = 3.5e-4f; // [g/L]

	// Computation parameters
	protected static int gridSide = 129; // multigrid grid side

	// protected static int gridSide = 65; // multigrid grid side

	protected static float kShov = 1.0f; // shoving parameter[dim/less]

	//
	// protected static float tubeRadius = 13000f; //[um]
	protected static float tubeRadius = 3500f; // [um]

	protected static float tubeLength = 1e6f; // [um]

	// protected static float flow = 1.021e-5f * 1e18f * 3600f; //[um3/h]
	protected static float flow = 3e14f; // [um3/h]

	protected static float viscosity = 1.0e-3f * 1e3f * 1e9f * 3600f; // [1e-15g/um/h]

	//
	protected static float systemSize = tubeRadius
			/ (0.5f - 1f / (float) gridSide) * 1.01f; // [um]

	// relativeMaximumRadius defines the maximum radius of the biomass particles
	// in relation to the system size
	// the maximum radius of a particle is rmax =
	// systemSize*relativeMaximumRadius
	protected static float relativeMaximumRadius = 0.003f;

	// Similarly to relativeMaximumRadius, relativeMinimumRadius defines the
	// minimum radius of a particle in the system
	protected static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	// Defines the thickness of the concentration boundary layer in the system.
	// Here, the thickness of the boundary layer is 0.1*2000 = 200 um
	protected static float relativeBoundaryLayer = 0.1f;

	// other model parameters

	// initial number of particles in the system (inoculum)
	protected static int initialParticleNumber = 200;

	protected static float compulsoryTimeStep = 20f;

	// protected static float compulsoryTimeStep = 2f; //TODO

	protected static float finishIterationTime = 10000f;

	// protected static float finishIterationTime = 370f; //TODO

	//protected float _detachmentShearStress = 5e23f;

	//protected float _detachmentShearStress = 5e16f;

	protected float _detachmentShearStress = 2e17f;

	// protected float _detachmentShearStress = 0f;

	//

	/**
	 * Constructor
	 */
	public MonoSpeciesTube() {
		super(tubeRadius, tubeLength, flow, viscosity);
	}

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// 1. Create the solutes
		// substrate
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		// set up the simplest type of bulk concentration: constant
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));
		// 2. Create the particulate species (soliids)
		// X active mass
		ParticulateSpecies activeX = new ParticulateSpecies("activeX",
				specificMassX, Color.blue);
		ParticulateSpecies inertX = new ParticulateSpecies("inertX",
				specificMassX, Color.gray);
		// array of fixed species that constitute speciesX (in this case,
		// speciesX is entirely constituted by active mass)
		ParticulateSpecies[] spX = { activeX, inertX };
		float[] fractionalVolumeCompositionH1 = { 1.0f, 0f };
		// 3. Create the biomass species
		BiomassSpecies speciesX = new BiomassSpecies("speciesX", spX,
				fractionalVolumeCompositionH1);
		speciesX.setActiveMass(activeX);
		speciesX.setInertMass(inertX);
		// speciesX.getColorFromGrowth();
		// 4. Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mS = new Saturation(substrate, KS);
		// The Saturation class creates a process factor with the form
		// Cs/(Cs+KS) where Cs is the concentration of substrate
		//
		// 5. Create the reactions
		// growth
		Reaction growth = new Reaction("growth", activeX, uMax, 1);
		growth.addFactor(mS);
		// This creates a growth rate that equals:
		// rX = uMax*Cs/(Cs+KS)*Cx
		// where Cx is the concentration of biomass
		Reaction decay = new Reaction("decay", activeX, maintenanceCoefficient,
				0);
		// 6. Assign reaction to the species through ReactionStoichiometries
		// active mass
		NetReaction rsXactive = new NetReaction(2);
		rsXactive.addReaction(growth, 1);
		rsXactive.addReaction(decay, -1);
		activeX.setProcesses(rsXactive);
		// This defines that biomass growth rate is 1*rX
		NetReaction rsXinert = new NetReaction(1);
		rsXinert.addReaction(decay, 1);
		inertX.setProcesses(rsXinert);
		//
		// assign reaction stoichiometry to the solutes
		// substrate
		NetReaction rsSubstrate = new NetReaction(1);
		rsSubstrate.addReaction(growth, -(1 / YXS));
		substrate.setProcesses(rsSubstrate);
		// This defines that substrate consumption rate is -(1 / YXS)*rX
		//
		// 7. add the solute species and the biomass species (which contain the
		// particulate species) to system
		addBiomassSpecies(speciesX);
		addSoluteSpecies(substrate);
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
		int[] nCells = { initialParticleNumber };
		inoculateRandomly(nCells);
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeDetachmentFunction() {
		// The detachment function is set here. However, in this case,
		// detachment is not considered since rdetach = 0
		// DetachmentSpeedFunction df = new TubeRadius2MassDetachment(1e-2f);
		// DetachmentSpeedFunction df = new TubeRadius2MassDetachment(5e-4f);
		// DetachmentSpeedFunction df = new TubeRadius2MassDetachment(0);
		ParticleBasedDetachmentInducedByShear df = new ParticleBasedDetachmentInducedByShear(
				_detachmentShearStress);
		// DetachmentInducedByShear df = new DetachmentInducedByShear(
		// _detachmentShearStress);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		// read output directory from the command line, unless
		// program is running localy (in my laptop)
		boolean runWithGraphics = false;
		if (args.length != 2) {
			throw new RuntimeException("Program arguments missing: "
					+ "2 program arguments should be supplied"
					+ " (1 - output directory,"
					+ " 2 - flag for running with graphics)");
		}
		// parse the input arguments
		outputDirectory = args[0];
		int arg1 = Integer.parseInt(args[1]);
		switch (arg1) {
		case 0:
			runWithGraphics = false;
			break;
		case 1:
			runWithGraphics = true;
			break;
		default:
			throw new RuntimeException("second program" + " argument must be 0"
					+ " (for running with no graphics) "
					+ "or 1 (for running with graphics)");
		}
		// //start
		MultigridVariable.setSteps(5, 50);
		// create a hande for the application, which will be decorated
		ApplicationComponent app = new MonoSpeciesTube();
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		// the thickness series
		VariableSeries thickness = new BiofilmMaximumThicknessSeries();
		// The following code will be omitted if no vizuals are desired
		if (runWithGraphics) {
			// The following code will be omitted if no vizuals are desired
			// start decorationg the application
			app = new BiomassVizualizer(app);
			// the total biomass visualizer
			app = new SeriesVizualizer(app, biomass);
			// detached biomass
			app = new DetachedBiomassVizualizer(app);
			// add vizualizer for solutes rates
			// app = new SoluteRateSeriesVizualizer(app);
			// add the vontroler
			app = new VizualModelControler(app);
		}
		try {
			// create the space
			app.setSystemSpaceParameters(geometry, systemSize,
					relativeMaximumRadius, relativeMinimumRadius,
					relativeBoundaryLayer, gridSide, kShov);
			// --- nothing to set in this case: constant bulk concentration
			// initialize
			app.initializeSystemSpace();
			app.intializeStateWriters(outputDirectory);
			app.addTimedStateWriter(new PovRayWriter());
			app.addTimedStateWriter(new SoluteConcentrationWriter());
			app.addTimedStateWriter(new SolidsConcentrationWriter());
			app.addTimedStateWriter(new ParticlePositionWriter());
			app.addTimedStateWriter(new FlowWriter());
			app.addTimedStateWriter(new TubeDetachmentInfoWriter());
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
			// add the time constraints writer
			app.addStateWriter(new TimeConstraintsWriter());
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			// and creates the biomass particle container
			app.initializeDetachmentFunction();
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
			// add the velocity profile to write
			app.addTimedStateWriter(new ProfileWriter(TubeVelocityField
					.getInstance().getVelocityRadialProfile()));
			app.addTimedStateWriter(new ProfileWriter(TubeVelocityField
					.getInstance().getShearRadialProfile()));
			app.writeState();
			//
			Model.model().setCompulsoryTimeStep(compulsoryTimeStep);
			Model.model().setFinishIterationTime(finishIterationTime);
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
		}
		try {
			// app.waitForStartIteratingRequest();
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		System.out.println("Simulation finished.");
	}
}