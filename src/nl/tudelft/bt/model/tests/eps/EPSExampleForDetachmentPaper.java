package nl.tudelft.bt.model.tests.eps;

import java.awt.Color;
import java.util.Iterator;
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
 * Example for model description paper
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class EPSExampleForDetachmentPaper extends ModelHandler {
	//output directory name
	private static String outputDirectory = "/Users/jxavier/results/";

	// geometry (2D or 3D)
	private static int geometry = 2;

	//
	//Solute species
	//Oxygen
	private static float oxygenBulkConcentration = 4e-3f; //[gO/l]

	private static float oxygenDiffusivity = 8.33e6f; //[um2/h]

	private static float KO = 3.5e-4f; //[gO/l]

	//private static float KOstar = 7e-4f; //[gO/l]

	//Substrate
	private static float substrateBulkConcentration = 0.1f; //[gCOD_S/l]

	private static float substrateDiffusivity = 4.17e6f; //[um2/h]

	private static float KS = 4e-3f; //[gCOD_S/l]

	//
	//Particulate species
	//Non PHB producing heterotroph H2
	private static float densityH = 200.0f; //[gCOD_X/l]

	//EPS
	private static float densityEPS = densityH / 6; //[gCOD_EPS/l]

	//Inert
	private static float densityI = densityH; //[gCOD_I/l]

	//
	//Yield coefficients
	private static float YSX = 0.2063f; //[gCOD_X/gCOS_S]

	private static float YSP = 0.2888f; //[gCOD_EPS/gCOS_S]

	//
	// Processes
	//substrate uptake
	private static float qMax = 0.952f; //[gCOD_S/gCOD_X/h]

	//Inert formation
	private static float kDecay = 0.0033f; //[gCOD_X/gCOD_X/h]

	//EPS hydrolysis rate coefficient
	//private static float kEps = 0.15f; //[gCOD_EPS/gCOD_EPS/h]
	private static float kEps = 0.014f; //[gCOD_EPS/gCOD_EPS/h]

	//private static float kEps = 0f; //[gCOD_EPS/gCOD_EPS/h]

	// Computation parameters
	private static float systemSize = 2000; // [um]

	private static float relativeMaximumRadius = 0.003f;

	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	private static float relativeBoundaryLayer = 0.1f;

	// other model parameters
	private static int gridSide = 65; // multigrid grid side

	private static float kShov = 1.0f; // shoving parameter[dim/less]

	private static float rdetach = 1e-3f;

	private static float simulationEnd = 2500f;
	private static float timeStep = 20.0f;
	
	

	// detachment constant[g/l/h]
	private static int initialCellNumber = 200;

	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// create the solutes
		//substrate
		SoluteSpecies substrate = new SoluteSpecies("substrate",
				substrateDiffusivity);
		substrate.setBulkConcentration(new ConstantBulkConcentration(
				substrateBulkConcentration));
		//oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		// create the particulates
		//H active mass
		ParticulateSpecies activeMass = new ParticulateSpecies("activeMass",
				densityH, Color.blue);
		//EPS
		ParticulateSpecies eps = new ParticulateSpecies("eps", densityEPS,
				Color.yellow);
		//Inert
		ParticulateSpecies inert = new ParticulateSpecies("inert", densityI,
				Color.red);
		// array of fixed species that constitute H
		ParticulateSpecies[] spH = { activeMass, eps, inert };
		float[] fractionalVolumeComposition = { 1.0f, 0, 0 };
		// create the biomass species
		BiomassSpecies heterotroph = new BiomassSpecies("heterotroph", spH,
				fractionalVolumeComposition);
		heterotroph.setActiveMass(activeMass);
		heterotroph.setInertMass(inert);
		heterotroph.setEpsMass(eps);
		//Create the Reaction factors, Monod and inhibition coefficients
		ProcessFactor mS = new Saturation(substrate, KS);
		ProcessFactor mO = new Saturation(oxygen, KO);
		//		ProcessFactor mOstar = new Saturation(oxygen, KOstar);
		// create the reactions
		//substrate uptake
		Reaction substrateUptake = new Reaction("substrateUptake", activeMass,
				qMax, 2);
		substrateUptake.addFactor(mS);
		substrateUptake.addFactor(mO);
		// growth
		//		Reaction growth = new Reaction("growth", activeMass,
		//				qMax*YSX, 1);
		//		growth.addFactor(mOstar);
		//decay
		Reaction decay = new Reaction("decayH2", activeMass, kDecay, 0);
		//hydrolysis of EPS
		Reaction hydrolysisEPS = new Reaction("hydrolysisEPS", eps, kEps, 0);
		//create the stoichiometries and add to species
		//H2 active mass
		NetReaction rsActive = new NetReaction(2);
		rsActive.addReaction(substrateUptake, YSX);
		rsActive.addReaction(decay, -1);
		activeMass.setProcesses(rsActive);
		//EPS
		NetReaction rsEps = new NetReaction(2);
		rsEps.addReaction(substrateUptake, YSP);
		//rsEps.addReaction(growth, -YSP/YSX);
		rsEps.addReaction(hydrolysisEPS, -1);
		eps.setProcesses(rsEps);
		//Inert
		NetReaction rsInert = new NetReaction(1);
		rsInert.addReaction(decay, 0.4f);
		inert.setProcesses(rsInert);
		// assign reaction stoichiometry to the solutes
		//substrate
		NetReaction rsSubstrate = new NetReaction(3);
		rsSubstrate.addReaction(substrateUptake, -1);
		rsSubstrate.addReaction(decay, 0.6f);
		rsSubstrate.addReaction(hydrolysisEPS, 1);
		substrate.setProcesses(rsSubstrate);
		//oxygen
		NetReaction rsOxygen = new NetReaction(1);
		rsOxygen.addReaction(substrateUptake, -(1 - YSX - YSP));
		//rsOxygen.addReaction(growth, -(YSP/YSX - 1));
		oxygen.setProcesses(rsOxygen);
		// add the species to system
		addBiomassSpecies(heterotroph);
		addSoluteSpecies(substrate);
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
			int[] nCells = { initialCellNumber };
			//int[] nCells = {initialCellNumber, initialCellNumber};
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
		// TODO set dtachment function here
		DetachmentSpeedFunction df = new Height2MassDetachment(rdetach);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation of 1000 iterative steps, storing results at each iteration
	 * 
	 * @param args -
	 *            argument should be the substrate bulk concentration it will be
	 *            checked if it is a number in the range 0.1 to 0.0008
	 */
	public static void main(String[] args) {
		//read the substrate concentration from command line
		//check argument
		if (args.length == 1) {
			substrateBulkConcentration = Float.parseFloat(args[0]);
			if ((substrateBulkConcentration > 0.1f)
					| (substrateBulkConcentration < 0.0008f))
				throw new ModelRuntimeException("concentration "
						+ substrateBulkConcentration + " in the wrong range");
			outputDirectory += "./C_S" + args[0];
			System.out.println("Simulation started with C_S = "
					+ substrateBulkConcentration);
		} else {
			throw new ModelRuntimeException("A value for the substrate"
					+ " concentration must be provided in the command line");
		}
		MultigridVariable.setSteps(50, 500);
		ApplicationComponent app = new EPSExampleForDetachmentPaper();
		// set the end of simulation
		Model.model().setFinishIterationTime(simulationEnd);
		// set compulsory time step
		Model.model().setCompulsoryTimeStep(timeStep);
		//
		app = new BiomassVizualizer(app);
		// the biomass thickness visualizer
		VaribleSeries thickness = new BiofilmMaximumThicknessSeries();
		app = new SeriesVizualizer(app, thickness);
		// the produced biomass
		ProducedBiomassSeries prod = new ProducedBiomassSeries();
		//app = new SeriesVizualizer(app, prod);
		// the biofilm total biomass
		FixedTotalBiomassSeries biomass = new FixedTotalBiomassSeries();
		//app = new SeriesVizualizer(app, biomass);
		// add vizualizer for solutes rates
		app = new SoluteRateSeriesVizualizer(app);
		// detached biomass
		app = new DetachedBiomassVizualizer(app);
		// bulk concentrations
		app = new BulkConcentrationVizualizer(app);
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
			app.intializeStateWriters(outputDirectory);
			//app.addStateWriter(new PovRayWriter());
			app.addStateWriter(new SoluteConcentrationWriter());
			app.addStateWriter(new SolidsConcentrationWriter());
			app.addStateWriter(new ParticlePositionWriter());
			// the simulation parameters writter
			SimulationResultsWriter spw = new SimulationResultsWriter();
			spw.addSeries(thickness);
			spw.addSeries(Model.detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
			spw.addSeries(Model.detachedBiomassContainer()
					.getErodedBiomassSeries());
			spw.addSeries(Model.detachedBiomassContainer()
					.getSloughedBiomassSeries());
			spw.addSeries(prod);
			spw.addSeries(biomass);
			app.addStateWriter(spw);
			// initialize
			app.initializeDiffusionReactionSystem(); // also innoculates
			//
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
			// write the header in parameters.txt file
			spw.initializeParametersWriting();
		} catch (ModelException e) {
			System.out.println(e);
			System.exit(-1);
		}
		try {
			//wait for user to press start iteration
			//app.waitForStartIteratingRequest();
			// start iterating cycle
			app.startIterating();
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		//terminate program
		System.out.println("done.");
		System.exit(1);
	}
}