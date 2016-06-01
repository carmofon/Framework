package nl.tudelft.bt.model.tests.phb;
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
 * Simple 2D monospecies system, can be used to study the effect of diffusion
 * limitation on biofilm structure (change growth rate and density, for
 * instance)
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class PHBProducerEPSProducerSameUptake extends ModelHandler {
	//output directory name
	private static String outputDirectory = "E:/results/eps/lixo/";
	//	private static String outputDirectory = "D:\\results\\phb"
	//+ "\\resultsForModelDescriptionPaper\\feastFamine\\";
	// geometry (2D or 3D)
	private static int geometry = 2;
	//
	//Solute species
	//Acetate
	private static float acetateBulkConcentration = 0.1f; //[gCOD_S/l]
	private static float carbonSourceDiffusivity = 4e6f; //[um2/h]
	private static float KHAc = 4e-3f; //[gO/l]
	private static float feastTime = 4f; //[h]
	private static float famineTime = 24f - feastTime; //[h]
	private static float initialFeastTime = feastTime; //[h]
	//Oxygen
	private static float oxygenBulkConcentration = 4e-3f; //[gO/l]
	private static float oxygenDiffusivity = 8e6f; //[um2/h]
	private static float KO = 3.5e-4f; //[gO/l]
	//This value produces large amounts of polymer
	//private static float KO_star = 1.72E-02f; //[gO/l]
	private static float KO_star = 7E-04f; //[gO/l]
	//
	//Particulate species
	//Non PHB producing heterotroph H2
	private static float densityH = 200.0f; //[gCOD_X/l]
	//EPS
	private static float densityEPS = densityH / 6; //[gCOD_EPS/l]
	//PHB
	private static float densityPHB = densityH * 5; //[gCOD_PHB/l]
	//Inert
	private static float densityI = densityH; //[gCOD_I/l]
	//
	//Yield coefficients
	private static float YSX = 0.495f; //[gCOD_X/gCOS_S]
	private static float YSP = 0.75f; //[gCOD_PHB/gCOS_S]
	private static float YPX = 0.624f; //[gCOD_X/gCOS_PHB]
	//
	// Processes
	//Acetate uptake
	private static float qMax = 0.952f; //[gCOD_S/gCOD_X/h]
	//PHB consumption
	private static float kPhb = 0.15f; //[gCOD_PHB/gCOD_X/h]
	// Polymer inhibition
	private static float KP = 1f; //[fP]
	//Inert formation
	private static float kDecay = 0.0033f; //[gCOD_X/gCOD_X/h]
	// Computation parameters
	private static float systemSize = 1000; // [um]
	private static float relativeMaximumRadius = 0.003f;
	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;
	private static float relativeBoundaryLayer = 0.1f;
	// other model parameters
	private static int gridSide = 33; // multigrid grid side
	private static float kShov = 1.0f; // shoving parameter[dim/less]
	private static float rdetach = 3e-3f;
	// detachment constant[g/l/h]
	private static int initialCellNumber = 31;
	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		//           ---- Define the solutes ----
		//---carbon source
		SoluteSpecies carbonSource = new SoluteSpecies("carbonSource",
				carbonSourceDiffusivity);
		carbonSource.setBulkConcentration(new ConstantBulkConcentration(
				acetateBulkConcentration));
		//The following code is used for feast/famine cycles to
		//keep acetate concentration average in time:
		//float ac = acetateBulkConcentration
		//				/ (feastTime / (feastTime + famineTime));
		//acetate.setBulkConcentration(new ItermittentBulkConcentration(ac,
		//				initialFeastTime, feastTime, famineTime));
		//---oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenBulkConcentration));
		//           ---- Create the particulates ----
		//---H-PHB active mass
		ParticulateSpecies activeH1 = new ParticulateSpecies(
				"H-PHB-active-mass", densityH, Color.blue);
		//---H-EPS active mass
		ParticulateSpecies activeH2 = new ParticulateSpecies(
				"H-EPS-active-mass", densityH, Color.red);
		//---PHB
		ParticulateSpecies phb = new ParticulateSpecies("phb", densityPHB,
				Color.yellow);
		//---EPS
		ParticulateSpecies eps = new ParticulateSpecies("eps", densityEPS,
				Color.yellow);
		//---Inert
		//Inert is descriminated by species to assess mass balances
		ParticulateSpecies inertH1 = new ParticulateSpecies("inert-H-PHB",
				densityI, Color.gray);
		ParticulateSpecies inertH2 = new ParticulateSpecies("inert-H-EPS",
				densityI, Color.gray);
		//           ---- Create the biomass species ----
		//----PHB producer
		//array of fixed species that constitute H-PHB biomass species
		ParticulateSpecies[] spH1 = {activeH1, phb, inertH1};
		float[] fractionalVolumeCompositionH1 = {1.0f, 0, 0};
		BiomassSpecies speciesH1 = new BiomassSpecies("PHB producer", spH1,
				fractionalVolumeCompositionH1);
		speciesH1.setActiveMass(activeH1);
		speciesH1.setInertMass(inertH1);
		//----EPS producer
		//array of fixed species that constitute H-PHB biomass species
		ParticulateSpecies[] spH2 = {activeH2, eps, inertH2};
		float[] fractionalVolumeCompositionH2 = {1.0f, 0, 0};
		BiomassSpecies speciesH2 = new BiomassSpecies("EPS producer", spH2,
				fractionalVolumeCompositionH2);
		speciesH2.setActiveMass(activeH2);
		speciesH2.setInertMass(inertH2);
		speciesH2.setEpsMass(eps); //defines the species that is modeled as
								   // capsule
		//           ---- Create the process terms ----
		//Monod and inhibition coefficients
		ProcessFactor mHAc = new Saturation(carbonSource, KHAc);
		ProcessFactor iHAc = new Inhibition(carbonSource, KHAc);
		ProcessFactor mO = new Saturation(oxygen, KO);
		ProcessFactor mO_star = new Saturation(oxygen, KO_star);
		ProcessFactor iPHB = new InhibitionFromFraction(phb, activeH1, KP);
		ProcessFactor iEPS = new InhibitionFromFraction(eps, activeH2, KP);
		//           ---- Create the reactions ----
		//---carbon source uptake by H-PHB
		Reaction sUptakeH1 = new Reaction(
				"carbon source uptake by H-PHB", activeH1, qMax, 3);
		sUptakeH1.addFactor(mHAc);
		sUptakeH1.addFactor(mO);
		sUptakeH1.addFactor(iPHB);
		//---acetate uptake by H-EPS
		Reaction sUptakeH2 = new Reaction(
				"carbon source uptake by H-EPS", activeH2, qMax, 3);
		sUptakeH2.addFactor(mHAc);
		sUptakeH2.addFactor(mO);
		sUptakeH2.addFactor(iEPS);
		//---Growth of H-PHB
		Reaction growthH1 = new Reaction("growthH1",
				activeH1, qMax * YSX, 2);
		growthH1.addFactor(mHAc);
		growthH1.addFactor(mO_star);
		//---Growth of H-EPS
		Reaction growthH2 = new Reaction("growthH2",
				activeH2, qMax * YSX, 2);
		growthH2.addFactor(mHAc);
		growthH2.addFactor(mO_star);
		//---decay of H-PHB
		Reaction decayH1 = new Reaction("decayH1", activeH1,
				kDecay, 0);
		//---decay H-EPS
		Reaction decayH2 = new Reaction("decayH2", activeH2,
				kDecay, 0);
		//---PHB consumption
		Reaction phbConsumption = new Reaction(
				"phbConsumption", phb, kPhb, 2);
		phbConsumption.addFactor(iHAc);
		phbConsumption.addFactor(mO);
		// 		---- Assign reaction to each species through NetRaction ----
		//---H-PHB active mass
		NetReaction rsH1active = new NetReaction(3); // 3 is number of reactions
		rsH1active.addReaction(growthH1, 1);
		rsH1active.addReaction(phbConsumption, YPX);
		rsH1active.addReaction(decayH1, -1);
		activeH1.setProcesses(rsH1active);
		//---H-EPS active mass
		NetReaction rsH2active = new NetReaction(2);
		rsH2active.addReaction(growthH2, 1);
		rsH2active.addReaction(decayH2, -1);
		activeH2.setProcesses(rsH2active);
		//---PHB
		NetReaction rsPhb = new NetReaction(3);
		rsPhb.addReaction(sUptakeH1, YSP);
		rsPhb.addReaction(growthH1, -YSP / YPX);
		rsPhb.addReaction(phbConsumption, -1);
		phb.setProcesses(rsPhb);
		//---EPS
		NetReaction rsEps = new NetReaction(2);
		rsEps.addReaction(sUptakeH2, YSP);
		rsEps.addReaction(growthH2, -YSP / YPX);
		eps.setProcesses(rsEps);
		//---InertH1
		NetReaction rsI1 = new NetReaction(1);
		rsI1.addReaction(decayH1, 0.4f);
		inertH1.setProcesses(rsI1);
		//---InertH2
		NetReaction rsI2 = new NetReaction(1);
		rsI2.addReaction(decayH2, 0.4f);
		inertH2.setProcesses(rsI2);
		//---carbon source
		NetReaction rsAcetate = new NetReaction(4);
		rsAcetate.addReaction(sUptakeH1, -1);
		rsAcetate.addReaction(sUptakeH2, -1);
		rsAcetate.addReaction(decayH1, +0.6f);
		rsAcetate.addReaction(decayH2, +0.6f);
		carbonSource.setProcesses(rsAcetate);
		//---oxygen
		NetReaction rsOxygen = new NetReaction(5);
		rsOxygen.addReaction(sUptakeH1, -(1 - YSP));
		rsOxygen.addReaction(growthH1, -(YSP / YSX - 1));
		rsOxygen.addReaction(phbConsumption, -(1 - YPX));
		rsOxygen.addReaction(sUptakeH2, -(1 - YSP));
		rsOxygen.addReaction(growthH2, -(YSP / YSX - 1));
		oxygen.setProcesses(rsOxygen);
		//           ---- Add the species to system ----
		addBiomassSpecies(speciesH1);
		addBiomassSpecies(speciesH2);
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
			int[] nCells = {initialCellNumber, initialCellNumber};
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
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length == 1)
			outputDirectory = args[0];
		MultigridVariable.setSteps(10, 100);
		ApplicationComponent app = new PHBProducerEPSProducerSameUptake();
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
			app.addStateWriter(new PovRayWriter());
			app.addStateWriter(new SoluteConcentrationWriter());
			app.addStateWriter(new SolidsConcentrationWriter());
			app.addStateWriter(new ParticlePositionWriter());
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
		//TODO remove this line
		System.out.println("done.");
	}
}