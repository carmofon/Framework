package nl.tudelft.bt.model.tests.sbr;

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
import nl.tudelft.bt.model.detachment.levelset.functions.Height2VolumetricDetachment;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.reaction.*;

/**
 * Example for model description paper
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public class ConstantFeedingStratification extends ModelHandler {
	// geometry (2D or 3D)
	private static int geometry = 2;

	//Feast famine cycles:
	private static float feastTime = 0.0417f; // 1 h feast

	private static float famineTime = 0.0833f; // 2 h famine

	//Solute species
	//Oxygen
	private static float oxygenFeastConcentration = 0f; //[gO/m3]

	private static float oxygenFamineConcentration = 4f; //[gO/m3]

	private static float oxygenDiffusivity = 0.0002f; //[m2/d]

	private static float KO_PAO = 0.1f; //[gO/m3]

	private static float KO_NH = 0.5f; //[gO/m3]

	private static float KO_NO = 0.5f; //[gO/m3]

	//acetate
	private static float acetateFeastConcentration = 0.2f; //400f[g-CODAc/m3]

	private static float acetateFamineConcentration = 0f; //[g-CODAc/m3]

	private static float acetateDiffusivity = 0.000138f; //[m2/d]

	private static float KAc_PAO = 20f; //[g-CODAc/m3]

	//Ammonium
	private static float ammoniumBulkConcentration = 50f; //[g-N/m3]

	private static float ammoniumDiffusivity = 0.001f; //[m2/d]

	private static float KNH4_PAO = 0.01f; //[g-N/m3]

	private static float KNH4_NH = 2.4f; //[g-N/m3]

	//Nitrate
	private static float no3BulkConcentration = 0f; //[g-N/m3]

	private static float no3Diffusivity = 0.001f; //[m2/d]

	private static float KNO3_PAO = 0.5f; //[g-N/m3]

	private static float etaNO3 = 0.8f; //[g-N/m3]

	//Nitrite
	private static float no2BulkConcentration = 0f; //[g-N/m3]

	private static float no2Diffusivity = 0.001f; //[m2/d]

	private static float KNO2_PAO = 0.5f; //[g-N/m3]

	private static float KNO2_NO = 0.5f; //[g-N/m3]

	private static float etaNO2 = 0.8f; //[g-N/m3]

	//Phosphate
	private static float po4FeastConcentration = 15f; //[g-P/m3]

	private static float po4FamineConcentration = 0f; //[g-P/m3]

	private static float po4Diffusivity = 0.001f; //[m2/d]

	private static float KPO4_PAO = 0.5f; //[g-P/m3]

	private static float KPO4_PAO_PP = 3.1f; //[g-P/m3]

	private static float KPO4_NH = 0.01f; //[g-P/m3]

	private static float KPO4_NO = 0.01f; //[g-P/m3]

	private static float g_PP = 0.9f; //[]

	//Particulate species
	// biomass density
	private static float densityX = 142000f; //[gCOD/m3]

	//PHB
	private static float densityPHB = 1.42e7f; //[gCOD/m3]

	private static float kFPhb = 0.33f; //[g-PHB/gCOD-XPAO]

	//PP
	private static float densityPP = densityPHB; //[gCOD/m3]

	private static float KPP_PAO = 0.01f; //[g-P/m3]

	private static float kFPp = 0.01f; //[g-P/gCOD-XPAO]

	private static float fmaxPp = 0.35f; //[g-P/gCOD-XPAO]

	//
	//Yield coefficients
	private static float f_XI = 0.1f; //	gCOD/gCOD

	private static float i_NBM = 0.08f; //	gN/gCOD

	private static float i_NXI = 0.03f; //	gN/gCOD

	private static float i_NXS = 0.04f; //	gN/gCOD

	private static float i_PBM = 0.02f; //	gP/gCOD

	private static float i_PXI = 0f; //	gP/gCOD

	private static float i_PXS = 0f; //	gP/gCOD

	private static float mATP = 0.014f; //	molATP/gCOD/d

	private static float Y_AM = 0.159f; //	gCOD/gN

	private static float Y_GLY = 0.5f; //	gCOD/gCOD

	private static float Y_GLYNO2 = 1.15f; //	gCOD/gCOD

	private static float Y_GLYNO3 = 1.15f; //	gCOD/gCOD

	private static float Y_GLYO = 1.12f; //	gCOD/gCOD

	private static float Y_NH = 0.15f; //	gCOD/gN

	private static float Y_NO = 0.041f; //	gCOD/gN

	private static float Y_PHB = 1.5f; //	gCOD/gCOD

	private static float Y_PHBNO2 = 1.7f; //	gCOD/gCOD

	private static float Y_PHBNO3 = 1.7f; //	gCOD/gCOD

	private static float Y_PHBO = 1.39f; //	gCOD/gCOD

	private static float Y_PO4 = 0.3f; //	gP/gCOD

	private static float Y_PPNO2 = 2.92f; //	gP/gCOD

	private static float Y_PPNO3 = 2.92f; //	gP/gCOD

	private static float Y_PPO = 4.37f; //	gP/gCOD

	//
	// Processes
	// PAO
	// Anaerobic
	//Storage of PHB
	private static float qMaxAc = 10f; //[gCOD_Ac/gCOD/d]

	//[gCOD_Ac/gCOD/d]
	// Aerobic
	//Consumption of PHB
	private static float kPhb = 7.55f; //[gCOD_PHB/gCOD/d]

	//Storage of PP
	private static float kPp = 0.45f; //[g-P/gCOD/d]

	//Autotrophs growth
	//NH
	private static float umax_XNH = 2.05f; //[gCOD-XNH/gCOD-XNH/d]

	//NO
	private static float umax_XNO = 2f; //[gCOD-XNH/gCOD-XNH/d]

	//
	// Computation parameters
	private static float systemSize = 0.001f; // [m]

	private static float relativeMaximumRadius = 0.006f;

	private static float relativeMinimumRadius = relativeMaximumRadius * 0.0001f;

	private static float relativeBoundaryLayer = 0.1f;

	// other model parameters
	private static int gridSide = 17; // multigrid grid side

	private static float kShov = 1.0f; // shoving parameter[dim/less]

	private static float rdetach = 70.423f; //detachment coef for computed prof

	//private static float rdetach = 5e9f; //for imposed O2 profiles
	// detachment constant[g/l/h]
	private static int initialCellNumber = 800;

	//private static int initialCellNumber = 31;
	/**
	 * Define the single bacteria species, the chemical species and the
	 * processes
	 */
	private void defineSpeciesAndReactions() throws ModelException {
		// create the solutes
		//oxygen
		SoluteSpecies oxygen = new SoluteSpecies("oxygen", oxygenDiffusivity);
		oxygen.setBulkConcentration(new ConstantBulkConcentration(
				oxygenFamineConcentration));
		//ammonium
		SoluteSpecies ammonium = new SoluteSpecies("ammonium",
				ammoniumDiffusivity);
		ammonium.setBulkConcentration(new ConstantBulkConcentration(
				ammoniumBulkConcentration));
		//nitrite (NO2)
		SoluteSpecies no2 = new SoluteSpecies("no2", no2Diffusivity);
		no2.setBulkConcentration(new ConstantBulkConcentration(
				no2BulkConcentration));
		// create the particulates
		// PAO - phosphate accumulating organisms
		//PAO active mass
		ParticulateSpecies activePAO = new ParticulateSpecies("PAOActiveMass",
				densityX, Color.yellow);
		// array of fixed species that constitute PAO organisms
		ParticulateSpecies[] spPao = { activePAO };
		float[] fractionalVolumeCompositionPao = { 1 };
		BiomassSpecies pao = new BiomassSpecies("pao", spPao,
				fractionalVolumeCompositionPao);
		// XNH - ammonium oxydizers
		//active mass
		ParticulateSpecies activeXnh = new ParticulateSpecies("XNHActiveMass",
				densityX, Color.blue);
		// array of fixed species that constitute XNH organisms
		ParticulateSpecies[] spXnh = { activeXnh };
		float[] fractionalVolumeCompositionXnh = { 1f };
		// create the biomass species
		BiomassSpecies xnh = new BiomassSpecies("xnh", spXnh,
				fractionalVolumeCompositionXnh);
		//Create the Reaction factors, Monod and inhibition coefficients
		// PAO
		ProcessFactor iO = new Inhibition(oxygen, KO_PAO);
		ProcessFactor mO = new Saturation(oxygen, KO_PAO);
		ProcessFactor mNo2 = new Saturation(no2, KNO2_PAO);
		// XNH
		ProcessFactor mO_XNH = new Saturation(oxygen, KO_NH);
		ProcessFactor mNh4_XNH = new Saturation(ammonium, KNH4_NH);
		// create the reactions
		//PAO
		//aerobic
		//consumption of PHB (growth)
//		Reaction aerobicConsumptionOfPhb = new Reaction(
//				"aerobicConsumptionOfPhb", activePAO, kPhb * 0.05f, 1);
//		//NOTE: 0.5 accounts for growth limitting PHB fraction
//		aerobicConsumptionOfPhb.addFactor(mO);
		//anoxic (NO2)
		//consumption of PHB (growth)
		Reaction anoxicNo2ConsumptionOfPhb = new Reaction(
				"anoxicNo2ConsumptionOfPhb", activePAO, kPhb * etaNO2 * 0.0001f, 2);
		anoxicNo2ConsumptionOfPhb.addFactor(mNo2);
		anoxicNo2ConsumptionOfPhb.addFactor(iO);
		//XNH - ammonium oxidizer organisms
		//growth (nitrification)
		Reaction xnhGrowth = new Reaction("xnhGrowth", activeXnh, umax_XNH, 2);
		xnhGrowth.addFactor(mO_XNH);
		xnhGrowth.addFactor(mNh4_XNH);
		// assign reaction to the species through ReactionStoichiometries
		//Particulate species
		//PAO active mass
		NetReaction rsActivePao = new NetReaction(1);
//		rsActivePao.addReaction(aerobicConsumptionOfPhb, 1 / Y_PHBO);
		rsActivePao.addReaction(anoxicNo2ConsumptionOfPhb, 1 / Y_PHBNO2);
		activePAO.setProcesses(rsActivePao);
		// XNH active mass
		NetReaction rsXnh = new NetReaction(1);
		rsXnh.addReaction(xnhGrowth, 1);
		activeXnh.setProcesses(rsXnh);
		//Particulate species
		//Oxygen
		NetReaction rsO = new NetReaction(1);
//		rsO.addReaction(aerobicConsumptionOfPhb, 1 / Y_PHB - 1);
		rsO.addReaction(xnhGrowth, 1 - 3.43f / Y_NH);
		oxygen.setProcesses(rsO);
		//NH4
		NetReaction rsNh4 = new NetReaction(1);
		rsNh4.addReaction(xnhGrowth, -1 / Y_NH - i_NBM);
		ammonium.setProcesses(rsNh4);
		//NO2
		NetReaction rsNo2 = new NetReaction(2);
		rsNo2.addReaction(anoxicNo2ConsumptionOfPhb, (1 - Y_PHBNO2) / 1.71f
				* Y_PHBNO2);
		rsNo2.addReaction(xnhGrowth, 1 / Y_NH);
		no2.setProcesses(rsNo2);
		// add the species to system
		addBiomassSpecies(pao);
		addBiomassSpecies(xnh);
		addSoluteSpecies(oxygen);
		addSoluteSpecies(ammonium);
		addSoluteSpecies(no2);
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
			int[] nCells = { initialCellNumber, initialCellNumber };
			//			int[] nCells = { initialCellNumber, initialCellNumber,
			//					initialCellNumber };
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
		DetachmentSpeedFunction df = new Height2VolumetricDetachment(rdetach);
		setDetachmentHandler(df);
	}

	/**
	 * Simulation of 1000 iterative steps, storing results at each iteration
	 * 
	 * @param args
	 */
	public static void main(String[] args) {
		//output directory name
		//String outputDirectory = "D:/results/"
		//	+ "laspidouModel/2D/computedO2_"+ rdetach +"/";
		String outputDirectory = "E:/results/SBR/simple1/lixo3/";
		//
		MultigridVariable.setSteps(100, 1000);
		ApplicationComponent app = new ConstantFeedingStratification();
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
			spw.addSeries(Model.detachedBiomassContainer()
					.getTotalDetachedBiomassSeries());
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