package nl.tudelft.bt.model.apps;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Iterator;
import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.output.StateWriter;
import nl.tudelft.bt.model.apps.output.StateWriterInterface;
import nl.tudelft.bt.model.detachment.DetachmentHandler;
import nl.tudelft.bt.model.detachment.levelset.functions.DetachmentSpeedFunction;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.*;
import nl.tudelft.bt.model.multigrid.boundary_conditions.BiofilmBoundaryConditions;
import nl.tudelft.bt.model.multigrid.boundary_layers.BoundaryLayer;
import nl.tudelft.bt.model.multigrid.boundary_layers.DilationFillHolesBoundaryLayer;
import nl.tudelft.bt.model.multigrid.boundary_layers.SphericalDilationBoundaryLayer;
import nl.tudelft.bt.model.particlebased.BiomassParticle;
import nl.tudelft.bt.model.particlebased.BiomassParticleContainer;
import nl.tudelft.bt.model.povray.Povray3DScene;

/**
 * Application interface - all applications will implement this interface
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class ModelHandler extends ApplicationComponent {
	protected Model _m;

	private boolean _spaceSystemIsEditable;

	protected boolean _physicalSystemIsEditable;

	// system size parameters
	protected int _geometry; // 2D or 3D

	protected float _systemSize;

	protected float _relativeMaximumRadius;

	protected float _relativeMinimumRadius;

	private float _relativeBoundaryLayerHeight;

	protected int gridSize; // for solutes concentration

	protected float _kShoving;

	protected static final float FSHOVING = 0.01f;

	// containers
	ArrayList _biomassSpecies;

	ArrayList _soluteSpecies;

	// container of state output objects
	ArrayList _stateWritters;

	ArrayList _timedStateWritters;

	protected BoundaryLayer _boundaryLayer;

	/**
	 * Constructor
	 */
	public ModelHandler() {
		_m = Model.model();
		reset();
	}

	/**
	 * Resets all model properties
	 */
	public void reset() {
		_m.reset();
		_biomassSpecies = new ArrayList();
		_soluteSpecies = new ArrayList();
		_stateWritters = new ArrayList();
		_timedStateWritters = new ArrayList();
		_spaceSystemIsEditable = true;
		_physicalSystemIsEditable = true;
	}

	/**
	 * Sets the system space parameters
	 * 
	 * @param d
	 *            system dimensionallity (2 or 3)
	 * @param s
	 *            system size
	 * @param rmax
	 *            relative maximum biomass particle radius
	 * @param rmin
	 *            relative minimum biomass particle radius
	 * @param h
	 *            relative height of boundary layer
	 * @param g
	 *            grid size for solute concentration
	 * @param k
	 *            shoving multiplier
	 * @param rdetach
	 *            detachment rate [10^-15g/um^2/h]
	 */
	public void setSystemSpaceParameters(int d, float s, float rmax,
			float rmin, float h, int g, float k) {
		// check if system can be edited
		if (!_spaceSystemIsEditable) {
			throw (new SystemEditViolationException("Illegal setting "
					+ "system parameters"));
		}
		// validate the grid resolution
		float gridElementSize = s / g;
		float maximumDiameter = rmax * s * 2;
		if (rmax > 1 / (g * 2.5f)) {
			throw new ModelRuntimeException("Grid resolution (" + g
					+ ") to high" + " for particle radius used\n"
					+ "maximum value for the maximum relative radius is " + 1
					/ (g * 2.5f));
		}
		// set the attributes
		_geometry = d;
		_systemSize = s;
		_relativeMaximumRadius = rmax;
		_relativeMinimumRadius = rmin;
		_relativeBoundaryLayerHeight = h;
		gridSize = g;
		_kShoving = k;
	}

	/*
	 * (non-Javadoc)
	 */
	public void initializeSystemSpace() throws InvalidValueException {
		_m.buildSystem(_geometry, _systemSize, gridSize);
		_spaceSystemIsEditable = false;
	}

	/**
	 * Initialize the system writters and define the output directory. NOTE:
	 * this operation is for all state writers
	 * 
	 * @param dir
	 *            directory to write to
	 */
	public void intializeStateWriters(String dir) {
		StateWriter.setModelHandle(this);
		try {
			StateWriter.setNewResultDirectoryOrDeleteContents(dir);
		} catch (IOException e) {
			throw new ModelIOException("Not able to create output directory");
		}
	}

	/**
	 * Restart the state writer for output. Set the value of the ouput directory
	 * to an existing directory and the handle to the model handler
	 * 
	 * @param dir
	 *            directory name
	 */
	public void restartStateWriters(String dir) {
		StateWriter.setModelHandle(this);
		StateWriter.reUseResultDirectory(dir);
	}

	/**
	 * Add a new bacterium species to system. Definition of processes associated
	 * to this species are left to the user classes.
	 * 
	 * @param s
	 *            bacterium species to add
	 * @throws SystemNotEditableException
	 *             if system is not in editable mode
	 */
	public void addBiomassSpecies(BiomassSpecies s)
			throws SystemEditViolationException {
		if (!_physicalSystemIsEditable) {
			throw (new SystemEditViolationException("Illegal addition of"
					+ " bacterium species"));
		}
		_biomassSpecies.add(s);
	}

	/**
	 * Add a new solute species to system. Definition of processes associated to
	 * this species are left to the user classes.
	 * 
	 * @param c
	 *            solute species to add
	 * @throws SystemNotEditableException
	 *             if system is not in editable mode
	 */
	public void addSoluteSpecies(SoluteSpecies c)
			throws SystemEditViolationException {
		if (!_physicalSystemIsEditable)
			throw (new SystemEditViolationException("Illegal addition of"
					+ " solute species"));
		_soluteSpecies.add(c);
	}

	/*
	 * (non-Javadoc)
	 */
	public void addStateWriter(StateWriterInterface w) {
		_stateWritters.add(w);
	}

	/*
	 * (non-Javadoc)
	 */
	public void addTimedStateWriter(StateWriterInterface w) {
		_timedStateWritters.add(w);
	}

	/**
	 * Initialize system. Also performs innoculation, using abstract method
	 * innoculate
	 * 
	 * @throws ModelException
	 */
	public void initializeDiffusionReactionSystem() throws ModelException {
		if (_spaceSystemIsEditable)
			throw new SystemEditViolationException(
					"Trying to set physical system" + " but space is not set");
		// convert the collections to arrays
		BiomassSpecies[] sp = new BiomassSpecies[_biomassSpecies.size()];
		int i = 0;
		for (Iterator iter = _biomassSpecies.iterator(); iter.hasNext();) {
			BiomassSpecies element = (BiomassSpecies) iter.next();
			sp[i++] = element;
		}
		SoluteSpecies[] cs = new SoluteSpecies[_soluteSpecies.size()];
		i = 0;
		for (Iterator iter = _soluteSpecies.iterator(); iter.hasNext();) {
			SoluteSpecies element = (SoluteSpecies) iter.next();
			cs[i++] = element;
		}
		// set the maximum biofilm height to leave enough free space
		// to boundary layer
		_m.setVerticalCutoffSize(_systemSize
				* (1 - _relativeBoundaryLayerHeight) * 0.9f);
		try {
			// create the boundary layer
			createBoundaryLayer(_systemSize * _relativeBoundaryLayerHeight);
			// create the
			_m.setupDiffusionReactionSystem(sp, cs, _boundaryLayer);
			// initialize the bacteria container
			builBiomassParticleContainer();
			// inoculate
			inoculate();
			// compute initial concentration for chemicals
			initializeSoluteConcentrations();
			// close system edit mode
			_physicalSystemIsEditable = false;
		} catch (MultigridSystemNotSetException e) {
			throw (new SystemEditViolationException(e.toString()));
		}
	}

	protected void createBoundaryLayer(float h)
			throws MultigridSystemNotSetException {
		// uncomment the type of boundary layer wanted
		// _boundaryLayer = new FrontBoundaryLayer(h);
		_boundaryLayer = new SphericalDilationBoundaryLayer(h);
		// create the boundary conditions
		MultigridVariable.setBoundaryConditions(new BiofilmBoundaryConditions());
	}

	protected void builBiomassParticleContainer() {
		_m.buildBacteriaContainer(_relativeMaximumRadius * _systemSize,
				_relativeMinimumRadius * _systemSize, _kShoving, FSHOVING);
	}

	/**
	 * Add randomly distributed particles. The number of particles of each
	 * species to inoculate the system is defined from array nCells, which must
	 * have the same length as the number of bacterium species in the system
	 * 
	 * @param nCells
	 *            number of cells of each species
	 * @throws NonMatchingNumberException
	 *             if length of nCells array does not match the number of
	 *             species in the system
	 * @throws SystemNotEditableException
	 *             if system is not in editable mode
	 */
	protected void inoculateRandomly(int[] nCells)
			throws NonMatchingNumberException, SystemEditViolationException {
		if (!_physicalSystemIsEditable) {
			throw (new SystemEditViolationException(
					"Tried to innoculate system" + " when not in edit mode"));
		}
		_m.inoculateRandomlyMultispecies(nCells);
	}

	/**
	 * Place a simgle cell of species s in center of computational volume
	 * 
	 * @param s
	 */
	protected void placeSingleCellInCenter(BiomassSpecies s) {
		_m.placeSingleCellInCenter(s);		
	}
	
	/**
	 * Place a quantity of bimass in a given postision
	 * 
	 * @param s
	 * @param x
	 * @param y
	 * @param z
	 * @return reference to particle added
	 */
	public BiomassParticle placeBiomass(BiomassSpecies s, float x, float y,
			float z) {
		return _m.placeBiomass(s, x, y, z);
	}

	/**
	 * Inoculates the system with biomass
	 */
	protected abstract void inoculate();

	/*
	 * (non-Javadoc)
	 */
	public void writeState() throws ModelException {
		// write state of all state writers
		if (!(_stateWritters == null))
			for (Iterator iter = _stateWritters.iterator(); iter.hasNext();) {
				StateWriterInterface element = (StateWriterInterface) iter
						.next();
				element.write();
			}
		// write state of all timed state writers
		if (!(_timedStateWritters == null) & (_m.writeTimedWriters()))
			for (Iterator iter = _timedStateWritters.iterator(); iter.hasNext();) {
				StateWriterInterface element = (StateWriterInterface) iter
						.next();
				element.write();
			}
	}

	/**
	 * Force writing only the timed state writers.
	 */
	public void forceWriteTimedStateWriters() throws ModelException {
		// write state of all state writers
		if (!(_timedStateWritters == null))
			for (Iterator iter = _timedStateWritters.iterator(); iter.hasNext();) {
				StateWriterInterface element = (StateWriterInterface) iter
						.next();
				element.dump();
			}
	}
	
	/**
	 * Force writing of all state writers to disk for debbuging puproses. Also
	 * writes the object to disk a a serialized object
	 */
	public void forceWriteStateWithoutSerializing() throws ModelException {
		// write state of all state writers
		if (!(_stateWritters == null))
			for (Iterator iter = _stateWritters.iterator(); iter.hasNext();) {
				StateWriterInterface element = (StateWriterInterface) iter
						.next();
				element.dump();
			}
	}

	
	/**
	 * Force writing of all state writers to disk for debbuging puproses. Also
	 * writes the object to disk a a serialized object
	 */
	public void forceWriteState() throws IOException {
		// write state of all state writers
		if (!(_stateWritters == null))
			for (Iterator iter = _stateWritters.iterator(); iter.hasNext();) {
				StateWriterInterface element = (StateWriterInterface) iter
						.next();
				element.dump();
			}
		// write the object to file
		File dir = new File(StateWriter.getOutputDirectoryName());
		ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(
				StateWriter.getOutputDirectoryName()
						+ "/modelhandler_iteration"
						+ Model.model().getFormatedIterationNumber() + ".dump"));
		out.writeObject(this);
		MultigridVariable.serializeStaticState(out);
		StateWriter.serializeStaticState(out);
		Povray3DScene.serializeStaticState(out);
		out.close(); // Also flushes output
	}

	/**
	 * Load a model handler from file
	 * 
	 * @param fileToRead
	 * @return th einstance of ModelHandler read from file
	 * @throws IOException
	 */
	public static ModelHandler loadFromDumpFile(String fileToRead)
			throws IOException, ClassNotFoundException {
		ObjectInputStream in = new ObjectInputStream(new FileInputStream(
				fileToRead));
		ModelHandler app = (ModelHandler) in.readObject();
		MultigridVariable.deserializeStaticState(in);
		StateWriter.deserializeStaticState(in);
		Povray3DScene.deserializeStaticState(in);
		return app;
	}

	/**
	 * Biomass growth and division
	 */
	public void performGrowthAndDivision() throws ModelException {
		_m.performGrowthAndDivision();
	}

	/**
	 * Single spreading step
	 */
	public void performSpreadingStep() {
		_m.performSpreadingStep();
	}

	/**
	 * Perform the biomass shrinking by pressure
	 */
	public void performSpreadingByPulling() {
		_m.performSpreadingByPulling();
	}

	/**
	 * Performs shoving steps until full relaxation
	 */
	public void spreadByShovingCompletely() {
		_m.spreadByShovingCompletely();
	}

	/**
	 * Performs shoving steps until full relaxation
	 */
	public void spreadCompletely() {
		_m.spread();
	}

	/**
	 * Detach biomass
	 */
	public void detach() throws ModelException {
		_m.detach();
	}

	/**
	 * Compute the solute concentration fields by solving the diffusion reaction
	 * PDE's
	 * 
	 * @throws MultigridSystemNotSetException
	 */
	public void computeSoluteConcentrations()
			throws MultigridSystemNotSetException {
		_m.updateSoluteConcentrations();
	}

	/**
	 * Initialize solute concentration fields by solving the diffusion reaction
	 * PDE's
	 * 
	 * @throws MultigridSystemNotSetException
	 */
	public void initializeSoluteConcentrations()
			throws MultigridSystemNotSetException {
		_m.initializeSoluteConcentrations();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		String s = "" + _m.getIterationCounter() + " - t = " + _m.getTime()
				+ " h; timestep - " + _m.getTimeStep() + " h; particles - ";
		for (Iterator iter = _biomassSpecies.iterator(); iter.hasNext();) {
			BiomassSpecies element = (BiomassSpecies) iter.next();
			s += element.getName() + ":" + _m.getNumberOfParticles(element)
					+ " time: " + System.currentTimeMillis() + "  ";
		}
		return s;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.ApplicationComponent#markForDetachment()
	 */
	public void markForDetachment() {
		((BiomassParticleContainer) Model.model().biomassContainer)
				.markParticlesForErosion();
	}

	/**
	 * Sets the detachment finction to use
	 * 
	 * @param df
	 */
	public void setDetachmentHandler(DetachmentHandler df) {
		((BiomassParticleContainer) Model.model().biomassContainer)
				.setDetachmentFunction(df);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.apps.ApplicationComponent#setFinishIterationSwitch()
	 */
	public void setFinishSimulationSwitch() {
		_m.setFinishSimulation();
	}
}