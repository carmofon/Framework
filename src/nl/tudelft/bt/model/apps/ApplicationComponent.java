package nl.tudelft.bt.model.apps;

import java.io.IOException;
import java.io.Serializable;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.Model;
import nl.tudelft.bt.model.apps.output.StateWriterInterface;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

/**
 * Interface that allows application components to be added to a guiven
 * application following the decorator design pattern
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class ApplicationComponent implements Serializable{
	private static boolean _simulationRunning;

	/**
	 * Resets all model properties
	 */
	public abstract void reset();

	/**
	 * Sets the system space parameters
	 * 
	 * @param d
	 *            system dimensionallity (2 or 3)
	 * @param s
	 *            system size
	 * @param r
	 *            relative maximum biomass particle radius
	 * @param h
	 *            relative height of boundary layer
	 * @param g
	 *            grid size for solute concentration
	 * @param k
	 *            shoving multiplier
	 * @param rdetach
	 *            detachment rate [10^-15g/um^2/h]
	 */
	public abstract void setSystemSpaceParameters(int d, float s, float rmax,
			float rmin, float h, int g, float k);

	/**
	 * Set model parameters relative to the reactor, for global mass balances
	 * 
	 * @param residenceTime
	 *            [h]
	 * @param carrierArea
	 *            [um^2]
	 * @param reactorVolume
	 *            [um^3]
	 */
	public void setReactorParameters(float residenceTime, float carrierArea,
			float reactorVolume) {
		Model.model().setResidenceTime(residenceTime);
		Model.model().setCarrierArea(carrierArea);
		Model.model().setReactorVolume(reactorVolume);
	}

	/**
	 * Initialize system
	 * 
	 * @throws InvalidValueException
	 */
	public abstract void initializeSystemSpace() throws InvalidValueException;

	/**
	 * Initialize the system writters and define the output directory
	 * 
	 * @param dir
	 *            output directory
	 */
	public abstract void intializeStateWriters(String dir);

	/**
	 * Add a new biomass species to system. Definition of processes associated
	 * to this species are left to the user classes.
	 * 
	 * @param s
	 *            biomass species to add
	 * @throws SystemNotEditableException
	 *             if system is not in editable mode
	 */
	public abstract void addBiomassSpecies(BiomassSpecies s)
			throws SystemEditViolationException;

	/**
	 * Add a new solute species to system. Definition of processes associated to
	 * this species are left to the user classes.
	 * 
	 * @param c
	 *            solute species to add
	 * @throws SystemNotEditableException
	 *             if system is not in editable mode
	 */
	public abstract void addSoluteSpecies(SoluteSpecies c)
			throws SystemEditViolationException;

	/**
	 * Adds a StateWriter that writes at each iteration of the simultaion cycle
	 * 
	 * @param w
	 *            system writter
	 * @throws SystemEditViolationException
	 */
	public abstract void addStateWriter(StateWriterInterface w)
			throws SystemEditViolationException;

	/**
	 * Adds a StateWriter that writes only at times defined by the
	 * _compulsoryTime atribute of Model
	 * 
	 * @param w
	 *            system writter
	 */
	public abstract void addTimedStateWriter(StateWriterInterface w);

	/**
	 * Initialize system. Also performs innoculation, using abstract method
	 * innoculate
	 * @throws ModelException If system may not be edited
	 */
	public abstract void initializeDiffusionReactionSystem() throws ModelException;

	public abstract void initializeDetachmentFunction();

	/**
	 * Wait for a notify all
	 * 
	 * @throws InterruptedException
	 */
	public synchronized void waitForStartIteratingRequest()
			throws InterruptedException {
		// stop here, and wait to be awaken by
		// other thread
		wait();
		_simulationRunning = false;
	}

	/**
	 * Start performing model iterations
	 */
	public synchronized void startIterating() throws InterruptedException {
		try {			
			// only start iterating if its not iterating already!
			if (!_simulationRunning) {
				_simulationRunning = true;
				//write iteration 0
				writeState();
				// simulation may be interrupted by other classes
				// running in other threads, by evoking interruptIteration
				while (!Model.model().endSimulation()) {
					if (!_simulationRunning) {
						wait();
					}
					performFullIteration();
					writeState();
					System.out.println(toString());
				}
				_simulationRunning = false;
			}
		} catch (ModelException e) {
			e.printStackTrace();
			System.out.println(e);
		}
	}

	/**
	 * Stop performing iterations
	 */
	public void stopIterating() {
		_simulationRunning = false;
	}

	/**
	 * continue performing iterations
	 */
	public synchronized void continueIterating() {
		_simulationRunning = true;
		notifyAll();
	}

	/**
	 * Requests a full iteration step (growth, division, full shoving
	 * relaxation, biomass detachment and psudo steady state determination) of
	 * solute concentration fields from the model
	 * 
	 * @throws ModelException
	 */
	public void performFullIteration() throws ModelException {
		performGrowthAndDivision();
		spreadCompletely();
		detach();
		computeSoluteConcentrations();
	}

	public abstract void setFinishSimulationSwitch();

	/**
	 * Performs all state writing operations (writing to disk)
	 * 
	 * @throws ModelException
	 */
	public abstract void writeState() throws ModelException;
	
	
	/**
	 * Force writing only the timed state writers.
	 */
	public abstract void forceWriteTimedStateWriters() throws ModelException;


	/**
	 * Force writing of all state writers to disk to use as last write
	 * 
	 * @throws ModelException
	 */
	public abstract void forceWriteStateWithoutSerializing() throws ModelException;

	
	/**
	 * Force writing of all state writers to disk for debbuging puproses
	 * @throws IOException TODO
	 */
	public abstract void forceWriteState() throws IOException;

	/**
	 * Biomass growth and division
	 */
	public abstract void performGrowthAndDivision() throws ModelException;

	/**
	 * Single shoving step
	 */
	public abstract void performSpreadingStep();

	/**
	 * Perform the biomass shrinking by pressure
	 */
	public abstract void performSpreadingByPulling();

	/**
	 * Performs shoving steps until full relaxation
	 */
	public abstract void spreadByShovingCompletely();

	/**
	 * Performs shoving steps until full relaxation
	 */
	public abstract void spreadCompletely();

	/**
	 * Detach biomass
	 */
	public abstract void detach() throws ModelException;

	public abstract void markForDetachment();

	/**
	 * Compute the solute concentration fields by solving the diffusion reaction
	 * PDE's
	 * 
	 * @throws MultigridSystemNotSetException
	 */
	public abstract void computeSoluteConcentrations()
			throws MultigridSystemNotSetException;

}