package nl.tudelft.bt.model.apps.components;

import java.io.IOException;
import java.util.Iterator;

import nl.tudelft.bt.model.BiomassSpecies;
import nl.tudelft.bt.model.apps.ApplicationComponent;
import nl.tudelft.bt.model.apps.output.StateWriter;
import nl.tudelft.bt.model.apps.output.StateWriterInterface;
import nl.tudelft.bt.model.exceptions.*;
import nl.tudelft.bt.model.multigrid.SoluteSpecies;

/**
 * Application decorator abstract class, for implementation of decorator design
 * pattern
 * 
 * @author Joao Xavier (j.xavier@tnw.tudelft.nl)
 */
public abstract class ApplicationDecorator extends ApplicationComponent {
	protected ApplicationComponent component;

	/**
	 * Creates instance, getting component from previous initialization
	 */
	public ApplicationDecorator(ApplicationComponent c) {
		component = c;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#reset()
	 */
	public void reset() {
		component.reset();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#setSystemSpaceParameters(int,
	 *      float, float, float, int, float, float)
	 */
	public void setSystemSpaceParameters(int d, float s, float rmax,
			float rmin, float h, int g, float k) {
		component.setSystemSpaceParameters(d, s, rmax, rmin, h, g, k);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#initializeSystemSpace()
	 */
	public void initializeSystemSpace() throws InvalidValueException {
		component.initializeSystemSpace();
	}

	/*
	 * (non-Javadoc)
	 */
	public void intializeStateWriters(String dir) throws ModelIOException {
		component.intializeStateWriters(dir);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#addBacteriumSpecies(org.photobiofilms.phlip.multigrid.BacteriaSpecies)
	 */
	public void addBiomassSpecies(BiomassSpecies s)
			throws SystemEditViolationException {
		component.addBiomassSpecies(s);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#addSoluteSpecies(org.photobiofilms.phlip.multigrid.ChemicalSpecies)
	 */
	public void addSoluteSpecies(SoluteSpecies c)
			throws SystemEditViolationException {
		component.addSoluteSpecies(c);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#initializeDiffusionReactionSystem()
	 */
	public void initializeDiffusionReactionSystem() throws ModelException {
		component.initializeDiffusionReactionSystem();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#writeState()
	 */
	public void writeState() throws ModelException {
		component.writeState();
	}

	
	/**
	 * Force writing only the timed state writers.
	 */
	public void forceWriteTimedStateWriters() throws ModelException {
		component.forceWriteTimedStateWriters();
	}

	
	/**
	 * Force writing of all state writers to disk to use as last write
	 * 
	 * @throws ModelException
	 */
	public void forceWriteStateWithoutSerializing() throws ModelException {
		component.forceWriteStateWithoutSerializing();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#writeState()
	 */
	public void forceWriteState() throws IOException {
		component.forceWriteState();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#performGrowthAndDivision()
	 */
	public void performGrowthAndDivision() throws ModelException {
		component.performGrowthAndDivision();
		output();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#performShovingStep()
	 */
	public void performSpreadingStep() {
		component.performSpreadingStep();
		output();
	}

	/**
	 * Perform the biomass shrinking by pressure
	 */
	public void performSpreadingByPulling() {
		component.performSpreadingByPulling();
		output();
	}

	/**
	 * Performs shoving steps until full relaxation
	 */
	public void spreadByShovingCompletely() {
		component.spreadByShovingCompletely();
		output();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#relaxCompletelyByShoving()
	 */
	public void spreadCompletely() {
		component.spreadCompletely();
		output();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#detach()
	 */
	public void detach() throws ModelException {
		component.detach();
		output();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.phlip.apps.ApplicationComponent#computeSoluteConcentrations()
	 */
	public void computeSoluteConcentrations()
			throws MultigridSystemNotSetException {
		component.computeSoluteConcentrations();
	}

	/**
	 * Refresh component output, which can either be refreshing display for
	 * gui's or terminal text interfaces, or outputting internal state to disk
	 */
	protected abstract void output();

	/*
	 * (non-Javadoc)
	 */
	public String toString() {
		return component.toString();
	}

	/*
	 * (non-Javadoc)
	 */
	public void addStateWriter(StateWriterInterface w)
			throws SystemEditViolationException {
		component.addStateWriter(w);
	}

	/*
	 * (non-Javadoc)
	 */
	public void addTimedStateWriter(StateWriterInterface w) {
		component.addTimedStateWriter(w);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.ApplicationComponent#markForDetachment()
	 */
	public void markForDetachment() {
		component.markForDetachment();
		output();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.photobiofilms.model.apps.ApplicationComponent#initializeDetachmentFunction()
	 */
	public void initializeDetachmentFunction() {
		component.initializeDetachmentFunction();
		output();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see nl.tudelft.bt.model.apps.ApplicationComponent#setFinishIterationSwitch()
	 */
	public void setFinishSimulationSwitch() {
		component.setFinishSimulationSwitch();

	}
}